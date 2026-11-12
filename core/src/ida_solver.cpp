/// @file ida_solver.cpp
/// @brief SUNDIALS IDA 自适应 DAE 求解器实现
#include <fluxcell/core/solver/ida_solver.hpp>

#if defined(FLUXCELL_WITH_SUNDIALS)

#include <ida/ida.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#if defined(FLUXCELL_WITH_SUNDIALS_KLU)
#include <sunlinsol/sunlinsol_klu.h>
#include <sunmatrix/sunmatrix_sparse.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace fluxcell::solver {

namespace {

/// 超过该规模的系统改用 GMRES 迭代求解器（dense LU 成本 ~ N³）
constexpr long kDenseSizeLimit = 300;

void check(SUNErrCode code, const char* what) {
    if (code != SUN_SUCCESS) {
        throw std::runtime_error(std::string("IdaSolver: ") + what
                                 + " 失败 (SUNErrCode "
                                 + std::to_string(static_cast<int>(code)) + ")");
    }
}

} // namespace

IdaSolver::IdaSolver(DaeProblem& problem, real_t t0, real_t rel_tol)
    : problem_(problem) {
    const int n = problem_.dae_size();
    if (n <= 0) {
        throw std::invalid_argument("IdaSolver: 问题维度必须为正");
    }

    check(SUNContext_Create(SUN_COMM_NULL, &ctx_), "SUNContext_Create");

    // N_Vector 直接包裹 Eigen 缓冲（零拷贝交换）
    y_eigen_.setZero(n);
    yp_eigen_.setZero(n);
    res_eigen_.setZero(n);
    y_ = N_VMake_Serial(n, y_eigen_.data(), ctx_);
    yp_ = N_VMake_Serial(n, yp_eigen_.data(), ctx_);
    if (y_ == nullptr || yp_ == nullptr) {
        throw std::runtime_error("IdaSolver: N_VMake_Serial 失败");
    }

    // 变量类型掩码与绝对容限（缓冲存为成员：N_VMake 引用其内存）
    mask_eigen_.resize(n);
    problem_.dae_mask(mask_eigen_);
    id_ = N_VMake_Serial(n, mask_eigen_.data(), ctx_);
    has_algebraic_ = (mask_eigen_.array() == 0.0).any();

    atol_eigen_.resize(n);
    problem_.dae_atol(atol_eigen_);
    atol_ = N_VMake_Serial(n, atol_eigen_.data(), ctx_);

    // 初始状态打包
    problem_.dae_pack_state(y_eigen_, yp_eigen_);

    ida_ = IDACreate(ctx_);
    if (ida_ == nullptr) {
        throw std::runtime_error("IdaSolver: IDACreate 失败");
    }

    if (IDAInit(ida_, &IdaSolver::residual_callback, t0, y_, yp_) != IDA_SUCCESS) {
        throw std::runtime_error("IdaSolver: IDAInit 失败");
    }
    if (IDASVtolerances(ida_, rel_tol, atol_) != IDA_SUCCESS) {
        throw std::runtime_error("IdaSolver: IDASVtolerances 失败");
    }
    if (IDASetUserData(ida_, this) != IDA_SUCCESS) {
        throw std::runtime_error("IdaSolver: IDASetUserData 失败");
    }
    // 微分/代数变量区分：代数变量不参与局部误差控制（其导数无意义）
    if (has_algebraic_) {
        if (IDASetId(ida_, id_) != IDA_SUCCESS
            || IDASetSuppressAlg(ida_, SUNTRUE) != IDA_SUCCESS) {
            throw std::runtime_error("IdaSolver: IDASetId/SuppressAlg 失败");
        }
    }

    // 线性求解器：小系统 dense LU（稳健）；大系统优先 KLU 稀疏直解
    // （PyBaMM IDAKLUSolver 同构：稀疏 Jacobian 差分填充 + KLU 直接分解，
    // 根治无预条件 GMRES 在厚膜长循环中的收敛崩溃）；KLU 不可用时
    // 回退 GMRES（无预条件，雅可比-向量积由 IDA 内部差分近似）
    if (n <= kDenseSizeLimit) {
        SUNMatrix A = SUNDenseMatrix(n, n, ctx_);
        SUNLinearSolver LS = SUNLinSol_Dense(y_, A, ctx_);
        if (A == nullptr || LS == nullptr
            || IDASetLinearSolver(ida_, LS, A) != IDA_SUCCESS) {
            throw std::runtime_error("IdaSolver: dense 线性求解器初始化失败");
        }
    }
#if defined(FLUXCELL_WITH_SUNDIALS_KLU)
    else if (setup_klu(t0)) {
        sparse_active_ = true;
    }
#endif
    else {
        SUNLinearSolver LS =
            SUNLinSol_SPGMR(y_, SUN_PREC_NONE, 40, ctx_);
        if (LS == nullptr || IDASetLinearSolver(ida_, LS, nullptr) != IDA_SUCCESS) {
            throw std::runtime_error("IdaSolver: GMRES 线性求解器初始化失败");
        }
    }

    // 输出点推进可能需要远超默认上限（500）的内部步数
    IDASetMaxNumSteps(ida_, 1000000);

    // 初始一致化迭代额度放宽（默认 10 次迭代 / 4 次雅可比）：
    // 大电流段切换时 Newton 初值离解更远，需要更多余量
    IDASetMaxNumItersIC(ida_, 20);
    IDASetMaxNumJacsIC(ida_, 8);
}

#if defined(FLUXCELL_WITH_SUNDIALS_KLU)

void IdaSolver::cleanup_klu() {
    if (LS_klu_ != nullptr) {
        SUNLinSolFree(LS_klu_);
        LS_klu_ = nullptr;
    }
    if (A_sp_ != nullptr) {
        SUNMatDestroy(A_sp_);
        A_sp_ = nullptr;
    }
    sparse_active_ = false;
}

bool IdaSolver::setup_klu(real_t t0) {
    const int n = problem_.dae_size();
    const auto t_start = std::chrono::steady_clock::now();

    // ── 1. 逐列扰动探测稀疏模式 ──────────────────────────────────
    // J = ∂F/∂y + cj·∂F/∂y'：分两遍探测（仅扰 y_j / 仅扰 yp_j）取
    // 支持集并集，避免单遍 cj=1 下对角项 M_jj − ∂f_j/∂y_j 抵消漏检。
    // 行尺度（区分真导数响应与浮点噪声）：微分行 F = y' − f(y) 的
    // 评估量级 ≈ |yp_i|（pack 后 y' = f(y)），代数行 ≈ |g(y)| = |F_i|，
    // 统一取两者之和；阈值 = 1e-10 × 行尺度（差分噪声 ~ eps/σ ≈
    // 1e-12 × 行尺度，留 100 倍余量；漏检仅当 |J_ij| < 1e-6 × 行尺度/σ，
    // 对 Newton 无影响）
    Vector y = y_eigen_;
    Vector yp = yp_eigen_;
    Vector base(n), probe(n), resp(n);
    problem_.dae_residual(t0, y, yp, base);

    // 行尺度阈值仅依赖行号，预计算（阈值 = 1e-10 × 行尺度）
    std::vector<real_t> row_tol(n);
    for (int i = 0; i < n; ++i) {
        const real_t row_scale = std::abs(yp_eigen_(i)) + std::abs(base(i));
        row_tol[i] = row_scale > 0.0 ? 1e-10 * row_scale : 1e-12;
    }

    // 单遍逐列扰动探测：响应对（行, 阈值, 支持集容器）由调用方给定，
    // 与 base 差分后按行尺度阈值判定非零
    auto probe_pass = [&](std::vector<std::vector<int>>& out_cols) {
        for (int j = 0; j < n; ++j) {
            // σ 为相对步长 + 逐变量 atol 下限：绝对 σ 对小量纲变量是相对
            // 巨扰动（L_sei~1.5e-6 m 被 σ=1e-4 扰 67 倍 → 膜压降 0.45→30 V
            // 爆进 BV 指数 → 差分商 ~1e261 → klu_factor 崩溃）；atol 即
            // 变量零点尺度（与 IDA 自身差分格式的 1/ewt 下限同思路）
            const real_t scale = std::max(std::abs(y_eigen_(j)), atol_eigen_(j));
            const real_t sigma = 1e-4 * (scale > 0.0 ? scale : real_t(1.0));
            // ∂F/∂y 支持集：仅扰动 y_j
            y(j) = y_eigen_(j) + sigma;
            problem_.dae_residual(t0, y, yp, probe);
            y(j) = y_eigen_(j);
            resp = probe - base;
            // ∂F/∂y' 支持集（质量矩阵结构）：仅扰动 yp_j
            yp(j) = yp_eigen_(j) + sigma;
            problem_.dae_residual(t0, y, yp, probe);
            yp(j) = yp_eigen_(j);
            probe -= base;
            for (int i = 0; i < n; ++i) {
                if (std::abs(resp(i)) > row_tol[i]
                    || std::abs(probe(i)) > row_tol[i]) {
                    out_cols[j].push_back(i);
                }
            }
        }
    };

    std::vector<std::vector<int>> cols(n);
    probe_pass(cols);

    // ── 1b. 增广外电路端口的 CV 行支持集 ──
    // 稀疏模式在构造时固定，而外电路行形态随控制模式切换：CC 行
    // I_cell − I_target 只依赖 I_cell 列，CV 行 V(y) − v_target 还依赖
    // i_surf/δ/T/c_e。若只按 CC 探测，CV 段 reinit 后外电路行的这些
    // 列条目缺失 → Jacobian 缺块 → Newton 方向错误 → corrector 反复
    // 失败（IDA_CONV_FAIL，实测 4.4V CCCV 在 CV 首步即崩）。切换到
    // Voltage 再探测一遍取并集；v_target 不进 Jacobian，取任意值。
    if (problem_.dae_port_augmented()) {
        problem_.dae_set_voltage(0.0);
        problem_.dae_pack_state(y_eigen_, yp_eigen_);
        y = y_eigen_;
        yp = yp_eigen_;
        problem_.dae_residual(t0, y, yp, base);
        for (int i = 0; i < n; ++i) {
            const real_t row_scale = std::abs(yp_eigen_(i)) + std::abs(base(i));
            row_tol[i] = row_scale > 0.0 ? 1e-10 * row_scale : 1e-12;
        }
        std::vector<std::vector<int>> cols_cv(n);
        probe_pass(cols_cv);
        for (int j = 0; j < n; ++j) {
            cols[j].insert(cols[j].end(), cols_cv[j].begin(),
                           cols_cv[j].end());
        }
        // 恢复恒流模式与构造时的自洽 pack（供后续 IDAInit）
        problem_.dae_set_current(problem_.dae_current_target());
        problem_.dae_pack_state(y_eigen_, yp_eigen_);
    }

    // ── 1c. 列支持集裁剪（模型可选）：移除经中间变量链差分出的
    // 病态巨元（如 SEI 厚度列的 j_sei ∝ 1/L 链），避免行对角劣势
    // 拖垮 KLU pivot；被删反馈由残差内不动点/高阶小量保证可删
    for (int j = 0; j < n; ++j) {
        problem_.dae_jac_col_prune(j, cols[j]);
    }

    // ── 2. 距离-1 贪婪列着色 ──────────────────────────────────
    // 同色列不得共享任何残差行，该行差分响应才能唯一归因到列
    std::vector<std::vector<int>> rows_of(n);
    for (int j = 0; j < n; ++j) {
        for (int i : cols[j]) {
            rows_of[i].push_back(j);
        }
    }
    sp_color_.assign(n, -1);
    int ncolors = 0;
    std::vector<char> used;
    for (int j = 0; j < n; ++j) {
        std::fill(used.begin(), used.end(), 0);
        for (int i : cols[j]) {
            for (int j2 : rows_of[i]) {
                const int c = sp_color_[j2];
                if (c >= 0 && c < static_cast<int>(used.size())) {
                    used[c] = 1;
                }
            }
        }
        int c = 0;
        while (c < static_cast<int>(used.size()) && used[c] != 0) {
            ++c;
        }
        if (c == static_cast<int>(used.size())) {
            used.push_back(0);
        }
        used[c] = 1;
        sp_color_[j] = c;
        ncolors = std::max(ncolors, c + 1);
    }

    // ── 3. CSC 模式（列内行升序去重）+ 槽位映射 + 颜色分组 ─────
    // CC/CV 两遍探测取并集后同一 (行, 列) 可能重复入列，去重后再建 CSC
    sp_rows_of_col_ = std::move(cols);
    for (auto& r : sp_rows_of_col_) {
        std::sort(r.begin(), r.end());
        r.erase(std::unique(r.begin(), r.end()), r.end());
    }

    sp_colptr_.assign(n + 1, 0);
    for (int j = 0; j < n; ++j) {
        sp_colptr_[j + 1] = sp_colptr_[j]
                          + static_cast<sunindextype>(sp_rows_of_col_[j].size());
    }
    const sunindextype nnz = sp_colptr_[n];
    sp_rowidx_.clear();
    sp_rowidx_.reserve(nnz);
    for (int j = 0; j < n; ++j) {
        for (int i : sp_rows_of_col_[j]) {
            sp_rowidx_.push_back(i);
        }
    }
    sp_slots_.assign(n, {});
    for (int j = 0; j < n; ++j) {
        sp_slots_[j].resize(sp_rows_of_col_[j].size());
        for (std::size_t k = 0; k < sp_slots_[j].size(); ++k) {
            sp_slots_[j][k] = static_cast<int>(sp_colptr_[j]
                                              + static_cast<sunindextype>(k));
        }
    }
    sp_groups_.assign(ncolors, {});
    for (int j = 0; j < n; ++j) {
        sp_groups_[sp_color_[j]].push_back(j);
    }
    jac_y_buf_.setZero(n);
    jac_yp_buf_.setZero(n);
    jac_f1_buf_.setZero(n);
    jac_sigma_buf_.setZero(n);

    // ── 4. SUNDIALS KLU 对象（全部成功后才挂接 IDA）────────────
    A_sp_ = SUNSparseMatrix(n, n, nnz, CSC_MAT, ctx_);
    if (A_sp_ != nullptr) {
        std::copy(sp_colptr_.begin(), sp_colptr_.end(),
                  SUNSparseMatrix_IndexPointers(A_sp_));
        std::copy(sp_rowidx_.begin(), sp_rowidx_.end(),
                  SUNSparseMatrix_IndexValues(A_sp_));
        LS_klu_ = SUNLinSol_KLU(y_, A_sp_, ctx_);
    }
    if (LS_klu_ == nullptr
        || IDASetLinearSolver(ida_, LS_klu_, A_sp_) != IDA_SUCCESS
        || IDASetJacFn(ida_, &IdaSolver::jacobian_callback) != IDA_SUCCESS) {
        cleanup_klu();
        return false;
    }

    const double ms =
        std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - t_start)
            .count();
    std::printf("[IdaSolver] KLU 稀疏直解：n=%d nnz=%lld 颜色数=%d "
                "(模式探测+着色 %.0f ms)\n",
                n, static_cast<long long>(nnz), ncolors, ms);
    return true;
}

int IdaSolver::jacobian_callback(real_t t, real_t cj, N_Vector y, N_Vector yp,
                                 N_Vector rr, SUNMatrix Jac, void* user_data,
                                 N_Vector /*tmp1*/, N_Vector /*tmp2*/,
                                 N_Vector /*tmp3*/) {
    auto* self = static_cast<IdaSolver*>(user_data);
    const int n = self->problem_.dae_size();
    const auto* y_ptr = N_VGetArrayPointer(y);
    const auto* yp_ptr = N_VGetArrayPointer(yp);
    const auto* rr_ptr = N_VGetArrayPointer(rr);
    auto* data = SUNSparseMatrix_Data(Jac);
    auto* colptr = SUNSparseMatrix_IndexPointers(Jac);
    auto* rowidx = SUNSparseMatrix_IndexValues(Jac);
    if (y_ptr == nullptr || yp_ptr == nullptr || rr_ptr == nullptr
        || data == nullptr || colptr == nullptr || rowidx == nullptr) {
        return -1; // 不可恢复错误
    }

    // IDA 在每次装配前对 DIRECT 求解器执行 SUNMatZero，而
    // SUNMatZero_Sparse 会把 CSC 的 colptr/rowidx 一并清零（见
    // sundials_sparse.c SUNMatZero_Sparse），因此每次回调都必须
    // 重建稀疏模式（SUNDIALS KLU 样例同做法）
    std::copy(self->sp_colptr_.begin(), self->sp_colptr_.end(), colptr);
    std::copy(self->sp_rowidx_.begin(), self->sp_rowidx_.end(), rowidx);

    const Eigen::Map<const Vector> y0(y_ptr, n);
    const Eigen::Map<const Vector> yp0(yp_ptr, n);
    const Eigen::Map<const Vector> f0(rr_ptr, n); // = F(t, y, yp)

    Vector& y_t = self->jac_y_buf_;
    Vector& yp_t = self->jac_yp_buf_;
    Vector& f1 = self->jac_f1_buf_;
    Vector& sigmas = self->jac_sigma_buf_;
    y_t = y0;
    yp_t = yp0;

    // 同色列一次差分：扰动 y_j += σ_j、yp_j += cj·σ_j，
    // J_ij = ΔF_i / σ_j（σ 与模式探测时同式，随当前状态缩放）
    for (const auto& group : self->sp_groups_) {
        for (int j : group) {
            const real_t scale =
                std::max(std::abs(y0(j)), self->atol_eigen_(j));
            const real_t sigma =
                1e-4 * (scale > 0.0 ? scale : real_t(1.0));
            sigmas(j) = sigma;
            y_t(j) += sigma;
            yp_t(j) += cj * sigma;
        }
        self->problem_.dae_residual(t, y_t, yp_t, f1);
        if (!f1.allFinite()) {
            return 1; // 可恢复：IDA 减步长重试
        }
        for (int j : group) {
            const real_t inv = 1.0 / sigmas(j);
            const auto& rows = self->sp_rows_of_col_[j];
            const auto& slots = self->sp_slots_[j];
            for (std::size_t k = 0; k < rows.size(); ++k) {
                data[slots[k]] = (f1(rows[k]) - f0(rows[k])) * inv;
            }
        }
        for (int j : group) {
            y_t(j) = y0(j);
            yp_t(j) = yp0(j);
        }
    }

    // 差分商异常（爆炸/非有限值）时返回可恢复失败，交由 IDA 减步长重试
    for (sunindextype k = 0; k < self->sp_colptr_[n]; ++k) {
        const real_t v = std::abs(data[k]);
        if (!std::isfinite(v) || v > 1e30) {
            return 1;
        }
    }
    return 0;
}

#endif // FLUXCELL_WITH_SUNDIALS_KLU

IdaSolver::~IdaSolver() {
#if defined(FLUXCELL_WITH_SUNDIALS_KLU)
    cleanup_klu();
#endif
    if (ida_ != nullptr) {
        IDAFree(&ida_);
    }
    if (y_ != nullptr) {
        N_VDestroy(y_);
    }
    if (yp_ != nullptr) {
        N_VDestroy(yp_);
    }
    if (id_ != nullptr) {
        N_VDestroy(id_);
    }
    if (atol_ != nullptr) {
        N_VDestroy(atol_);
    }
    if (ctx_ != nullptr) {
        SUNContext_Free(&ctx_);
    }
}

int IdaSolver::residual_callback(real_t t, N_Vector y, N_Vector yp,
                                 N_Vector res, void* user_data) {
    auto* self = static_cast<IdaSolver*>(user_data);
    const auto* y_ptr = N_VGetArrayPointer(y);
    const auto* yp_ptr = N_VGetArrayPointer(yp);
    auto* res_ptr = N_VGetArrayPointer(res);
    if (y_ptr == nullptr || yp_ptr == nullptr || res_ptr == nullptr) {
        return -1; // 不可恢复错误
    }

    const int n = self->problem_.dae_size();
    // 借用 N_Vector 内存构造 Eigen 视图（避免每步拷贝）
    Eigen::Map<const Vector> y_map(y_ptr, n);
    Eigen::Map<const Vector> yp_map(yp_ptr, n);
    Eigen::Map<Vector> res_map(res_ptr, n);

    self->problem_.dae_residual(t, y_map, yp_map, res_map);
    return res_map.allFinite() ? 0 : 1; // 非有限值 = 可恢复错误（IDA 减步长重试）
}

void IdaSolver::advance_to(real_t t_out) {
    if (IDASetStopTime(ida_, t_out) != IDA_SUCCESS) {
        throw std::runtime_error("IdaSolver: IDASetStopTime 失败");
    }

    real_t t_reached = 0.0;
    const int flag = IDASolve(ida_, t_out, &t_reached, y_, yp_, IDA_NORMAL);
    // IDA_TSTOP_RETURN：到达 IDASetStopTime 设定的输出点（与 IDA_SUCCESS 同为成功）
    if (flag != IDA_SUCCESS && flag != IDA_TSTOP_RETURN) {
        const char* reason = flag == IDA_TOO_MUCH_WORK       ? "输出区间内步数超限"
                             : flag == IDA_CONV_FAIL         ? "Newton 收敛失败"
                             : flag == IDA_LSETUP_FAIL       ? "线性求解器装配失败"
                             : flag == IDA_LSOLVE_FAIL       ? "线性求解器求解失败"
                             : flag == IDA_CONSTR_FAIL       ? "约束校验失败"
                             : flag == IDA_ILL_INPUT         ? "非法输入"
                                                              : "未知错误";
        // 诊断：失败前最后成功点的状态摘要（定位 corrector 失败的状态区）
        problem_.dae_dump_state(t_reached,
                                Eigen::Map<const Vector>(
                                    N_VGetArrayPointer(y_),
                                    problem_.dae_size()));
        throw std::runtime_error("IdaSolver: 推进到 t=" + std::to_string(t_out)
                                 + " 失败：" + reason + " (IDA flag "
                                 + std::to_string(flag) + ")");
    }

    sunrealtype h_last = 0.0;
    IDAGetLastStep(ida_, &h_last);
    last_step_ = h_last;
    ++total_steps_;

    // 状态写回模型（导出量同步重算）
    problem_.dae_load_state(t_reached, y_eigen_);
}

void IdaSolver::reinit(real_t t, real_t dt_hint) {
    // 用当前模型状态重新打包（导数 = f(y)，保证 yp 与新残差一致）
    problem_.dae_pack_state(y_eigen_, yp_eigen_);

    if (IDAReInit(ida_, t, y_, yp_) != IDA_SUCCESS) {
        throw std::runtime_error("IdaSolver: IDAReInit 失败");
    }

    // 段电流切换后代数变量（i_surf/δ）不满足新约束，需一致化；
    // 纯 ODE 问题（SPM/SPMe）初值天然一致，无需 CalcIC
    if (has_algebraic_ && dt_hint > 0.0) {
        // 大电流跳变（如 2C 充电）下 Newton 初值离解较远，CalcIC 可能不收敛。
        // 逐级缩小 tout1 重试：更小的时间尺度 → 更保守的初始步长 →
        // 微分变量近似冻结，代数方程（BV + 守恒）线性化良好，易于收敛。
        // 每次重试前重新打包/ReInit，丢弃上次失败 CalcIC 对初值的污染
        int flag = IDA_ILL_INPUT;
        real_t hint = dt_hint;
        for (int attempt = 0; attempt < 4; ++attempt) {
            problem_.dae_pack_state(y_eigen_, yp_eigen_);
            if (IDAReInit(ida_, t, y_, yp_) != IDA_SUCCESS) {
                throw std::runtime_error("IdaSolver: IDAReInit 失败");
            }
            flag = IDACalcIC(ida_, IDA_YA_YDP_INIT, t + hint);
            if (flag == IDA_SUCCESS) {
                break;
            }
            hint *= 0.01;
        }
        if (flag != IDA_SUCCESS) {
            throw std::runtime_error(
                "IdaSolver: IDACalcIC 初始一致化失败 (IDA flag "
                + std::to_string(flag) + ")");
        }
        // CalcIC 校正后的状态写回模型
        problem_.dae_load_state(t, y_eigen_);
    }
}

real_t IdaSolver::time() const {
    sunrealtype t = 0.0;
    IDAGetCurrentTime(ida_, &t);
    return t;
}

real_t IdaSolver::last_step() const {
    return last_step_;
}

long IdaSolver::total_steps() const {
    return total_steps_;
}

} // namespace fluxcell::solver

#endif // FLUXCELL_WITH_SUNDIALS
