/// @file ida_solver.hpp
/// @brief SUNDIALS IDA 自适应 DAE 求解器包装（仅 FLUXCELL_WITH_SUNDIALS 编译）
///
/// 封装 IDA 内存管理、N_Vector/Eigen 交换与输出点推进：
/// - advance_to(t_out)：内部自适应步长，由 IDASetStopTime 保证精确停在输出点
/// - reinit(t, dt_hint)：段电流切换等残差突变后重置内部历史，
///   含代数变量时用 IDACalcIC 重新一致化（IDARerinit 后历史失效）
#pragma once

#include <fluxcell/core/solver/dae_problem.hpp>
#include <fluxcell/core/types.hpp>

#if defined(FLUXCELL_WITH_SUNDIALS)

#include <sundials/sundials_context.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>

#include <vector>

namespace fluxcell::solver {

/// IDA 自适应求解器（RAII 包装；同一实例仅允许单线程使用）
class IdaSolver {
public:
    /// @param problem  DAE 问题（生命周期须长于求解器）
    /// @param t0       初始时间 [s]
    /// @param rel_tol  相对误差容限
    IdaSolver(DaeProblem& problem, real_t t0, real_t rel_tol);
    ~IdaSolver();

    IdaSolver(const IdaSolver&) = delete;
    IdaSolver& operator=(const IdaSolver&) = delete;

    /// 自适应推进到 t_out（内部步长自动缩放，输出点精确停住）
    /// @throws std::runtime_error 求解失败（收敛失败等）
    void advance_to(real_t t_out);

    /// 重置内部历史（保留状态）。含代数变量时执行 IDACalcIC 一致化
    /// @param t       重置时刻（通常为当前模型时间）
    /// @param dt_hint 首步时间尺度提示（一般取输出步长）
    void reinit(real_t t, real_t dt_hint);

    /// 当前求解时间 [s]
    real_t time() const;

    /// 最近一个内部步长 [s]（诊断：观察自适应行为）
    real_t last_step() const;

    /// 累计内部步数（诊断：自适应效率）
    long total_steps() const;

private:
    /// IDA 残差回调（转发到 DaeProblem::dae_residual）
    static int residual_callback(real_t t, N_Vector y, N_Vector yp,
                                 N_Vector res, void* user_data);

#if defined(FLUXCELL_WITH_SUNDIALS_KLU)
    /// IDA Jacobian 回调（KLU 路径：按颜色组同时扰动多列差分填充 CSC）
    static int jacobian_callback(real_t t, real_t cj, N_Vector y, N_Vector yp,
                                 N_Vector rr, SUNMatrix Jac, void* user_data,
                                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

    /// 构建稀疏模式（逐列扰动探测）+ 距离-1 着色 + 挂接 KLU 线性求解器；
    /// 任一环节失败返回 false（调用方回退 GMRES）
    bool setup_klu(real_t t0);

    /// 释放 KLU 稀疏路径对象（析构与回退共用）
    void cleanup_klu();
#endif

    DaeProblem& problem_;

    SUNContext ctx_ = nullptr;   ///< SUNDIALS 上下文（每实例独立）
    void* ida_ = nullptr;        ///< IDA 内存块
    N_Vector y_ = nullptr;       ///< 状态向量（包裹 y_eigen_ 内存）
    N_Vector yp_ = nullptr;      ///< 导数向量（包裹 yp_eigen_ 内存）
    N_Vector id_ = nullptr;      ///< 微分/代数变量掩码（无代数变量时为空）
    N_Vector atol_ = nullptr;    ///< 绝对容限向量

    // N_VMake 包裹的交换缓冲（N_Vector 直接读写其内存，零拷贝；
    // 生命周期与求解器一致，指针必须保持稳定）
    Vector y_eigen_;
    Vector yp_eigen_;
    Vector res_eigen_;
    Vector mask_eigen_; ///< 变量类型掩码缓冲（id_ 引用其内存）
    Vector atol_eigen_; ///< 绝对容限缓冲（atol_ 引用其内存）

    bool has_algebraic_ = false; ///< 是否含代数变量（决定 CalcIC 策略）
    real_t last_step_ = 0.0;
    long total_steps_ = 0;

#if defined(FLUXCELL_WITH_SUNDIALS_KLU)
    // ── KLU 稀疏直解路径（PyBaMM IDAKLUSolver 同构）────────────
    //  Jacobian 由差分填充：逐列扰动探测出模式后，同色列（不共享
    //  残差行）可一次评估同时填充，成本 ~ 颜色数次残差评估
    bool sparse_active_ = false;         ///< 是否已启用 KLU 稀疏直解
    SUNMatrix A_sp_ = nullptr;           ///< CSC 稀疏 Jacobian 矩阵
    SUNLinearSolver LS_klu_ = nullptr;   ///< KLU 线性求解器
    std::vector<sunindextype> sp_colptr_;           ///< CSC 列指针（n+1）
    std::vector<sunindextype> sp_rowidx_;           ///< CSC 行索引（列内升序）
    std::vector<std::vector<int>> sp_rows_of_col_;  ///< 每列依赖的残差行
    std::vector<std::vector<int>> sp_slots_;        ///< 每列各行在 data 中的槽位
    std::vector<int> sp_color_;                     ///< 每列颜色（距离-1 着色）
    std::vector<std::vector<int>> sp_groups_;       ///< 每颜色的列集合
    // 差分回调的工作缓冲（避免每次回调堆分配）
    Vector jac_y_buf_, jac_yp_buf_, jac_f1_buf_, jac_sigma_buf_;
#endif
};

} // namespace fluxcell::solver

#endif // FLUXCELL_WITH_SUNDIALS
