/// @file dfn_model.cpp
/// @brief 全阶 DFN (P2D) 模型实现：2D 固相颗粒 + 电解液 1D 传输 + 逐格 BV 反应分布
#include <fluxcell/core/model/dfn_model.hpp>

#include <fluxcell/core/discretization/discretization.hpp>
#include <fluxcell/core/physics/physics.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fluxcell::model {

namespace {
using fluxcell::physics::exchange_current_density_chen2020_neg;
using fluxcell::physics::exchange_current_density_chen2020_pos;
using fluxcell::physics::graphite_chen2020_ocp;
using fluxcell::physics::nmc_chen2020_ocp;

constexpr real_t kFaraday = 96485.0;  // C/mol
constexpr real_t kPi = 3.14159265358979323846;
constexpr real_t kGasConst = 8.314;   // J/(mol·K)
constexpr real_t kCeFloor = 1.0;      // 电解液浓度数值下限 [mol/m³]
constexpr real_t kKappaFloor = 1e-6;  // κ 数值下限 [S/m]（盐耗尽 c_e→0 防除零）

/// BV 的 sinh 有限防护：bv_coef·η 超过 ~700 时 sinh→inf 污染残差/Gummel；
/// 钳位值下电流已远超任何工作电流，只影响远离解的试探步
inline real_t safe_sinh(real_t x) {
    return std::sinh(std::clamp(x, -700.0, 700.0));
}

/// 产热计算的场输入（CN 路径用成员场、DAE 残差用局部状态向量共用）
struct HeatSourceArgs {
    int m = 0;             ///< 厚度方向格数
    int m_neg = 0;         ///< 负极格数
    int pos_start = 0;     ///< 正极区起始格索引
    real_t dx = 0.0;       ///< 厚度方向步长 [m]
    real_t area = 0.0;     ///< 电极面积 A [m²]
    real_t a_n = 0.0;      ///< 负极比表面积 [m²/m³]
    real_t a_p = 0.0;      ///< 正极比表面积 [m²/m³]
    real_t sigma_n_eff = 0.0;
    real_t sigma_p_eff = 0.0;
    real_t i_app = 0.0;    ///< 端口面电流密度 [A/m²]
    real_t T = 0.0;        ///< 电芯温度 [K]
    real_t entropic_n = 0.0; ///< dU_n/dT 常数 [V/K]（LUT 缺失时 fallback）
    real_t entropic_p = 0.0; ///< dU_p/dT 常数 [V/K]
    const EntropyLut* lut_n = nullptr; ///< dU/dT(sto) 表（非空时优先于常数）
    const EntropyLut* lut_p = nullptr;
    const Vector* c_surf_n = nullptr; ///< 负极各格表面浓度 [mol/m³]（LUT 求值）
    const Vector* c_surf_p = nullptr; ///< 正极各格表面浓度 [mol/m³]
    real_t cs_max_n = 0.0; ///< 负极最大嵌入浓度 [mol/m³]（sto 归一用）
    real_t cs_max_p = 0.0; ///< 正极最大嵌入浓度 [mol/m³]
    const Vector* i_e_face = nullptr; ///< 液相面电流密度（m+1）[A/m²]
    const Vector* phi_e = nullptr;    ///< 液相电位分布 [V]
    const Vector* i_surf_n = nullptr; ///< 负极孔壁电流密度 [A/m²]
    const Vector* i_surf_p = nullptr;
    const Vector* eta_n = nullptr;    ///< 负极过电位分布 [V]
    const Vector* eta_p = nullptr;
    const Vector* j_sei_n = nullptr;  ///< 负极 SEI 副反应电流 [A/m²]（可空）
    const Vector* eta_sei_n = nullptr; ///< 负极 SEI 反应过电位 [V]（可空）
};

/// 瞬时产热功率（PyBaMM base_thermal 三项口径，逐格累乘 dx·A）[W]：
/// - Q_ohm：固相 i_s²/σ^eff + 液相 −i_e·∂φ_e/∂x（后者已含浓差扩散贡献）
/// - Q_rxn：不可逆界面热 a·i_surf·η
/// - Q_rev：可逆熵热 a·i_surf·T·dU/dT
/// 集流体欧姆热忽略（Chen2020 集流体面电导 ~700 S，贡献 < 1e-9 W）
void compute_heat_source(const HeatSourceArgs& a,
                         real_t& q_ohm, real_t& q_rxn, real_t& q_rev) {
    q_ohm = 0.0;
    q_rxn = 0.0;
    q_rev = 0.0;
    for (int i = 0; i < a.m; ++i) {
        real_t q = 0.0;
        if (i < a.m_neg) {
            const real_t i_s = a.i_app - (*a.i_e_face)(i);
            q += i_s * i_s / a.sigma_n_eff;
        } else if (i >= a.pos_start) {
            const real_t i_s = a.i_app - (*a.i_e_face)(i);
            q += i_s * i_s / a.sigma_p_eff;
        }
        if (i > 0) {
            q += -(*a.i_e_face)(i)
                 * ((*a.phi_e)(i) - (*a.phi_e)(i - 1)) / a.dx;
        }
        if (i < a.m_neg) {
            const real_t a_j = a.a_n * (*a.i_surf_n)(i); // 体积电流密度 [A/m³]
            q_rxn += a_j * (*a.eta_n)(i) * a.dx;
            // 可逆熵热：dU/dT 按各格表面嵌入率 sto 插值（LUT 空时用常数）
            const real_t dudt = (a.lut_n != nullptr && a.cs_max_n > 0.0)
                                    ? a.lut_n->eval((*a.c_surf_n)(i) / a.cs_max_n)
                                    : a.entropic_n;
            q_rev += a_j * a.T * dudt * a.dx;
            // SEI 副反应不可逆热（小量，等温对标不影响；lumped 时并入 q_rxn）
            if (a.j_sei_n != nullptr) {
                q_rxn += a.a_n * (*a.j_sei_n)(i) * (*a.eta_sei_n)(i) * a.dx;
            }
        } else if (i >= a.pos_start) {
            const int k = i - a.pos_start;
            const real_t a_j = a.a_p * (*a.i_surf_p)(k);
            q_rxn += a_j * (*a.eta_p)(k) * a.dx;
            const real_t dudt = (a.lut_p != nullptr && a.cs_max_p > 0.0)
                                    ? a.lut_p->eval((*a.c_surf_p)(k) / a.cs_max_p)
                                    : a.entropic_p;
            q_rev += a_j * a.T * dudt * a.dx;
        }
        q_ohm += q * a.dx;
    }
    q_ohm *= a.area;
    q_rxn *= a.area;
    q_rev *= a.area;
}
} // namespace

DFNModel::DFNModel(const Parameters& params, real_t initial_current)
    : params_(params) {
    const auto& solid = params.solid;

    n_ = solid.particle_cells;
    if (n_ < 3) {
        throw std::invalid_argument("DFNModel: particle_cells must be >= 3");
    }
    m_ = params.electrolyte_cells;
    if (m_ < 10) {
        throw std::invalid_argument("DFNModel: electrolyte_cells must be >= 10");
    }
    if (solid.temperature <= 0.0) {
        throw std::invalid_argument("DFNModel: temperature must be positive");
    }
    if (solid.initial_soc < 0.0 || solid.initial_soc > 1.0) {
        throw std::invalid_argument("DFNModel: initial_soc must be in [0, 1]");
    }
    if (params.porosity_neg <= 0.0 || params.porosity_sep <= 0.0
        || params.porosity_pos <= 0.0) {
        throw std::invalid_argument("DFNModel: porosity must be positive");
    }
    if (params.electrolyte_conductivity <= 0.0
        || params.electrolyte_diffusivity <= 0.0) {
        throw std::invalid_argument("DFNModel: electrolyte transport must be positive");
    }
    params.validate_transport_table();
    if (params.neg_conductivity <= 0.0 || params.pos_conductivity <= 0.0) {
        throw std::invalid_argument("DFNModel: solid conductivity must be positive");
    }
    if (params.max_iterations < 1 || params.tolerance <= 0.0) {
        throw std::invalid_argument("DFNModel: invalid iteration settings");
    }
    if (params.thermal_type == ThermalType::Lumped) {
        const auto& th = params.thermal_params;
        if (!(th.t_ambient > 0.0) || !(th.t_initial > 0.0)
            || th.h_total < 0.0 || th.cooling_area < 0.0
            || th.cell_volume <= 0.0
            || th.rho_c_p_neg <= 0.0 || th.rho_c_p_sep <= 0.0
            || th.rho_c_p_pos <= 0.0 || th.rho_c_p_cc_neg <= 0.0
            || th.rho_c_p_cc_pos <= 0.0
            || th.cc_thickness_neg < 0.0 || th.cc_thickness_pos < 0.0) {
            throw std::invalid_argument("DFNModel: invalid thermal parameters");
        }
    }
    if (params.degradation.enabled) {
        const auto& deg = params.degradation;
        if (!(deg.l_sei_0 >= 0.0) || deg.v_bar_sei <= 0.0 || deg.z_sei <= 0.0
            || deg.r_sei < 0.0 || deg.k_sei <= 0.0 || deg.c_ec_0 <= 0.0
            || deg.d_ec <= 0.0 || deg.t_ref <= 0.0
            || !(deg.alpha_sei > 0.0 && deg.alpha_sei < 1.0)) {
            throw std::invalid_argument("DFNModel: invalid degradation parameters");
        }
    }

    area_ = solid.electrode_height * solid.electrode_width;
    a_n_ = 3.0 * solid.neg_active_fraction / solid.neg_particle_radius;
    a_p_ = 3.0 * solid.pos_active_fraction / solid.pos_particle_radius;

    // ── x 方向网格：均匀 dx，区域边界对齐最近格边界（同 SPMe）──
    const real_t L = solid.neg_thickness + solid.sep_thickness + solid.pos_thickness;
    dx_ = L / m_;

    const auto cells_for = [&](real_t thickness) {
        return std::clamp(static_cast<int>(std::lround(thickness / dx_)), 1, m_ - 2);
    };
    sep_start_ = cells_for(solid.neg_thickness);
    pos_start_ = std::clamp(sep_start_ + cells_for(solid.sep_thickness),
                            sep_start_ + 1, m_ - 1);
    m_neg_ = sep_start_;
    m_pos_ = m_ - pos_start_;

    eps_e_ = Vector::Constant(m_, params.porosity_sep);
    eps_e_.segment(0, sep_start_).setConstant(params.porosity_neg);
    eps_e_.segment(pos_start_, m_ - pos_start_).setConstant(params.porosity_pos);

    // 有效输运参数（固相电导缓存：孔隙率与温度不变；液相 κ/D 浓度相关，
    // 在电位积分/DAE 残差内逐格求值）
    const real_t brug = params.bruggeman_exponent;
    sigma_n_eff_ = params.neg_conductivity
                 * std::pow(solid.neg_active_fraction, brug);
    sigma_p_eff_ = params.pos_conductivity
                 * std::pow(solid.pos_active_fraction, brug);

    // ── 热模型初始化（须先于初始反应分布求解：BV/i0 依赖 t_cell_）──
    // 等温：温度恒为参数温度；集总：状态变量，初值 = 初始温度
    rho_c_p_eff_ = rho_c_p_eff();
    t_cell_ = (params.thermal_type == ThermalType::Lumped)
            ? params.thermal_params.t_initial : solid.temperature;

    // ── 颗粒径向网格 + CN 算子素材（同 SPM）──
    auto mesh_neg = discretization::SphericalMesh::create(solid.neg_particle_radius, n_);
    auto mesh_pos = discretization::SphericalMesh::create(solid.pos_particle_radius, n_);

    D_neg_ = discretization::assemble_spherical_diffusion(mesh_neg, solid.neg_diffusion);
    D_pos_ = discretization::assemble_spherical_diffusion(mesh_pos, solid.pos_diffusion);

    surf_area_n_ = mesh_neg.face_areas()[n_];
    surf_vol_n_ = mesh_neg.volumes()[n_ - 1];
    surf_area_p_ = mesh_pos.face_areas()[n_];
    surf_vol_p_ = mesh_pos.volumes()[n_ - 1];
    vol_neg_.assign(mesh_neg.volumes().begin(), mesh_neg.volumes().end());
    vol_pos_.assign(mesh_pos.volumes().begin(), mesh_pos.volumes().end());

    // ── 2D 初始固相浓度：每列颗粒相同，由初始 SOC 线性映射端点（同 SPM）──
    const real_t x0 = solid.neg_stoich_at_soc0
                    + solid.initial_soc
                          * (solid.neg_stoich_at_soc1 - solid.neg_stoich_at_soc0);
    const real_t y0 = solid.pos_stoich_at_soc0
                    + solid.initial_soc
                          * (solid.pos_stoich_at_soc1 - solid.pos_stoich_at_soc0);
    const Vector init_n = Vector::Constant(n_, x0 * solid.neg_max_conc);
    const Vector init_p = Vector::Constant(n_, y0 * solid.pos_max_conc);
    c_s_neg_.resize(n_, m_neg_);
    c_s_pos_.resize(n_, m_pos_);
    for (int k = 0; k < m_neg_; ++k) {
        c_s_neg_.col(k) = init_n;
    }
    for (int k = 0; k < m_pos_; ++k) {
        c_s_pos_.col(k) = init_p;
    }

    // ── 液相初始状态与场变量 ──
    c_e_.setConstant(m_, params.c_e_init);
    phi_e_.setZero(m_);
    phi_s_n_.setZero(m_neg_);
    phi_s_p_.setZero(m_pos_);
    phi_s_full_.setZero(m_);
    i_e_face_.setZero(m_ + 1);
    i_surf_n_.setZero(m_neg_);
    i_surf_p_.setZero(m_pos_);
    eta_n_.setZero(m_neg_);
    eta_p_.setZero(m_pos_);
    c_surf_n_.setZero(m_neg_);
    c_surf_p_.setZero(m_pos_);
    u_n_vec_.setZero(m_neg_);
    u_p_vec_.setZero(m_pos_);

    dae_current_ = initial_current;
    i_port_ = initial_current; // 增广端口初值（pack 的 I_cell）

    cn_dt_ = -1.0;

    // SEI 降解状态：初始厚度逐格 = l_sei_0，副反应电流场在反应分布
    // 求解的迭代内建立；此处先零初始化（j_sei 进入 Gummel 初值）
    if (params_.degradation.enabled) {
        l_sei_.setConstant(m_neg_, params_.degradation.l_sei_0);
        j_sei_n_.setZero(m_neg_);
        eta_sei_n_.setZero(m_neg_);
    }

    refresh_surface_quantities();
    solve_reaction_distribution(initial_current / area_);

    // 初始产热缓存（t = 0 快照诊断用；等温模式同样计算，仅不更新温度）
    refresh_heat(initial_current / area_);
}

DFNModel::~DFNModel() = default;

Vector DFNModel::electrolyte_div_flux(const Vector& c) const {
    // 守恒通量形式三对角散度（同 SPMe）：D_e^eff = D(c_e)·ε^bruggeman（逐格），
    // 界面算术平均，两端零通量（离散与常量算子版本一致）
    Vector deff(m_);
    for (int i = 0; i < m_; ++i) {
        deff(i) = params_.diffusivity_at(c(i))
                * std::pow(eps_e_(i), params_.bruggeman_exponent);
    }
    const real_t inv_dx2 = 1.0 / (dx_ * dx_);
    Vector r = Vector::Zero(m_);
    for (int i = 0; i < m_; ++i) {
        const real_t d_left = (i > 0) ? 0.5 * (deff(i - 1) + deff(i)) : 0.0;
        const real_t d_right = (i < m_ - 1) ? 0.5 * (deff(i) + deff(i + 1)) : 0.0;
        const real_t flux_left = (i > 0) ? d_left * (c(i - 1) - c(i)) : 0.0;
        const real_t flux_right = (i < m_ - 1) ? d_right * (c(i + 1) - c(i)) : 0.0;
        r(i) = (flux_left + flux_right) * inv_dx2;
    }
    return r;
}

real_t DFNModel::kappa_eff_at(const Eigen::Ref<const Vector>& c_e, int i) const {
    // κ^eff = κ(c_e)·ε^bruggeman；κ 下限防盐耗尽（c_e → 0）欧姆项除零
    return std::max(params_.conductivity_at(c_e(i)), kKappaFloor)
         * std::pow(eps_e_(i), params_.bruggeman_exponent);
}

void DFNModel::ensure_cn(real_t dt) {
    const real_t half_dt = 0.5 * dt;

    if (cn_dt_ != dt) {
        // ── 颗粒算子（每电极一份，各列共用，同 SPM；dt 不变时复用）──
        auto build_particle = [this, half_dt](const SparseMatrix& D, SparseMatrix& plus,
                                              Eigen::SparseLU<SparseMatrix>& solve) {
            SparseMatrix lhs = -half_dt * D;
            for (int i = 0; i < n_; ++i) {
                lhs.coeffRef(i, i) += 1.0;
            }
            plus = half_dt * D;
            for (int i = 0; i < n_; ++i) {
                plus.coeffRef(i, i) += 1.0;
            }
            solve.compute(lhs);
            if (solve.info() != Eigen::Success) {
                throw std::runtime_error("DFNModel: particle CN factorization failed");
            }
        };
        build_particle(D_neg_, cn_plus_neg_, cn_solve_neg_);
        build_particle(D_pos_, cn_plus_pos_, cn_solve_pos_);

        cn_dt_ = dt;
    }

    // ── 液相算子：ε dc/dt = L(c_e)·c + s，L 为守恒通量形式三对角（同 SPMe）──
    // D_e 依赖 c_e（浓度相关表），每步以当前浓度重建（滞后显式；
    // IDA 路径在残差内全隐式求值），界面取算术平均，两端零通量
    Vector deff(m_);
    for (int i = 0; i < m_; ++i) {
        deff(i) = params_.diffusivity_at(c_e_(i))
                * std::pow(eps_e_(i), params_.bruggeman_exponent);
    }
    const real_t inv_dx2 = 1.0 / (dx_ * dx_);

    TripletList lhs_triplets;
    TripletList plus_triplets;
    lhs_triplets.reserve(3 * static_cast<std::size_t>(m_));
    plus_triplets.reserve(3 * static_cast<std::size_t>(m_));

    for (int i = 0; i < m_; ++i) {
        const real_t d_left = (i > 0) ? 0.5 * (deff(i - 1) + deff(i)) : 0.0;
        const real_t d_right = (i < m_ - 1) ? 0.5 * (deff(i) + deff(i + 1)) : 0.0;

        lhs_triplets.emplace_back(i, i, eps_e_(i)
                                            + half_dt * (d_left + d_right) * inv_dx2);
        plus_triplets.emplace_back(i, i, eps_e_(i)
                                             - half_dt * (d_left + d_right) * inv_dx2);
        if (i > 0) {
            lhs_triplets.emplace_back(i, i - 1, -half_dt * d_left * inv_dx2);
            plus_triplets.emplace_back(i, i - 1, half_dt * d_left * inv_dx2);
        }
        if (i < m_ - 1) {
            lhs_triplets.emplace_back(i, i + 1, -half_dt * d_right * inv_dx2);
            plus_triplets.emplace_back(i, i + 1, half_dt * d_right * inv_dx2);
        }
    }

    le_lhs_.resize(m_, m_);
    le_plus_.resize(m_, m_);
    le_lhs_.setFromTriplets(lhs_triplets.begin(), lhs_triplets.end());
    le_plus_.setFromTriplets(plus_triplets.begin(), plus_triplets.end());
    le_solve_.compute(le_lhs_);
    if (le_solve_.info() != Eigen::Success) {
        throw std::runtime_error("DFNModel: electrolyte CN factorization failed");
    }
}

void DFNModel::dae_set_current(real_t current) {
    if (port_control_ == solver::PortControl::Current
        && current == dae_current_) {
        return; // 恒流段内重复调用：构造/上次切换已一致化，短路
    }
    port_control_ = solver::PortControl::Current;
    dae_current_ = current;
    // 电流切换时解析重估代数变量（i_surf 分布与 δ）：大电流跳变
    //（如 2C 充电）时 IDACalcIC 的 Newton 从旧状态出发难以收敛
    //（初值距解太远，线搜索陷入非有限残差区），而半解析 Gummel
    // 迭代对跳变稳健（与固定步长路径共用同一逻辑）。
    // 重估后 dae_pack_state 打包的 y 即满足代数约束，CalcIC 零修正收敛
    refresh_surface_quantities();
    solve_reaction_distribution(current / area_);
    refresh_heat(current / area_); // 产热缓存与重估状态同步
}

void DFNModel::dae_set_voltage(real_t voltage) {
    if (port_control_ == solver::PortControl::Voltage
        && voltage == v_target_) {
        return; // 恒压段内重复调用：已一致化，短路
    }
    // 纯切换不重估：当前状态是 IDA 刚推进完的自洽解（i_surf 满足
    // BV/守恒行 @ I_cell），比冻结浓度的 Gummel 半解析解更自洽；
    // CC→CV 切换发生在电压截止（V 已 ≈ 目标），外电路行初残差
    // 仅 mV 级，CalcIC 小修正即可收敛。探测场景（IdaSolver 构造）
    // 同样只需行形态切换
    port_control_ = solver::PortControl::Voltage;
    v_target_ = voltage;
}

void DFNModel::dae_dump_state(real_t t,
                              const Eigen::Ref<const Vector>& y) const {
    const int s_n = n_ * m_neg_;
    const int s_p = n_ * m_pos_;
    const int off_e = s_n + s_p;
    const int off_isn = off_e + m_;
    const int off_isp = off_isn + m_neg_;
    const int off_d = off_isp + m_pos_;
    const bool lumped = params_.thermal_type == ThermalType::Lumped;
    const int off_sei = off_d + 2 + (lumped ? 1 : 0);

    std::fprintf(stderr,
                 "[dfn-dump] t=%.3f  c_n=[%.6g, %.6g]  c_p=[%.6g, %.6g]  "
                 "c_e=[%.6g, %.6g]\n",
                 t, y.segment(0, s_n).minCoeff(), y.segment(0, s_n).maxCoeff(),
                 y.segment(s_n, s_p).minCoeff(),
                 y.segment(s_n, s_p).maxCoeff(),
                 y.segment(off_e, m_).minCoeff(),
                 y.segment(off_e, m_).maxCoeff());
    std::fprintf(stderr,
                 "[dfn-dump]   is_n=[%.6g, %.6g]  is_p=[%.6g, %.6g]  "
                 "d_n=%.6f d_p=%.6f\n",
                 y.segment(off_isn, m_neg_).minCoeff(),
                 y.segment(off_isn, m_neg_).maxCoeff(),
                 y.segment(off_isp, m_pos_).minCoeff(),
                 y.segment(off_isp, m_pos_).maxCoeff(),
                 y(off_d), y(off_d + 1));
    std::fprintf(stderr, "[dfn-dump]   I_cell=%.6g A  control=%s\n",
                 y(y.size() - 1),
                 port_control_ == solver::PortControl::Voltage ? "CV"
                                                               : "CC");
    if (params_.degradation.enabled) {
        std::fprintf(stderr, "[dfn-dump]   L_sei=[%.6e, %.6e]\n",
                     y.segment(off_sei, m_neg_).minCoeff(),
                     y.segment(off_sei, m_neg_).maxCoeff());
    }
    if (lumped) {
        std::fprintf(stderr, "[dfn-dump]   T=%.4f K\n", y(off_d + 2));
    }
}

void DFNModel::dae_jac_col_prune(int j, std::vector<int>& rows) const {
    // SEI 厚度列：仅保留 SEI 行对角。K 饱和区 j_sei ∝ 1/L 的强依赖经
    // 守恒行（a·dx·∂j_sei/∂L ~ 3e4/列）、BV 行（φ 链 ~ 3e5/列）与热行
    // （反应热 ~ 2e4/列）差分出巨元，而热行/守恒行的对角仅 ~1e-3/
    // 结构性 0 → 行对角劣势 1e-6~1e-9 → KLU pivot 失效（corrector
    // 全灭，lumped+SEI 200 圈第 2 圈 CV 段 flag -4 实证）。被删反馈
    // 对 Newton 的实际影响 = 巨元×dz_L（~1e-11）≈ 1e-6，高阶小量；
    // j_sei↔φ 回路由残差内不动点精确处理，裁剪后 quasi-Newton 收敛
    if (!params_.degradation.enabled) {
        return;
    }
    const bool lumped = params_.thermal_type == ThermalType::Lumped;
    const int off_sei = n_ * (m_neg_ + m_pos_) + m_ + m_neg_ + m_pos_ + 2
                      + (lumped ? 1 : 0);
    if (j >= off_sei && j < off_sei + m_neg_) {
        rows.assign(1, j);
    }
}

void DFNModel::refresh_surface_quantities() {
    const auto& solid = params_.solid;
    for (int k = 0; k < m_neg_; ++k) {
        const real_t c_last = c_s_neg_(n_ - 1, k);
        const real_t c_prev = c_s_neg_(n_ - 2, k);
        c_surf_n_(k) = std::clamp(1.5 * c_last - 0.5 * c_prev,
                                  0.0, solid.neg_max_conc);
        u_n_vec_(k) = solid.neg_ocp_at(c_surf_n_(k) / solid.neg_max_conc);
    }
    for (int k = 0; k < m_pos_; ++k) {
        const real_t c_last = c_s_pos_(n_ - 1, k);
        const real_t c_prev = c_s_pos_(n_ - 2, k);
        c_surf_p_(k) = std::clamp(1.5 * c_last - 0.5 * c_prev,
                                  0.0, solid.pos_max_conc);
        u_p_vec_(k) = solid.pos_ocp_at(c_surf_p_(k) / solid.pos_max_conc);
    }
}

void DFNModel::solve_reaction_distribution(real_t i_app) {
    const auto& solid = params_.solid;
    const real_t T = t_cell_; // 集总热模型：Arrhenius 与 BV 系数随电芯温度
    const real_t bv_coef = kFaraday / (2.0 * kGasConst * T);
    const real_t ref = std::max(std::abs(i_app), 0.1);
    const real_t tol_abs = params_.tolerance * ref;

    // 交换电流密度（当前浓度状态冻结；Arrhenius 在 i0 函数内部）
    Vector i0_n(m_neg_);
    Vector i0_p(m_pos_);
    for (int k = 0; k < m_neg_; ++k) {
        i0_n(k) = exchange_current_density_chen2020_neg(
            c_e_(k), c_surf_n_(k), solid.neg_max_conc, T);
    }
    for (int k = 0; k < m_pos_; ++k) {
        i0_p(k) = exchange_current_density_chen2020_pos(
            c_e_(pos_start_ + k), c_surf_p_(k), solid.pos_max_conc, T);
    }

    // 逐格 BV + 自洽 SEI 膜压降（distributed film resistance）：
    //   i = 2·i0·sinh(bv·(a − lr·i))，a = η + δ − j_sei·lr，lr = L·R_sei。
    // 代换 u = a − lr·i（膜后过电位）得 u + c·sinh(bv·u) = a（c = 2·i0·lr）。
    // g(u) 严格递增且 g' ≥ 1、Newton 步长 ≤ 1/bv 有界：从 u0 = a 出发
    // Newton 全局稳健收敛（3~6 步）—— i_surf ↔ η_eff 的强反馈环被
    // 逐格精确解死，外层 Gummel 只剩 φ 场 ↔ i_surf 弱耦合，长循环
    // 大 L（膜压降 V 级）下不再依赖阻尼收敛。返回 {i, di/dδ}
    //（di/dδ 解析：∂F/∂δ / |∂F/∂i|）
    auto bv_film_pair = [&](real_t a, real_t i0_loc,
                            real_t lr) -> std::pair<real_t, real_t> {
        if (lr <= 0.0) {
            // 无膜（降解关闭 / R_sei = 0 / 正极）：直接 BV
            const real_t x = bv_coef * a;
            return {2.0 * i0_loc * safe_sinh(x),
                    2.0 * i0_loc * std::cosh(std::clamp(x, -700.0, 700.0))
                        * bv_coef};
        }
        const real_t c = 2.0 * i0_loc * lr;
        real_t u = a;
        for (int it = 0; it < 40; ++it) {
            const real_t sh = std::sinh(bv_coef * u);
            const real_t g = u + c * sh - a;
            const real_t gp = 1.0 + c * bv_coef * std::cosh(bv_coef * u);
            const real_t step = g / gp;
            if (!std::isfinite(step)) {
                break; // sinh 溢出等病态：保持当前 u，由外层发散检测兜底
            }
            u -= step;
            if (std::abs(step) <= 1e-12 * (1.0 + std::abs(u))) {
                break;
            }
        }
        // di/dδ：a 对 δ 偏导 1；du/da = 1/g'；i = (a − u)/lr
        const real_t duds = 1.0
            / (1.0 + c * bv_coef * std::cosh(bv_coef * u));
        return {(a - u) / lr, (1.0 - duds) / lr};
    };

    // δ 标量 Newton：在当前 η 分布上加偏移 δ，使 Σ a·dx·i(δ) = i_app。
    // η 对 δ 偏导 = 1（φ_s 整体平移），R(δ) 单调增；i(δ) 由 eval 回调
    // 逐格给出（含自洽膜压降的 BV 解与解析导数）；初值取 −mean(η)
    //（δ 已含在 φ 场中，d* 恒为小修正；零电流时精确），Newton 2~4 步收敛
    auto newton_delta = [&](const auto& eval_i, int n_cells, real_t w,
                            real_t i_target, real_t delta0)
                        -> std::pair<real_t, bool> {
        real_t delta = delta0;
        bool ok = false;
        for (int it = 0; it < 25; ++it) {
            real_t R = 0.0;
            real_t dR = 0.0;
            for (int k = 0; k < n_cells; ++k) {
                const auto [ik, dik] = eval_i(k, delta);
                R += w * ik;
                dR += w * dik;
            }
            R -= i_target;
            if (std::abs(R) <= tol_abs) {
                ok = true;
                break;
            }
            if (std::abs(dR) < 1e-30) {
                break;
            }
            const real_t step = R / dR;
            if (!std::isfinite(step)) {
                break; // 病态：保持当前 δ，由外层发散检测兜底
            }
            // 单步限幅：0.3V 下 sinh(5.9)≈180 已远超任何工作电流密度；
            // 更大的步会从 BV 线性区（|η|<25 mV）直接跳进 sinh 指数爆炸区
            // ——静置平衡 η≈0 时线性导数小、满步 2 V 越过整个非线性区，
            // 25 步都退不回来，实测触发整轮 Gummel 发散
            delta -= std::clamp(step, -0.3, 0.3);
        }
        return {delta, ok};
    };

    Vector i_bv_n(m_neg_);
    Vector i_bv_p(m_pos_);
    bool converged = false;

    // 自适应阻尼：SEI 膜电阻反馈环（i_surf → η_eff → i_bv → i_surf）
    // 的环增益随 SEI 厚化（L·R_sei 增大）与温升（i0 Arrhenius 增大）
    // 逐圈增大，固定阻尼的定点迭代在长循环 + 段电流大跳变时振荡
    // 失效（表现为“迭代未收敛”）。残差较历史最优上升视为振荡，
    // 阻尼减半（下限 kMinDamp 保证仍可缓慢收敛）；持续改善则
    // 缓步恢复用户阻尼，兼顾稳健性与收敛速度
    constexpr real_t kMinDamp = 0.05;
    real_t damp = params_.damping;
    real_t best_err = std::numeric_limits<real_t>::infinity();

    const bool deg_on = params_.degradation.enabled;
    const auto& deg = params_.degradation;

    for (int iter = 0; iter < params_.max_iterations; ++iter) {
        // ── 1. 由当前反应场积分电位分布 ──
        integrate_fields(i_app);

        // ── 1.5 SEI 副反应电流（Gummel 不动点：用当前 i_surf/φ/L/T）──
        // j_sei 不动点残差计入收敛判据：滞后一拍的 j_sei 以 ~0.1 mV 级
        // 偏差进入膜压降与 φ 场源，使 Gummel 解不严格满足 dae_residual
        // 的代数行（实测 rb/守恒行 ~1e-3 A/m² 残差），破坏 CN/IDA/
        // DAE 残差三者的方程一致性
        real_t fix_sei = 0.0;
        if (deg_on) {
            for (int k = 0; k < m_neg_; ++k) {
                real_t g_tmp = 0.0;
                real_t q_tmp = 0.0;
                const real_t j_new = sei_kinetics_neg(k, l_sei_(k),
                                                      phi_s_n_(k) - phi_e_(k),
                                                      i_surf_n_(k), T,
                                                      g_tmp, q_tmp);
                fix_sei = std::max(fix_sei, std::abs(j_new - j_sei_n_(k)));
                j_sei_n_(k) = j_new;
                eta_sei_n_(k) = (phi_s_n_(k) - phi_e_(k)) - deg.u_sei
                              - i_surf_n_(k) * l_sei_(k) * deg.r_sei;
            }
        }

        // ── 2. 局部 BV（当前 φ 场）──
        for (int k = 0; k < m_neg_; ++k) {
            eta_n_(k) = phi_s_n_(k) - phi_e_(k) - u_n_vec_(k);
        }
        for (int k = 0; k < m_pos_; ++k) {
            eta_p_(k) = phi_s_p_(k) - phi_e_(pos_start_ + k) - u_p_vec_(k);
        }

        // ── 3. δ Newton（负极目标 = i_app − SEI 副反应分流；正极 −i_app）──
        // 负极逐格 BV 含自洽 SEI 膜压降 −(i + j_sei)·L·R_sei
        //（PyBaMM base_kinetics distributed 口径），i_surf ↔ η_eff
        // 强反馈环由逐格精确解解死（见 bv_film_pair 注释）
        const real_t j_sei_sum = deg_on ? a_n_ * dx_ * j_sei_n_.sum() : 0.0;
        auto eval_n = [&](int k, real_t d) {
            const real_t lr = deg_on ? l_sei_(k) * deg.r_sei : 0.0;
            const real_t a = eta_n_(k) + d
                           - (deg_on ? j_sei_n_(k) * lr : 0.0);
            return bv_film_pair(a, i0_n(k), lr);
        };
        auto eval_p = [&](int k, real_t d) {
            return bv_film_pair(eta_p_(k) + d, i0_p(k), 0.0); // 正极无 SEI 膜
        };
        const auto [d_n, ok_n] = newton_delta(eval_n, m_neg_, a_n_ * dx_,
                                              i_app - j_sei_sum,
                                              -eta_n_.mean());
        if (!ok_n) {
            // δ 反解失败：不更新 δ 基准/反应场（污染会使后续轮 φ 场
            // 崩坏成 inf/NaN），交给外层未收敛失败路径
            break;
        }
        delta_n_ += d_n;
        const auto [d_p, ok_p] = newton_delta(eval_p, m_pos_, a_p_ * dx_, -i_app,
                                              -eta_p_.mean());
        if (!ok_p) {
            break;
        }
        delta_p_ += d_p;

        // ── 4. 最终 δ 下的 BV 电流、不动点残差与电流守恒残差 ──
        // 收敛要求两者同时满足：仅守恒（首轮 Newton 即精确）不足以保证
        // φ 场与 j 场自洽，需 |i_bv − i_surf| → 0（Gummel 不动点）。
        real_t R_n = 0.0;
        real_t R_p = 0.0;
        real_t fix_n = 0.0;
        real_t fix_p = 0.0;
        for (int k = 0; k < m_neg_; ++k) {
            i_bv_n(k) = eval_n(k, d_n).first;
            R_n += a_n_ * dx_ * i_bv_n(k);
            fix_n = std::max(fix_n, std::abs(i_bv_n(k) - i_surf_n_(k)));
        }
        for (int k = 0; k < m_pos_; ++k) {
            i_bv_p(k) = eval_p(k, d_p).first;
            R_p += a_p_ * dx_ * i_bv_p(k);
            fix_p = std::max(fix_p, std::abs(i_bv_p(k) - i_surf_p_(k)));
        }
        R_n -= i_app - j_sei_sum; // 主反应守恒目标含 SEI 副反应分流
        R_p += i_app;

        // 电流守恒（R）与 φ–j 自洽（fix）同时满足；分母含电流分布幅值，
        // 避免静置（i_app → 0，残余分布电流 ≫ 端口电流）时判据过严。
        // j_sei 残差按自身幅值归一（并进主 err 会因 i_surf 量级差异
        // 被淹没）：直接驱动膜压降与寿命量，需同量级相对精度
        const real_t i_max = std::max(i_bv_n.cwiseAbs().maxCoeff(),
                                      i_bv_p.cwiseAbs().maxCoeff());
        const real_t err = std::max({std::abs(R_n), std::abs(R_p),
                                     fix_n, fix_p})
                         / std::max(ref, i_max);
        const real_t err_kin = deg_on
            ? fix_sei / std::max(j_sei_n_.cwiseAbs().maxCoeff(), 1e-9)
            : 0.0;
        const real_t err_total = std::max(err, err_kin);
        // 数值防护：BV 电流溢出（sinh 爆炸）立即判失败，避免 NaN 传入状态
        if (!i_bv_n.allFinite() || !i_bv_p.allFinite()) {
            throw std::runtime_error("DFN: 反应分布迭代发散（BV 电流溢出，可减小时间步长）");
        }
        if (err_total < params_.tolerance) {
            // 收敛：取完整 BV 解（不用阻尼），保证 φ 场与 j 场自洽
            i_surf_n_ = i_bv_n;
            i_surf_p_ = i_bv_p;
            for (int k = 0; k < m_neg_; ++k) {
                eta_n_(k) += d_n;
            }
            for (int k = 0; k < m_pos_; ++k) {
                eta_p_(k) += d_p;
            }
            converged = true;
            break;
        }

        // 自适应阻尼更新（见上方说明；升降判定基于历史最优残差，
        // 单次抖动不触发阻尼回撤）
        if (err < best_err) {
            best_err = err;
            damp = std::min(params_.damping, damp * 1.25);
        } else {
            damp = std::max(kMinDamp, damp * 0.5);
        }
        i_surf_n_ += damp * (i_bv_n - i_surf_n_);
        i_surf_p_ += damp * (i_bv_p - i_surf_p_);
    }

    if (!converged) {
        throw std::runtime_error("DFN: 反应分布迭代未收敛（可减小时间步长或电流）");
    }

    // 最终自洽积分：用收敛的 i_surf 重积分 φ 场并输出电压
    finalize_fields(i_app);
}

void DFNModel::finalize_fields(real_t i_app) {
    integrate_fields(i_app);
    for (int k = 0; k < m_neg_; ++k) {
        eta_n_(k) = phi_s_n_(k) - phi_e_(k) - u_n_vec_(k);
    }
    for (int k = 0; k < m_pos_; ++k) {
        eta_p_(k) = phi_s_p_(k) - phi_e_(pos_start_ + k) - u_p_vec_(k);
    }
    eta_e_ = phi_e_(m_ - 1);
    voltage_ = phi_s_p_(m_pos_ - 1) - delta_n_;
}

void DFNModel::integrate_fields(real_t i_app) {
    // 液相面电流密度：∂i_e/∂x = a·(i_surf + j_sei)，两端集流体处为 0
    //（SEI 副反应的电子也经固相/液相传导，电流守恒下 i_e 末端仍 = i_app）
    const bool deg_on = params_.degradation.enabled;
    i_e_face_(0) = 0.0;
    for (int i = 0; i < m_; ++i) {
        real_t src = 0.0;
        if (i < sep_start_) {
            src = a_n_ * (i_surf_n_(i) + (deg_on ? j_sei_n_(i) : 0.0)) * dx_;
        } else if (i >= pos_start_) {
            src = a_p_ * i_surf_p_(i - pos_start_) * dx_;
        }
        i_e_face_(i + 1) = i_e_face_(i) + src;
    }

    // 液相电位：∂φ_e/∂x = −i_e/κ^eff + (2RT/F)(1−t₊)·∂ln c_e/∂x，φ_e(0) = 0
    //（κ^eff = κ(c_e)·ε^bruggeman 浓度相关，逐格求值；RT 取电芯温度）
    const real_t conc_coef = 2.0 * kGasConst * t_cell_ / kFaraday
                           * (1.0 - params_.transference_number);
    phi_e_(0) = 0.0;
    for (int i = 1; i < m_; ++i) {
        const real_t kappa_face = 0.5 * (kappa_eff_at(c_e_, i - 1) + kappa_eff_at(c_e_, i));
        phi_e_(i) = phi_e_(i - 1)
                  - i_e_face_(i) / kappa_face * dx_
                  + conc_coef * (std::log(c_e_(i)) - std::log(c_e_(i - 1)));
    }

    // 固相电位：∂φ_s/∂x = −i_s/σ^eff，i_s = i_app − i_e（总电流守恒）
    // 负极基准 φ_s_n(0) = δ_n；正极基准 φ_s_p(0) = δ_p（界面处）
    phi_s_n_(0) = delta_n_;
    for (int i = 1; i < m_neg_; ++i) {
        phi_s_n_(i) = phi_s_n_(i - 1)
                    - (i_app - i_e_face_(i)) / sigma_n_eff_ * dx_;
    }
    phi_s_p_(0) = delta_p_;
    for (int i = 1; i < m_pos_; ++i) {
        phi_s_p_(i) = phi_s_p_(i - 1)
                    - (i_app - i_e_face_(pos_start_ + i)) / sigma_p_eff_ * dx_;
    }

    // 全厚度固相电位（隔膜区无固相，取负极末端延伸，仅供绘图/快照）
    phi_s_full_.segment(0, m_neg_) = phi_s_n_;
    phi_s_full_.segment(sep_start_, pos_start_ - sep_start_)
        .setConstant(phi_s_n_(m_neg_ - 1));
    phi_s_full_.segment(pos_start_, m_pos_) = phi_s_p_;
}

void DFNModel::step(real_t current, real_t dt) {
    if (!(dt > 0.0)) {
        throw std::invalid_argument("DFNModel::step: dt must be positive");
    }
    // CN 表面源用步前 i_surf 冻结值（显式处理），而表面浓度–BV 局部回路
    // 的弛豫时间尺度 τ 可远小于 dt（OCP 陡峭段 + 静置大步长最典型），
    // dt > 2τ 时显式源发散振荡 → 按 plan_substeps 细分保证子步内衰减。
    const int sub_steps = plan_substeps(dt);
    const real_t dt_sub = dt / static_cast<real_t>(sub_steps);
    for (int k = 0; k < sub_steps; ++k) {
        step_once(current, dt_sub);
    }
}

void DFNModel::step_once(real_t current, real_t dt) {
    ensure_cn(dt);

    // 端口电流 → 面电流密度 [A/m²]（反应分布求解内部统一用面电流语义）
    const real_t i_app = current / area_;
    const bool lumped = params_.thermal_type == ThermalType::Lumped;

    // 1. 用当前浓度状态解代数反应分布（本步的源项与初值）
    solve_reaction_distribution(i_app);

    // 步前产热（Q 滞后显式，与液相 D_e 滞后显式风格一致；等温仅作诊断缓存）
    const real_t q_old = refresh_heat(i_app);

    // 2. 固相推进：各列独立 CN，表面源 = −A_surf/V_surf·i_surf/F（逐列）
    const real_t s_n_coef = -dt * surf_area_n_ / surf_vol_n_ / kFaraday;
    Eigen::MatrixXd S_n = Eigen::MatrixXd::Zero(n_, m_neg_);
    S_n.row(n_ - 1) = s_n_coef * i_surf_n_.transpose();
    Eigen::MatrixXd B_n = cn_plus_neg_ * c_s_neg_ + S_n;
    c_s_neg_ = cn_solve_neg_.solve(B_n);
    c_s_neg_ = c_s_neg_.cwiseMax(0.0).cwiseMin(params_.solid.neg_max_conc);

    const real_t s_p_coef = -dt * surf_area_p_ / surf_vol_p_ / kFaraday;
    Eigen::MatrixXd S_p = Eigen::MatrixXd::Zero(n_, m_pos_);
    S_p.row(n_ - 1) = s_p_coef * i_surf_p_.transpose();
    Eigen::MatrixXd B_p = cn_plus_pos_ * c_s_pos_ + S_p;
    c_s_pos_ = cn_solve_pos_.solve(B_p);
    c_s_pos_ = c_s_pos_.cwiseMax(0.0).cwiseMin(params_.solid.pos_max_conc);

    // 3. 液相推进：逐格源 (1−t₊)/F·a·(i_surf + j_sei)
    //（SEI 消耗电解液锂离子，Yang2017 口径；固相边界通量仅主反应）
    const bool deg_on = params_.degradation.enabled;
    Vector s_vec = Vector::Zero(m_);
    for (int k = 0; k < m_neg_; ++k) {
        s_vec(k) = (1.0 - params_.transference_number) / kFaraday
                 * a_n_ * (i_surf_n_(k) + (deg_on ? j_sei_n_(k) : 0.0));
    }
    for (int k = 0; k < m_pos_; ++k) {
        s_vec(pos_start_ + k) = (1.0 - params_.transference_number) / kFaraday
                               * a_p_ * i_surf_p_(k);
    }
    Vector B_e = le_plus_ * c_e_ + dt * s_vec;
    c_e_ = le_solve_.solve(B_e);
    if (le_solve_.info() != Eigen::Success) {
        throw std::runtime_error("DFNModel: electrolyte CN solve failed");
    }
    c_e_ = c_e_.cwiseMax(kCeFloor);

    // 4. 刷新表面量，并用新状态重解代数系统（快照电压与最终状态一致）
    refresh_surface_quantities();
    solve_reaction_distribution(i_app);

    // 4.5 SEI 厚度推进（退化启用）：步末自洽 φ/T/i_surf 冻结，隐式推进；
    // SEI 生长慢（步内变化 ≪ L），下一步反应分布迭代会以新厚度自洽
    if (params_.degradation.enabled) {
        advance_sei_cn(dt);
    }

    // 5. 步末产热 + 集总温度更新（Crank–Nicolson）：
    // T1·(1 + α) = T0·(1 − α) + dt·(q_old + q_new)/2/C_th + dt·hA·T_amb/C_th，
    // α = dt·hA/(2·C_th)（散热半隐式，无条件稳定；产热梯形平均保二阶）
    const real_t q_new = refresh_heat(i_app);
    if (lumped) {
        const auto& th = params_.thermal_params;
        const real_t c_th = thermal_capacitance();
        const real_t cool = th.h_total * th.cooling_area;
        const real_t half_cool_dt = 0.5 * dt * cool / c_th;
        t_cell_ = (t_cell_ * (1.0 - half_cool_dt)
                   + dt * (0.5 * (q_old + q_new) + cool * th.t_ambient) / c_th)
                 / (1.0 + half_cool_dt);
    }

    time_ += dt;
}

int DFNModel::plan_substeps(real_t dt) const {
    const auto& solid = params_.solid;
    const real_t T = t_cell_;
    const real_t bv_coef = kFaraday / (2.0 * kGasConst * T);
    real_t lambda_max = 0.0;

    // 局部回路：dc_surf/dt = −(A/V)·i_bv(c_surf)/F，i_bv 经 i0(c) 与 U(c)
    // 依赖 c_surf（φ 场与 η 在步内冻结）。对 i_bv(c) 中心差分得
    // ∂i_bv/∂c_surf，λ = (A/V)·max(∂i_bv/∂c_surf, 0)/F 为回路弛豫率；
    // 静置（η ≈ 0）时刚性仅来自 OCP 陡峭段（|dU/dc| 大）。
    auto scan = [&](const Vector& c_surf, const Vector& eta, real_t c_max,
                    real_t av, int c_e_offset,
                    real_t (*i0_fn)(real_t, real_t, real_t, real_t),
                    real_t (*ocp_fn)(real_t)) {
        const real_t h = 1e-4 * c_max;
        for (int k = 0; k < c_surf.size(); ++k) {
            const real_t c = c_surf(k);
            const real_t c_e = c_e_(k + c_e_offset);
            const real_t drive = eta(k) + ocp_fn(c / c_max);
            auto bv_at = [&](real_t cc) {
                const real_t i0 = i0_fn(c_e, cc, c_max, T);
                return 2.0 * i0 * std::sinh(bv_coef * (drive - ocp_fn(cc / c_max)));
            };
            const real_t di_dc = (bv_at(c + h) - bv_at(c - h)) / (2.0 * h);
            lambda_max = std::max(lambda_max, av * std::max(di_dc, 0.0) / kFaraday);
        }
    };
    scan(c_surf_n_, eta_n_, solid.neg_max_conc, surf_area_n_ / surf_vol_n_, 0,
         exchange_current_density_chen2020_neg, graphite_chen2020_ocp);
    scan(c_surf_p_, eta_p_, solid.pos_max_conc, surf_area_p_ / surf_vol_p_,
         pos_start_, exchange_current_density_chen2020_pos, nmc_chen2020_ocp);

    if (lambda_max <= 0.0) {
        return 1;
    }
    // 子步长 ≤ 0.4/λ → 显式源放大因子 |1 − dt_sub·λ| ≤ 0.6（衰减）
    const int sub = static_cast<int>(std::ceil(dt * lambda_max / 0.4));
    return std::clamp(sub, 1, 1000);
}

real_t DFNModel::rho_c_p_eff() const {
    // PyBaMM rho_c_p_eff：各层（含集流体）ρc_p 按厚度加权
    const auto& th = params_.thermal_params;
    const auto& solid = params_.solid;
    const real_t l_total = solid.neg_thickness + solid.sep_thickness
                         + solid.pos_thickness
                         + th.cc_thickness_neg + th.cc_thickness_pos;
    return (th.rho_c_p_cc_neg * th.cc_thickness_neg
            + th.rho_c_p_neg * solid.neg_thickness
            + th.rho_c_p_sep * solid.sep_thickness
            + th.rho_c_p_pos * solid.pos_thickness
            + th.rho_c_p_cc_pos * th.cc_thickness_pos) / l_total;
}

real_t DFNModel::thermal_capacitance() const {
    return rho_c_p_eff_ * params_.thermal_params.cell_volume;
}

real_t DFNModel::refresh_heat(real_t i_app) {
    // CN 路径：成员场已是自洽解（调用点在 solve_reaction_distribution 之后）
    HeatSourceArgs a;
    a.m = m_;
    a.m_neg = m_neg_;
    a.pos_start = pos_start_;
    a.dx = dx_;
    a.area = area_;
    a.a_n = a_n_;
    a.a_p = a_p_;
    a.sigma_n_eff = sigma_n_eff_;
    a.sigma_p_eff = sigma_p_eff_;
    a.i_app = i_app;
    a.T = t_cell_;
    a.entropic_n = params_.thermal_params.entropic_change_neg;
    a.entropic_p = params_.thermal_params.entropic_change_pos;
    a.lut_n = &params_.thermal_params.entropic_lut_neg;
    a.lut_p = &params_.thermal_params.entropic_lut_pos;
    a.c_surf_n = &c_surf_n_; // refresh_surface_quantities 已刷新
    a.c_surf_p = &c_surf_p_;
    a.cs_max_n = params_.solid.neg_max_conc;
    a.cs_max_p = params_.solid.pos_max_conc;
    a.i_e_face = &i_e_face_;
    a.phi_e = &phi_e_;
    a.i_surf_n = &i_surf_n_;
    a.i_surf_p = &i_surf_p_;
    a.eta_n = &eta_n_;
    a.eta_p = &eta_p_;
    if (params_.degradation.enabled) {
        a.j_sei_n = &j_sei_n_;
        a.eta_sei_n = &eta_sei_n_;
    }
    compute_heat_source(a, q_ohm_, q_rxn_, q_rev_);
    q_total_ = q_ohm_ + q_rxn_ + q_rev_;
    return q_total_;
}

real_t DFNModel::total_lithium_mol() const {
    const auto& solid = params_.solid;

    // 固相：每列颗粒相同几何（格体积 A·dx·ε_s），颗粒数按列一致
    const real_t particle_vol_n = (4.0 / 3.0) * kPi
                                * solid.neg_particle_radius * solid.neg_particle_radius
                                * solid.neg_particle_radius;
    const real_t n_per_col_n = area_ * dx_ * solid.neg_active_fraction / particle_vol_n;
    const real_t particle_vol_p = (4.0 / 3.0) * kPi
                                * solid.pos_particle_radius * solid.pos_particle_radius
                                * solid.pos_particle_radius;
    const real_t n_per_col_p = area_ * dx_ * solid.pos_active_fraction / particle_vol_p;

    real_t mol_n = 0.0;
    for (int k = 0; k < m_neg_; ++k) {
        for (int i = 0; i < n_; ++i) {
            mol_n += vol_neg_[i] * c_s_neg_(i, k);
        }
    }
    real_t mol_p = 0.0;
    for (int k = 0; k < m_pos_; ++k) {
        for (int i = 0; i < n_; ++i) {
            mol_p += vol_pos_[i] * c_s_pos_(i, k);
        }
    }

    // 液相：Σ ε·c·A·dx（零通量 + 源守恒 → 恒定）
    real_t liquid_mol = 0.0;
    for (int i = 0; i < m_; ++i) {
        liquid_mol += eps_e_(i) * c_e_(i) * dx_;
    }

    return n_per_col_n * mol_n + n_per_col_p * mol_p + area_ * liquid_mol;
}

SolutionSnapshot DFNModel::snapshot(real_t current) const {
    SolutionSnapshot snap;
    snap.time = time_;
    snap.voltage = voltage_;
    snap.current = current;
    snap.temperature = t_cell_;
    snap.c_e = c_e_;
    snap.phi_e = phi_e_;
    snap.phi_s = phi_s_full_;
    snap.sep_start = sep_start_;
    snap.pos_start = pos_start_;
    snap.c_s_surf_neg = Vector::Constant(1, c_surf_n_.mean());
    snap.c_s_surf_pos = Vector::Constant(1, c_surf_p_.mean());
    snap.c_s_neg = c_s_neg_.rowwise().mean(); // 平均颗粒径向分布（GUI 兼容）
    snap.c_s_pos = c_s_pos_.rowwise().mean();
    snap.x_surf = c_surf_n_.mean() / params_.solid.neg_max_conc;
    snap.y_surf = c_surf_p_.mean() / params_.solid.pos_max_conc;
    snap.eta_neg = eta_n_.mean();
    snap.eta_pos = eta_p_.mean();
    snap.eta_e = eta_e_;
    snap.u_neg = u_n_vec_.mean();
    snap.u_pos = u_p_vec_.mean();

    // DFN 分布输出
    snap.c_surf_neg_x = c_surf_n_;
    snap.c_surf_pos_x = c_surf_p_;
    snap.i_surf_neg = i_surf_n_;
    snap.i_surf_pos = i_surf_p_;
    snap.eta_neg_x = eta_n_;
    snap.eta_pos_x = eta_p_;
    snap.u_neg_x = u_n_vec_;
    snap.u_pos_x = u_p_vec_;
    snap.i_e_face = i_e_face_;

    // 产热诊断（step/构造后缓存与当前状态一致）
    snap.q_ohm = q_ohm_;
    snap.q_rxn = q_rxn_;
    snap.q_rev = q_rev_;
    snap.q_total = q_total_;

    // 机理寿命诊断（退化未启用时恒 0）
    snap.li_loss_ah = lithium_loss_ah();
    snap.degraded_capacity_ah = degraded_capacity_ah();
    snap.sei_thickness = sei_thickness_avg();
    return snap;
}

// ══════════════════════════════════════════════════════════════
//  DAE 视图（SUNDIALS IDA）：指数-1 DAE
//  微分变量：c_s（每列颗粒）、c_e；代数变量：i_surf、δ_n、δ_p
//  代数约束：逐格 BV 方程（φ 场由 i_surf/δ 解析积分）+ 电流守恒
// ══════════════════════════════════════════════════════════════

int DFNModel::dae_size() const {
    // lumped 热模型在 δ 后追加温度状态（微分变量）；
    // SEI 降解在温度后追加 m_neg 个 SEI 厚度状态（微分变量）；
    // 端口电流 I_cell 恒在末位（增广外电路端口的代数变量）
    const int extra = (params_.thermal_type == ThermalType::Lumped ? 1 : 0)
                    + (params_.degradation.enabled ? m_neg_ : 0)
                    + 1;
    return n_ * (m_neg_ + m_pos_) + m_ + m_neg_ + m_pos_ + 2 + extra;
}

void DFNModel::dae_mask(Vector& id) const {
    const int off = n_ * (m_neg_ + m_pos_) + m_; // 浓度场之后全为代数变量
    const bool lumped = params_.thermal_type == ThermalType::Lumped;
    const bool deg_on = params_.degradation.enabled;
    const int off_d = n_ * (m_neg_ + m_pos_) + m_ + m_neg_ + m_pos_;
    id.resize(dae_size());
    id.segment(0, off).setOnes();
    id.segment(off, dae_size() - off).setZero();
    if (lumped) {
        id(off_d + 2) = 1.0; // 温度为微分变量
    }
    if (deg_on) {
        id.segment(off_d + 2 + (lumped ? 1 : 0), m_neg_).setOnes();
    }
}

void DFNModel::dae_atol(Vector& atol) const {
    const int s_n = n_ * m_neg_;
    const int s_p = n_ * m_pos_;
    const int off_e = s_n + s_p;
    const int off_isn = off_e + m_;
    const int off_isp = off_isn + m_neg_;
    const int off_d = off_isp + m_pos_;
    const bool lumped = params_.thermal_type == ThermalType::Lumped;
    const int off_sei = off_d + 2 + (lumped ? 1 : 0);

    atol.resize(dae_size());
    atol.segment(0, s_n).setConstant(0.5);       // c_s_n [mol/m³]
    atol.segment(s_n, s_p).setConstant(1.0);     // c_s_p [mol/m³]
    atol.segment(off_e, m_).setConstant(0.05);   // c_e [mol/m³]
    atol.segment(off_isn, m_neg_ + m_pos_).setConstant(0.005); // i_surf [A/m²]
    atol(off_d) = 1e-5;                          // δ_n [V]
    atol(off_d + 1) = 1e-5;                      // δ_p [V]
    if (params_.degradation.enabled) {
        atol.segment(off_sei, m_neg_).setConstant(1e-12); // L_sei [m]
    }
    if (lumped) {
        atol(off_d + 2) = 1e-4;                  // T_cell [K]
    }
    atol(dae_size() - 1) = 0.005;                // I_cell [A]（仅影响差分
                                                 // σ/CalcIC；SuppressAlg
                                                 // 下不参与误差控制）
}

void DFNModel::dae_pack_state(Vector& y, Vector& yp) const {
    const int s_n = n_ * m_neg_;
    const int s_p = n_ * m_pos_;
    const int off_e = s_n + s_p;
    const int off_isn = off_e + m_;
    const int off_isp = off_isn + m_neg_;
    const int off_d = off_isp + m_pos_;
    const bool lumped = params_.thermal_type == ThermalType::Lumped;
    const bool deg_on = params_.degradation.enabled;
    const int off_sei = off_d + 2 + (lumped ? 1 : 0);

    y.resize(dae_size());
    yp.resize(dae_size());
    Eigen::Map<Eigen::MatrixXd>(y.data(), n_, m_neg_) = c_s_neg_;
    Eigen::Map<Eigen::MatrixXd>(y.data() + s_n, n_, m_pos_) = c_s_pos_;
    y.segment(off_e, m_) = c_e_;
    y.segment(off_isn, m_neg_) = i_surf_n_;
    y.segment(off_isp, m_pos_) = i_surf_p_;
    y(off_d) = delta_n_;
    y(off_d + 1) = delta_p_;
    y(dae_size() - 1) = i_port_; // 增广端口电流（代数变量）

    // 导数：浓度场 = f(y)（代数变量导数为 0，IDA 只使用微分变量导数）
    Eigen::Map<Eigen::MatrixXd>(yp.data(), n_, m_neg_) =
        D_neg_ * c_s_neg_;
    Eigen::Map<Eigen::MatrixXd>(yp.data(), n_, m_neg_).row(n_ - 1) -=
        (surf_area_n_ / surf_vol_n_ / kFaraday) * i_surf_n_.transpose();
    Eigen::Map<Eigen::MatrixXd>(yp.data() + s_n, n_, m_pos_) =
        D_pos_ * c_s_pos_;
    Eigen::Map<Eigen::MatrixXd>(yp.data() + s_n, n_, m_pos_).row(n_ - 1) -=
        (surf_area_p_ / surf_vol_p_ / kFaraday) * i_surf_p_.transpose();
    Vector src = Vector::Zero(m_);
    // 液相源项与 dae_residual 同口径：负极含 SEI 副反应 j_sei 的电解液
    // 锂盐消耗（Yang2017）。pack 的 yp 是 IDA CalcIC 的一致化初值，
    // 漏 j_sei 会使液相行初值残差 ~(1−t₊)·a·j_sei/F（1C 厚膜 ~1.7
    // mol/m³/s，实测负极区 49 格逐格非零），每次 reinit/CalcIC 都
    // 从带病初值出发
    for (int k = 0; k < m_neg_; ++k) {
        src(k) = (1.0 - params_.transference_number) / kFaraday
               * a_n_ * (i_surf_n_(k) + (deg_on ? j_sei_n_(k) : 0.0));
    }
    for (int k = 0; k < m_pos_; ++k) {
        src(pos_start_ + k) = (1.0 - params_.transference_number) / kFaraday
                            * a_p_ * i_surf_p_(k);
    }
    yp.segment(off_e, m_) =
        (electrolyte_div_flux(c_e_) + src).cwiseQuotient(eps_e_);
    yp.segment(off_isn, m_neg_ + m_pos_ + 2).setZero();
    yp(dae_size() - 1) = 0.0; // I_cell：代数变量导数恒零
    if (deg_on) {
        // SEI 厚度状态与导数初值（IDA CalcIC 会进一步一致化）
        y.segment(off_sei, m_neg_) = l_sei_;
        const real_t dldt_coef = -params_.degradation.v_bar_sei
                               / (kFaraday * params_.degradation.z_sei);
        for (int k = 0; k < m_neg_; ++k) {
            yp(off_sei + k) = dldt_coef * j_sei_n_(k);
        }
    }
    if (lumped) {
        // 温度 yp 初值：用最近产热缓存估计（IDA CalcIC 会进一步一致化）
        const auto& th = params_.thermal_params;
        const real_t c_th = rho_c_p_eff_ * th.cell_volume;
        const real_t cool = th.h_total * th.cooling_area;
        y(off_d + 2) = t_cell_;
        yp(off_d + 2) = (q_total_ - cool * (t_cell_ - th.t_ambient)) / c_th;
    }
}

void DFNModel::dae_residual(real_t t, const Eigen::Ref<const Vector>& y,
                                const Eigen::Ref<const Vector>& yp,
                                Eigen::Ref<Vector> res) {
    (void)t;
    const auto& solid = params_.solid;

    const int s_n = n_ * m_neg_;
    const int s_p = n_ * m_pos_;
    const int off_e = s_n + s_p;
    const int off_isn = off_e + m_;
    const int off_isp = off_isn + m_neg_;
    const int off_d = off_isp + m_pos_;
    const int off_t = off_d + 2; // lumped 温度状态（微分变量）
    const bool lumped = params_.thermal_type == ThermalType::Lumped;
    const bool deg_on = params_.degradation.enabled;
    const auto& deg = params_.degradation;
    const int off_sei = off_d + 2 + (lumped ? 1 : 0);

    // lumped：BV 系数与液相浓差项用电芯温度（状态变量），等温用参数温度
    const real_t T = lumped ? y(off_t) : solid.temperature;
    const real_t bv_coef = kFaraday / (2.0 * kGasConst * T);
    // 端口电流来自状态（增广外电路端口）：内部方程（固相欧姆、守恒行、
    // 热行）对外电路电流的敏感度由此进入 Jacobian，恒压行才能由
    // Newton 与内部方程一致求解（pybamm.external_circuit 同构）
    const real_t i_cell = y(dae_size() - 1);
    const real_t i_app = i_cell / area_;

    Eigen::Map<const Eigen::MatrixXd> c_n(y.data(), n_, m_neg_);
    Eigen::Map<const Eigen::MatrixXd> c_p(y.data() + s_n, n_, m_pos_);
    Eigen::Map<const Vector> c_e(y.data() + off_e, m_);
    Eigen::Map<const Vector> is_n(y.data() + off_isn, m_neg_);
    Eigen::Map<const Vector> is_p(y.data() + off_isp, m_pos_);
    const real_t delta_n = y(off_d);
    const real_t delta_p = y(off_d + 1);

    Eigen::Map<const Eigen::MatrixXd> cd_n(yp.data(), n_, m_neg_);
    Eigen::Map<const Eigen::MatrixXd> cd_p(yp.data() + s_n, n_, m_pos_);
    Eigen::Map<const Vector> cd_e(yp.data() + off_e, m_);

    Eigen::Map<Eigen::MatrixXd> r_n(res.data(), n_, m_neg_);
    Eigen::Map<Eigen::MatrixXd> r_p(res.data() + s_n, n_, m_pos_);
    Eigen::Map<Vector> r_e(res.data() + off_e, m_);
    Eigen::Map<Vector> rb_n(res.data() + off_isn, m_neg_);
    Eigen::Map<Vector> rb_p(res.data() + off_isp, m_pos_);

    // SEI 降解：局部状态向量（L 来自 y 而非成员；j_sei/eta_sei 由残差
    // 自洽计算，供 BV 行/守恒行/液相源/SEI 行/产热共用）
    Vector j_sei_ln = Vector::Zero(m_neg_);
    Vector eta_sei_ln = Vector::Zero(m_neg_);
    Vector l_sei_y = Vector::Zero(m_neg_);
    if (deg_on) {
        l_sei_y = Eigen::Map<const Vector>(y.data() + off_sei, m_neg_);
    }

    // ── 微分行（固相）：c' = D·c − (A_surf/V_surf)·i_surf/F
    //（SEI 副反应锂来自电解液，固相边界通量不含 j_sei）──
    r_n = cd_n - (D_neg_ * c_n);
    r_n.row(n_ - 1) += (surf_area_n_ / surf_vol_n_ / kFaraday)
                       * is_n.transpose();
    r_p = cd_p - (D_pos_ * c_p);
    r_p.row(n_ - 1) += (surf_area_p_ / surf_vol_p_ / kFaraday)
                       * is_p.transpose();

    // ── 代数行：φ 场由 i_surf/j_sei/δ 解析积分 ──
    // j_sei 依赖 φ（eta_SEI 的膜电阻项）而 φ 场又含 j_sei 的体电流，
    // 弱耦合回路（增益 ≪ 1）：不动点迭代 2~4 轮收敛；初值用成员
    // j_sei_n_（残差连续调用间热启动）
    Vector i_e(m_ + 1), phi_e(m_), phi_sn(m_neg_), phi_sp(m_pos_);
    auto integrate_phi = [&](const Vector& j_sei_loc) {
        i_e(0) = 0.0;
        for (int i = 0; i < m_; ++i) {
            real_t face_src = 0.0;
            if (i < sep_start_) {
                face_src = a_n_ * (is_n(i) + j_sei_loc(i)) * dx_;
            } else if (i >= pos_start_) {
                face_src = a_p_ * is_p(i - pos_start_) * dx_;
            }
            i_e(i + 1) = i_e(i) + face_src;
        }
        const real_t conc_coef = 2.0 * kGasConst * T / kFaraday
                               * (1.0 - params_.transference_number);
        phi_e(0) = 0.0;
        for (int i = 1; i < m_; ++i) {
            const real_t kappa_face =
                0.5 * (kappa_eff_at(c_e, i - 1) + kappa_eff_at(c_e, i));
            phi_e(i) = phi_e(i - 1)
                     - i_e(i) / kappa_face * dx_
                     + conc_coef * (std::log(c_e(i)) - std::log(c_e(i - 1)));
        }
        phi_sn(0) = delta_n;
        for (int i = 1; i < m_neg_; ++i) {
            phi_sn(i) = phi_sn(i - 1)
                      - (i_app - i_e(i)) / sigma_n_eff_ * dx_;
        }
        phi_sp(0) = delta_p;
        for (int i = 1; i < m_pos_; ++i) {
            phi_sp(i) = phi_sp(i - 1)
                      - (i_app - i_e(pos_start_ + i)) / sigma_p_eff_ * dx_;
        }
    };

    Vector j_sei_it = deg_on ? Vector(j_sei_n_) : Vector::Zero(m_neg_);
    integrate_phi(j_sei_it);
    if (deg_on) {
        // 膜环（j_sei ↔ φ 场）负反馈回路的不动点迭代。迭代收敛后结果
        // 与初值无关（热启动仅加速）；但温度升高经 Arrhenius 放大 K、
        // 环增益变强后固定轮数会截断在中途 —— 同一状态两次求值
        //（残差/差分 Jacobian 列）得到不同结果，corrector 的残差
        // 一致性被破坏。改为自适应阻尼迭代到机器级收敛，保证残差是
        // y 的确定函数。
        real_t omega = 1.0;
        real_t last_err = std::numeric_limits<real_t>::infinity();
        Vector j_new = Vector::Zero(m_neg_);
        for (int it = 0; it < 48; ++it) {
            // err 每轮重置：本轮修正量范数才是收敛判据（若跨轮累计
            // 历史最大值，首轮后的 1e-12 判据永假，烧满轮数上限）
            real_t err = 0.0;
            for (int k = 0; k < m_neg_; ++k) {
                real_t g_tmp = 0.0;
                real_t q_tmp = 0.0;
                j_new(k) = sei_kinetics_neg(k, l_sei_y(k),
                                            phi_sn(k) - phi_e(k),
                                            is_n(k), T, g_tmp, q_tmp);
                err = std::max(err, std::abs(j_new(k) - j_sei_it(k)));
            }
            // 振荡检测：环误差不降则加强欠松弛（阻尼负反馈指数环）
            if (err > last_err) {
                omega = std::max(omega * 0.5, 0.05);
            }
            last_err = err;
            j_sei_it += omega * (j_new - j_sei_it);
            if (err <= 1e-12) {
                break;
            }
            integrate_phi(j_sei_it);
        }
        for (int k = 0; k < m_neg_; ++k) {
            j_sei_ln(k) = j_sei_it(k);
            eta_sei_ln(k) = (phi_sn(k) - phi_e(k)) - deg.u_sei
                          - is_n(k) * l_sei_y(k) * deg.r_sei;
        }
    }

    // 液相微分行：ε c' = L(c)·c + (1−t₊)/F·a·(i_surf + j_sei)
    //（SEI 消耗电解液锂离子，Yang2017 口径）
    Vector src = Vector::Zero(m_);
    for (int k = 0; k < m_neg_; ++k) {
        src(k) = (1.0 - params_.transference_number) / kFaraday
               * a_n_ * (is_n(k) + j_sei_ln(k));
    }
    for (int k = 0; k < m_pos_; ++k) {
        src(pos_start_ + k) = (1.0 - params_.transference_number) / kFaraday
                            * a_p_ * is_p(k);
    }
    r_e = eps_e_.cwiseProduct(cd_e) - (electrolyte_div_flux(c_e) + src);

    // BV 行：i_surf − 2·i0·sinh(F·η_eff/2RT) = 0，
    // η_eff = φ_s − φ_e − U − (j + j_sei)·L·R_sei（SEI 膜电阻压降）
    Vector eta_ln(m_neg_); // η 分布（lumped 产热计算复用）
    Vector eta_lp(m_pos_);
    Vector cs_ln(m_neg_);  // 表面浓度（lumped 产热 dU/dT LUT 求值复用）
    Vector cs_lp(m_pos_);
    for (int k = 0; k < m_neg_; ++k) {
        const real_t c_surf = std::clamp(1.5 * c_n(n_ - 1, k)
                                             - 0.5 * c_n(n_ - 2, k),
                                         0.0, solid.neg_max_conc);
        cs_ln(k) = c_surf;
        const real_t i0 = exchange_current_density_chen2020_neg(
            c_e(k), c_surf, solid.neg_max_conc, T);
        const real_t u = solid.neg_ocp_at(c_surf / solid.neg_max_conc);
        const real_t eta = phi_sn(k) - phi_e(k) - u;
        eta_ln(k) = eta;
        if (deg_on) {
            const real_t eta_eff = eta - (is_n(k) + j_sei_ln(k))
                                         * l_sei_y(k) * deg.r_sei;
            rb_n(k) = is_n(k) - 2.0 * i0 * safe_sinh(bv_coef * eta_eff);
        } else {
            rb_n(k) = is_n(k) - 2.0 * i0 * safe_sinh(bv_coef * eta);
        }
    }
    for (int k = 0; k < m_pos_; ++k) {
        const real_t c_surf = std::clamp(1.5 * c_p(n_ - 1, k)
                                             - 0.5 * c_p(n_ - 2, k),
                                         0.0, solid.pos_max_conc);
        cs_lp(k) = c_surf;
        const real_t i0 = exchange_current_density_chen2020_pos(
            c_e(pos_start_ + k), c_surf, solid.pos_max_conc, T);
        const real_t u = solid.pos_ocp_at(c_surf / solid.pos_max_conc);
        const real_t eta = phi_sp(k) - phi_e(pos_start_ + k) - u;
        eta_lp(k) = eta;
        rb_p(k) = is_p(k) - 2.0 * i0 * safe_sinh(bv_coef * eta);
    }

    // 电流守恒行：Σ a·dx·(i_surf + j_sei) = ±i_app（含 SEI 副反应分流）
    res(off_d) = a_n_ * dx_ * (is_n.sum() + j_sei_ln.sum()) - i_app;
    res(off_d + 1) = a_p_ * dx_ * is_p.sum() + i_app;

    // 外电路行（增广端口末行）：CC 锁端口电流 / CV 锁端电压；
    // V 与 finalize_fields 同口径（正极集流体端固相电位 − 负极基准 δ_n）
    if (port_control_ == solver::PortControl::Voltage) {
        res(dae_size() - 1) = (phi_sp(m_pos_ - 1) - delta_n) - v_target_;
    } else {
        res(dae_size() - 1) = i_cell - dae_current_;
    }

    // SEI 厚度微分行：L' = −V_bar/(F·z)·j_sei（j_sei < 0 → 生长）
    if (deg_on) {
        const real_t dldt_coef = -deg.v_bar_sei / (kFaraday * deg.z_sei);
        for (int k = 0; k < m_neg_; ++k) {
            res(off_sei + k) = yp(off_sei + k) - dldt_coef * j_sei_ln(k);
        }
    }

    // 热行（lumped）：T' = (Q − hA(T − T_amb)) / C_th（全隐式，IDA 自适应）
    if (lumped) {
        Vector is_n_v = is_n; // Map → Vector（compute_heat_source 入参）
        Vector is_p_v = is_p;
        HeatSourceArgs hs;
        hs.m = m_;
        hs.m_neg = m_neg_;
        hs.pos_start = pos_start_;
        hs.dx = dx_;
        hs.area = area_;
        hs.a_n = a_n_;
        hs.a_p = a_p_;
        hs.sigma_n_eff = sigma_n_eff_;
        hs.sigma_p_eff = sigma_p_eff_;
        hs.i_app = i_app;
        hs.T = T;
        hs.entropic_n = params_.thermal_params.entropic_change_neg;
        hs.entropic_p = params_.thermal_params.entropic_change_pos;
        hs.i_e_face = &i_e;   // 局部积分场（同 integrate_fields 口径）
        hs.phi_e = &phi_e;
        hs.i_surf_n = &is_n_v;
        hs.i_surf_p = &is_p_v;
        hs.eta_n = &eta_ln;
        hs.eta_p = &eta_lp;
        hs.lut_n = &params_.thermal_params.entropic_lut_neg;
        hs.lut_p = &params_.thermal_params.entropic_lut_pos;
        hs.c_surf_n = &cs_ln;
        hs.c_surf_p = &cs_lp;
        hs.cs_max_n = solid.neg_max_conc;
        hs.cs_max_p = solid.pos_max_conc;
        if (deg_on) {
            hs.j_sei_n = &j_sei_ln;
            hs.eta_sei_n = &eta_sei_ln;
        }
        real_t q_ohm = 0.0;
        real_t q_rxn = 0.0;
        real_t q_rev = 0.0;
        compute_heat_source(hs, q_ohm, q_rxn, q_rev);
        const auto& th = params_.thermal_params;
        const real_t c_th = rho_c_p_eff_ * th.cell_volume;
        const real_t cool = th.h_total * th.cooling_area;
        res(off_t) = yp(off_t)
                   - (q_ohm + q_rxn + q_rev - cool * (T - th.t_ambient)) / c_th;
    }
}

void DFNModel::dae_load_state(real_t t, const Eigen::Ref<const Vector>& y) {
    const int s_n = n_ * m_neg_;
    const int s_p = n_ * m_pos_;
    const int off_e = s_n + s_p;
    const int off_isn = off_e + m_;
    const int off_isp = off_isn + m_neg_;
    const int off_d = off_isp + m_pos_;
    const bool lumped = params_.thermal_type == ThermalType::Lumped;
    const bool deg_on = params_.degradation.enabled;
    const int off_sei = off_d + 2 + (lumped ? 1 : 0);

    Eigen::Map<const Eigen::MatrixXd> c_n(y.data(), n_, m_neg_);
    Eigen::Map<const Eigen::MatrixXd> c_p(y.data() + s_n, n_, m_pos_);
    Eigen::Map<const Vector> c_e(y.data() + off_e, m_);
    Eigen::Map<const Vector> is_n(y.data() + off_isn, m_neg_);
    Eigen::Map<const Vector> is_p(y.data() + off_isp, m_pos_);

    c_s_neg_ = c_n.cwiseMax(0.0).cwiseMin(params_.solid.neg_max_conc);
    c_s_pos_ = c_p.cwiseMax(0.0).cwiseMin(params_.solid.pos_max_conc);
    c_e_ = c_e.cwiseMax(kCeFloor);
    i_surf_n_ = is_n;
    i_surf_p_ = is_p;
    delta_n_ = y(off_d);
    delta_p_ = y(off_d + 1);
    i_port_ = y(dae_size() - 1); // 增广端口电流（CV 模式下为 Newton 解）
    // 导出量同步的电流口径：CC 模式下用目标电流（与 Gummel 重估的
    // i_surf 一致，维持既有行为）；CV 模式下端口电流是解的一部分，
    // 必须用状态中的 I_cell
    const real_t i_eff = (port_control_ == solver::PortControl::Voltage)
                       ? i_port_ : dae_current_;
    if (lumped) {
        t_cell_ = y(off_d + 2); // 先于 finalize_fields（φ_e 浓差项依赖 T）
    }
    if (deg_on) {
        l_sei_ = y.segment(off_sei, m_neg_).cwiseMax(0.0);
    }
    time_ = t;

    refresh_surface_quantities();
    finalize_fields(i_eff / area_);
    if (deg_on) {
        refresh_sei_fields(); // j_sei/eta_sei 与载入状态同步
    }
    refresh_heat(i_eff / area_); // 产热缓存与载入状态同步
}

// ══════════════════════════════════════════════════════════════
//  SEI 降解（ec reaction limited；PyBaMM sei_growth.py Yang2017）
// ══════════════════════════════════════════════════════════════

real_t DFNModel::sei_kinetics_neg(int k, real_t L, real_t eta_bare,
                                  real_t i_surf_k, real_t T,
                                  real_t& g, real_t& q) const {
    const auto& deg = params_.degradation;
    const real_t eta_sei = eta_bare - deg.u_sei - i_surf_k * L * deg.r_sei;
    const real_t K = deg.k_sei
        * std::exp(-deg.alpha_sei * kFaraday / (kGasConst * T) * eta_sei);
    const real_t arrh = std::exp(deg.e_sei / kGasConst
                                 * (1.0 / deg.t_ref - 1.0 / T));
    const real_t denom = 1.0 + L * K / deg.d_ec;
    // PyBaMM 口径：arrh 乘整体（分母不含）
    const real_t j_sei = -kFaraday * deg.c_ec_0 * K * arrh / denom;
    // dL/dt = −V_bar/(F·z)·j_sei = g/(1 + q·L)
    g = deg.v_bar_sei * deg.c_ec_0 * K * arrh / deg.z_sei;
    q = K / deg.d_ec;
    return j_sei;
}

void DFNModel::refresh_sei_fields() {
    if (!params_.degradation.enabled || m_neg_ == 0) {
        return;
    }
    for (int k = 0; k < m_neg_; ++k) {
        const real_t eta_bare = phi_s_n_(k) - phi_e_(k);
        real_t g = 0.0;
        real_t q = 0.0;
        j_sei_n_(k) = sei_kinetics_neg(k, l_sei_(k), eta_bare,
                                       i_surf_n_(k), t_cell_, g, q);
        eta_sei_n_(k) = eta_bare - params_.degradation.u_sei
                      - i_surf_n_(k) * l_sei_(k) * params_.degradation.r_sei;
    }
}

void DFNModel::advance_sei_cn(real_t dt) {
    const auto& deg = params_.degradation;
    for (int k = 0; k < m_neg_; ++k) {
        const real_t eta_bare = phi_s_n_(k) - phi_e_(k);
        real_t g = 0.0;
        real_t q = 0.0;
        j_sei_n_(k) = sei_kinetics_neg(k, l_sei_(k), eta_bare,
                                       i_surf_n_(k), t_cell_, g, q);
        eta_sei_n_(k) = eta_bare - deg.u_sei - i_surf_n_(k) * l_sei_(k) * deg.r_sei;
        const real_t L = l_sei_(k);
        // dL/dt = g/(1 + q·L)：完全隐式一元二次根（步内 φ/T/i_surf 冻结）
        if (q * L > 1e-12) {
            const real_t b = 1.0 - q * L;
            const real_t c = -(L + dt * g);
            l_sei_(k) = (-b + std::sqrt(b * b - 4.0 * q * c)) / (2.0 * q);
        } else {
            l_sei_(k) = L + dt * g;
        }
        // 新厚度下刷新副反应电流/过电位（产热与快照诊断一致）
        j_sei_n_(k) = sei_kinetics_neg(k, l_sei_(k), eta_bare,
                                       i_surf_n_(k), t_cell_, g, q);
    }
}

real_t DFNModel::sei_thickness_avg() const {
    return l_sei_.size() > 0 ? l_sei_.mean() : 0.0;
}

real_t DFNModel::lithium_loss_ah() const {
    // PyBaMM Q_sei = Δc_sei·V_k·z，c_sei = L·a/V_bar（逐格体平均）→
    // Q = z/V_bar·A·dx·Σ a·(L − L0)（只计初始厚度之后的增长）
    if (!params_.degradation.enabled || l_sei_.size() == 0) {
        return 0.0;
    }
    const auto& deg = params_.degradation;
    real_t dL_sum = 0.0;
    for (int k = 0; k < m_neg_; ++k) {
        dL_sum += l_sei_(k) - deg.l_sei_0;
    }
    return deg.z_sei / deg.v_bar_sei * a_n_ * dx_ * area_ * dL_sum
         * kFaraday / 3600.0;
}

real_t DFNModel::degraded_capacity_ah() const {
    return params_.solid.capacity - lithium_loss_ah();
}

} // namespace fluxcell::model
