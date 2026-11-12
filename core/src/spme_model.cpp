/// @file spme_model.cpp
/// @brief 扩展单颗粒模型 (SPMe) 实现：固相颗粒 + 电解液相 1D 传输
#include <fluxcell/core/model/spme_model.hpp>

#include <fluxcell/core/discretization/discretization.hpp>
#include <fluxcell/core/physics/physics.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace fluxcell::model {

namespace {
using fluxcell::physics::exchange_current_density_chen2020_neg;
using fluxcell::physics::exchange_current_density_chen2020_pos;
using fluxcell::physics::graphite_chen2020_ocp;
using fluxcell::physics::nmc_chen2020_ocp;
using fluxcell::physics::overpotential_from_bv;

constexpr real_t kFaraday = 96485.0;  // C/mol
constexpr real_t kPi = 3.14159265358979323846;
constexpr real_t kGasConst = 8.314;   // J/(mol·K)
constexpr real_t kCeFloor = 1.0;      // 电解液浓度数值下限 [mol/m³]
constexpr real_t kKappaFloor = 1e-6;  // κ 数值下限 [S/m]（盐耗尽 c_e→0 防除零）
} // namespace

// ══════════════════════════════════════════════════════════════
//  Crank-Nicolson 预分解算子（颗粒，同 SPM）
// ══════════════════════════════════════════════════════════════

struct SPMeModel::CNOperators {
    Eigen::SparseLU<SparseMatrix> solve_neg;
    Eigen::SparseLU<SparseMatrix> solve_pos;
    SparseMatrix plus_neg;
    SparseMatrix plus_pos;
    real_t surf_vol_neg = 0.0;
    real_t surf_area_neg = 0.0;
    real_t surf_vol_pos = 0.0;
    real_t surf_area_pos = 0.0;
};

// ══════════════════════════════════════════════════════════════
//  参数：浓度相关电解液输运
// ══════════════════════════════════════════════════════════════

real_t SPMeModel::Parameters::conductivity_at(real_t c_e) const {
    if (!conductivity_lut.empty()) {
        return fluxcell::physics::lut_interpolate(electrolyte_conc_lut,
                                                  conductivity_lut, c_e);
    }
    return electrolyte_conductivity;
}

real_t SPMeModel::Parameters::diffusivity_at(real_t c_e) const {
    if (!diffusivity_lut.empty()) {
        return fluxcell::physics::lut_interpolate(electrolyte_conc_lut,
                                                  diffusivity_lut, c_e);
    }
    return electrolyte_diffusivity;
}

void SPMeModel::Parameters::validate_transport_table() const {
    const auto check = [&](const std::vector<real_t>& values, const char* name) {
        if (values.empty()) {
            return;
        }
        if (electrolyte_conc_lut.size() < 2
            || values.size() != electrolyte_conc_lut.size()) {
            throw std::invalid_argument(
                std::string("SPMeModel: electrolyte transport table (")
                + name + ") 需与浓度节点等长且 ≥ 2 点");
        }
        for (std::size_t k = 1; k < electrolyte_conc_lut.size(); ++k) {
            if (!(electrolyte_conc_lut[k] > electrolyte_conc_lut[k - 1])) {
                throw std::invalid_argument(
                    "SPMeModel: electrolyte transport table 浓度节点须严格升序");
            }
        }
    };
    check(conductivity_lut, "conductivity");
    check(diffusivity_lut, "diffusivity");
}

SPMeModel::SPMeModel(const Parameters& params, real_t initial_current)
    : params_(params) {
    const auto& solid = params.solid;

    n_ = solid.particle_cells;
    if (n_ < 3) {
        throw std::invalid_argument("SPMeModel: particle_cells must be >= 3");
    }
    m_ = params.electrolyte_cells;
    if (m_ < 10) {
        throw std::invalid_argument("SPMeModel: electrolyte_cells must be >= 10");
    }
    if (solid.temperature <= 0.0) {
        throw std::invalid_argument("SPMeModel: temperature must be positive");
    }
    if (solid.initial_soc < 0.0 || solid.initial_soc > 1.0) {
        throw std::invalid_argument("SPMeModel: initial_soc must be in [0, 1]");
    }
    if (params.porosity_neg <= 0.0 || params.porosity_sep <= 0.0
        || params.porosity_pos <= 0.0) {
        throw std::invalid_argument("SPMeModel: porosity must be positive");
    }
    if (params.electrolyte_conductivity <= 0.0
        || params.electrolyte_diffusivity <= 0.0) {
        throw std::invalid_argument("SPMeModel: electrolyte transport must be positive");
    }
    params.validate_transport_table();

    area_ = solid.electrode_height * solid.electrode_width;

    // 电极比表面积 a = 3·ε_s/R [m²/m³]
    a_neg_ = 3.0 * solid.neg_active_fraction / solid.neg_particle_radius;
    a_pos_ = 3.0 * solid.pos_active_fraction / solid.pos_particle_radius;

    // 颗粒径向网格 + 扩散算子（同 SPM）
    auto mesh_neg = discretization::SphericalMesh::create(solid.neg_particle_radius, n_);
    auto mesh_pos = discretization::SphericalMesh::create(solid.pos_particle_radius, n_);

    D_neg_ = discretization::assemble_spherical_diffusion(mesh_neg, solid.neg_diffusion);
    D_pos_ = discretization::assemble_spherical_diffusion(mesh_pos, solid.pos_diffusion);

    cn_ = std::make_unique<CNOperators>();
    cn_->surf_vol_neg = mesh_neg.volumes()[n_ - 1];
    cn_->surf_area_neg = mesh_neg.face_areas()[n_];
    cn_->surf_vol_pos = mesh_pos.volumes()[n_ - 1];
    cn_->surf_area_pos = mesh_pos.face_areas()[n_];

    vol_neg_.assign(mesh_neg.volumes().begin(), mesh_neg.volumes().end());
    vol_pos_.assign(mesh_pos.volumes().begin(), mesh_pos.volumes().end());

    // 初始固相浓度：初始 SOC 线性映射化学计量端点（同 SPM）
    const real_t x0 = solid.neg_stoich_at_soc0
                    + solid.initial_soc
                          * (solid.neg_stoich_at_soc1 - solid.neg_stoich_at_soc0);
    const real_t y0 = solid.pos_stoich_at_soc0
                    + solid.initial_soc
                          * (solid.pos_stoich_at_soc1 - solid.pos_stoich_at_soc0);
    c_s_neg_.setConstant(n_, x0 * solid.neg_max_conc);
    c_s_pos_.setConstant(n_, y0 * solid.pos_max_conc);

    // ── 电解液 1D 网格：均匀 dx，区域边界对齐最近格边界 ──
    const real_t L = solid.neg_thickness + solid.sep_thickness + solid.pos_thickness;
    dx_ = L / m_;

    const auto cells_for = [&](real_t thickness) {
        return std::clamp(static_cast<int>(std::lround(thickness / dx_)), 1,
                          m_ - 2);
    };
    sep_start_ = cells_for(solid.neg_thickness);
    pos_start_ = std::clamp(sep_start_ + cells_for(solid.sep_thickness),
                            sep_start_ + 1, m_ - 1);

    eps_e_ = Vector::Constant(m_, params.porosity_sep);
    eps_e_.segment(0, sep_start_).setConstant(params.porosity_neg);
    eps_e_.segment(pos_start_, m_ - pos_start_).setConstant(params.porosity_pos);

    c_e_.setConstant(m_, params.c_e_init);
    phi_e_.setZero(m_);

    dae_current_ = initial_current;

    cn_dt_ = -1.0;

    update_voltage(initial_current);
}

Vector SPMeModel::electrolyte_div_flux(const Vector& c) const {
    // 守恒通量形式三对角散度：D_e^eff = D(c_e)·ε^bruggeman（逐格），
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

real_t SPMeModel::kappa_eff_at(int i) const {
    // κ^eff = κ(c_e)·ε^bruggeman；κ 下限防盐耗尽（c_e → 0）欧姆项除零
    return std::max(params_.conductivity_at(c_e_(i)), kKappaFloor)
         * std::pow(eps_e_(i), params_.bruggeman_exponent);
}

SPMeModel::~SPMeModel() = default;

real_t SPMeModel::neg_surface_conc() const {
    // FVM 边界浓度：二阶外推（同 SPM，与 PyBaMM edge 约定对齐）
    const real_t c_last = c_s_neg_(n_ - 1);
    const real_t c_prev = c_s_neg_(n_ - 2);
    return std::clamp(1.5 * c_last - 0.5 * c_prev, 0.0, params_.solid.neg_max_conc);
}

real_t SPMeModel::pos_surface_conc() const {
    const real_t c_last = c_s_pos_(n_ - 1);
    const real_t c_prev = c_s_pos_(n_ - 2);
    return std::clamp(1.5 * c_last - 0.5 * c_prev, 0.0, params_.solid.pos_max_conc);
}

void SPMeModel::ensure_cn(real_t dt) {
    const real_t half_dt = 0.5 * dt;

    if (cn_dt_ != dt) {
        // ── 颗粒算子（同 SPM，dt 不变时复用）──
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
                throw std::runtime_error("SPMeModel: particle CN factorization failed");
            }
        };
        build_particle(D_neg_, cn_->plus_neg, cn_->solve_neg);
        build_particle(D_pos_, cn_->plus_pos, cn_->solve_pos);

        cn_dt_ = dt;
    }

    // ── 电解液算子：ε dc/dt = L(c_e)·c + S，L 为守恒通量形式三对角 ──
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

        // L 行：[d_left, -(d_left+d_right), d_right] / dx²
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
        throw std::runtime_error("SPMeModel: electrolyte CN factorization failed");
    }
}

void SPMeModel::step(real_t current, real_t dt) {
    if (!(dt > 0.0)) {
        throw std::invalid_argument("SPMeModel::step: dt must be positive");
    }
    ensure_cn(dt);
    const auto& solid = params_.solid;

    // ── 固相颗粒（同 SPM：CN + 表面源）──
    const real_t N_neg = current / (area_ * a_neg_ * solid.neg_thickness) / kFaraday;
    const real_t N_pos = -current / (area_ * a_pos_ * solid.pos_thickness) / kFaraday;

    Vector b_neg = cn_->plus_neg * c_s_neg_;
    b_neg(n_ - 1) -= dt * cn_->surf_area_neg / cn_->surf_vol_neg * N_neg;
    Vector b_pos = cn_->plus_pos * c_s_pos_;
    b_pos(n_ - 1) -= dt * cn_->surf_area_pos / cn_->surf_vol_pos * N_pos;

    c_s_neg_ = cn_->solve_neg.solve(b_neg);
    c_s_pos_ = cn_->solve_pos.solve(b_pos);

    c_s_neg_ = c_s_neg_.cwiseMax(0.0).cwiseMin(solid.neg_max_conc);
    c_s_pos_ = c_s_pos_.cwiseMax(0.0).cwiseMin(solid.pos_max_conc);

    // ── 液相电解液（CN + 区域源项）──
    // 源 = (1−t₊)/F · ∂i_e/∂x /A = ±(1−t₊)·i/(F·A·L_region)，区域长度取实际网格长度
    Vector s_vec = Vector::Zero(m_);
    if (current != 0.0) {
        const real_t L_neg = sep_start_ * dx_;
        const real_t L_pos = (m_ - pos_start_) * dx_;
        const real_t s_n = (1.0 - params_.transference_number) * current
                           / (kFaraday * area_ * L_neg);
        const real_t s_p = -s_n * (L_neg / L_pos);
        s_vec.segment(0, sep_start_).setConstant(s_n);
        s_vec.segment(pos_start_, m_ - pos_start_).setConstant(s_p);
    }

    Vector b_e = le_plus_ * c_e_ + dt * s_vec;
    c_e_ = le_solve_.solve(b_e);
    if (le_solve_.info() != Eigen::Success) {
        throw std::runtime_error("SPMeModel: electrolyte CN solve failed");
    }
    c_e_ = c_e_.cwiseMax(kCeFloor);

    time_ += dt;
    update_voltage(current);
}

void SPMeModel::update_voltage(real_t current) {
    const auto& solid = params_.solid;
    const real_t T = solid.temperature;

    // 颗粒表面电流密度 [A/m²]（正 = 锂流出颗粒，放电方向）
    const real_t i_neg = current / (area_ * a_neg_ * solid.neg_thickness);
    const real_t i_pos = -current / (area_ * a_pos_ * solid.pos_thickness);

    // 交换电流密度使用电极区域 x 平均电解液浓度（PyBaMM SPMe 约定）
    const real_t c_e_neg = c_e_.segment(0, sep_start_).mean();
    const real_t c_e_pos = c_e_.segment(pos_start_, m_ - pos_start_).mean();

    const real_t i0_neg = exchange_current_density_chen2020_neg(
        c_e_neg, neg_surface_conc(), solid.neg_max_conc, T);
    const real_t i0_pos = exchange_current_density_chen2020_pos(
        c_e_pos, pos_surface_conc(), solid.pos_max_conc, T);

    eta_neg_ = overpotential_from_bv(i_neg, i0_neg, T);
    eta_pos_ = overpotential_from_bv(i_pos, i0_pos, T);

    u_neg_ = params_.solid.neg_ocp_at(neg_stoichiometry());
    u_pos_ = params_.solid.pos_ocp_at(pos_stoichiometry());

    // ── 液相过电位（逐格解析积分，κ^eff = κ(c_e)·ε^bruggeman 浓度相关）──
    // η_e = φ_e(L) − φ_e(0) = −∫ i_e/κ^eff dx + (2RT/F)(1−t₊)·ln(c_e(L)/c_e(0))
    // i_e 为液相面电流密度 [A/m²]，分段线性：负极 0→i/A、隔膜 i/A、正极 i/A→0
    const real_t L_neg = sep_start_ * dx_;
    const real_t L_sep = (pos_start_ - sep_start_) * dx_;
    const real_t L_pos = (m_ - pos_start_) * dx_;
    const real_t L = L_neg + L_sep + L_pos;
    const real_t j_e = current / area_; // 液相面电流密度 [A/m²]（放电 > 0）

    // 液相电位分布（φ_e(0) = 0 基准，逐格累加）；末值即 η_e，与标量严格一致
    phi_e_(0) = 0.0;
    const real_t dphi_conc_coef = 2.0 * kGasConst * T / kFaraday
                                  * (1.0 - params_.transference_number);
    for (int i = 1; i < m_; ++i) {
        const real_t x_face = i * dx_;
        real_t i_e;
        if (x_face <= L_neg) {
            i_e = j_e * (x_face / L_neg);
        } else if (x_face <= L_neg + L_sep) {
            i_e = j_e;
        } else {
            i_e = j_e * (L - x_face) / L_pos;
        }
        // 界面 κ^eff：两侧格算术平均（浓度相关）
        const real_t kappa_face = 0.5 * (kappa_eff_at(i - 1) + kappa_eff_at(i));
        const real_t dphi_ohmic = -i_e / kappa_face * dx_;
        const real_t dphi_conc = dphi_conc_coef * (std::log(c_e_(i)) - std::log(c_e_(i - 1)));
        phi_e_(i) = phi_e_(i - 1) + dphi_ohmic + dphi_conc;
    }
    eta_e_ = phi_e_(m_ - 1);

    voltage_ = u_pos_ - u_neg_ + eta_pos_ - eta_neg_ + eta_e_;
}

real_t SPMeModel::total_lithium_mol() const {
    const auto& solid = params_.solid;

    // 固相（同 SPM）
    auto electrode_mol = [&](const Vector& c, const std::vector<real_t>& volumes,
                             real_t len, real_t eps_s, real_t radius) {
        const real_t particle_volume = (4.0 / 3.0) * kPi * radius * radius * radius;
        const real_t n_particles = area_ * len * eps_s / particle_volume;
        real_t mol_per_particle = 0.0;
        for (int i = 0; i < n_; ++i) {
            mol_per_particle += volumes[i] * c(i);
        }
        return n_particles * mol_per_particle;
    };

    // 液相：Σ ε·c·A·dx（全边界零通量 + 对称源 → 恒守恒）
    real_t liquid_mol = 0.0;
    for (int i = 0; i < m_; ++i) {
        liquid_mol += eps_e_(i) * c_e_(i) * dx_;
    }

    return electrode_mol(c_s_neg_, vol_neg_, solid.neg_thickness,
                         solid.neg_active_fraction, solid.neg_particle_radius)
         + electrode_mol(c_s_pos_, vol_pos_, solid.pos_thickness,
                         solid.pos_active_fraction, solid.pos_particle_radius)
         + area_ * liquid_mol;
}

SolutionSnapshot SPMeModel::snapshot(real_t current) const {
    SolutionSnapshot snap;
    snap.time = time_;
    snap.voltage = voltage_;
    snap.current = current;
    snap.temperature = params_.solid.temperature;
    snap.c_e = c_e_;
    snap.phi_e = phi_e_;
    snap.sep_start = sep_start_;
    snap.pos_start = pos_start_;
    snap.eta_e = eta_e_;
    snap.c_s_surf_neg = Vector::Constant(1, neg_surface_conc());
    snap.c_s_surf_pos = Vector::Constant(1, pos_surface_conc());
    snap.c_s_neg = c_s_neg_;
    snap.c_s_pos = c_s_pos_;
    snap.x_surf = neg_stoichiometry();
    snap.y_surf = pos_stoichiometry();
    snap.eta_neg = eta_neg_;
    snap.eta_pos = eta_pos_;
    snap.u_neg = u_neg_;
    snap.u_pos = u_pos_;
    return snap;
}

// ══════════════════════════════════════════════════════════════
//  DAE 视图（SUNDIALS IDA）：纯 ODE，状态 = [c_s_n, c_s_p, c_e]
//  ODE 右端：c_s' = D·c_s − e_surf·N；ε c_e' = L·c_e + s（区域恒源）
// ══════════════════════════════════════════════════════════════

int SPMeModel::dae_size() const {
    return 2 * n_ + m_;
}

void SPMeModel::dae_mask(Vector& id) const {
    id.setOnes(dae_size()); // 全部为微分变量
}

void SPMeModel::dae_atol(Vector& atol) const {
    atol.resize(dae_size());
    atol.segment(0, n_).setConstant(0.5);      // c_s_n [mol/m³]
    atol.segment(n_, n_).setConstant(1.0);     // c_s_p [mol/m³]
    atol.segment(2 * n_, m_).setConstant(0.05); // c_e [mol/m³]
}

Vector SPMeModel::dae_electrolyte_source() const {
    // 源 = (1−t₊)/F · ∂i_e/∂x /A（解析电流分布，与 CN 路径同一公式）
    Vector s = Vector::Zero(m_);
    if (dae_current_ != 0.0) {
        const real_t L_neg = sep_start_ * dx_;
        const real_t L_pos = (m_ - pos_start_) * dx_;
        const real_t s_n = (1.0 - params_.transference_number) * dae_current_
                         / (kFaraday * area_ * L_neg);
        s.segment(0, sep_start_).setConstant(s_n);
        s.segment(pos_start_, m_ - pos_start_).setConstant(-s_n * (L_neg / L_pos));
    }
    return s;
}

void SPMeModel::dae_pack_state(Vector& y, Vector& yp) const {
    const int N = dae_size();
    y.resize(N);
    yp.resize(N);
    y.segment(0, n_) = c_s_neg_;
    y.segment(n_, n_) = c_s_pos_;
    y.segment(2 * n_, m_) = c_e_;

    const auto& solid = params_.solid;
    const real_t N_neg = dae_current_
                       / (area_ * a_neg_ * solid.neg_thickness) / kFaraday;
    const real_t N_pos = -dae_current_
                        / (area_ * a_pos_ * solid.pos_thickness) / kFaraday;
    yp.segment(0, n_) = D_neg_ * c_s_neg_;
    yp(n_ - 1) -= cn_->surf_area_neg / cn_->surf_vol_neg * N_neg;
    yp.segment(n_, n_) = D_pos_ * c_s_pos_;
    yp(2 * n_ - 1) -= cn_->surf_area_pos / cn_->surf_vol_pos * N_pos;
    yp.segment(2 * n_, m_) =
        (electrolyte_div_flux(c_e_) + dae_electrolyte_source()).cwiseQuotient(eps_e_);
}

void SPMeModel::dae_residual(real_t t, const Eigen::Ref<const Vector>& y,
                                 const Eigen::Ref<const Vector>& yp,
                                 Eigen::Ref<Vector> res) {
    (void)t;
    const auto& solid = params_.solid;
    const real_t N_neg = dae_current_
                       / (area_ * a_neg_ * solid.neg_thickness) / kFaraday;
    const real_t N_pos = -dae_current_
                        / (area_ * a_pos_ * solid.pos_thickness) / kFaraday;

    Eigen::Map<const Vector> c_n(y.data(), n_);
    Eigen::Map<const Vector> c_p(y.data() + n_, n_);
    Eigen::Map<const Vector> c_e(y.data() + 2 * n_, m_);
    Eigen::Map<const Vector> d_n(yp.data(), n_);
    Eigen::Map<const Vector> d_p(yp.data() + n_, n_);
    Eigen::Map<const Vector> d_e(yp.data() + 2 * n_, m_);
    Eigen::Map<Vector> r_n(res.data(), n_);
    Eigen::Map<Vector> r_p(res.data() + n_, n_);
    Eigen::Map<Vector> r_e(res.data() + 2 * n_, m_);

    r_n = d_n - (D_neg_ * c_n);
    r_n(n_ - 1) += cn_->surf_area_neg / cn_->surf_vol_neg * N_neg;
    r_p = d_p - (D_pos_ * c_p);
    r_p(n_ - 1) += cn_->surf_area_pos / cn_->surf_vol_pos * N_pos;
    // 液相微分行：ε·y' − (L(c)·c + s)，D_e 浓度相关
    r_e = eps_e_.cwiseProduct(d_e)
        - (electrolyte_div_flux(c_e) + dae_electrolyte_source());
}

void SPMeModel::dae_load_state(real_t t, const Eigen::Ref<const Vector>& y) {
    const auto& solid = params_.solid;
    Eigen::Map<const Vector> c_n(y.data(), n_);
    Eigen::Map<const Vector> c_p(y.data() + n_, n_);
    Eigen::Map<const Vector> c_e(y.data() + 2 * n_, m_);
    c_s_neg_ = c_n.cwiseMax(0.0).cwiseMin(solid.neg_max_conc);
    c_s_pos_ = c_p.cwiseMax(0.0).cwiseMin(solid.pos_max_conc);
    c_e_ = c_e.cwiseMax(kCeFloor);
    time_ = t;
    update_voltage(dae_current_);
}

} // namespace fluxcell::model
