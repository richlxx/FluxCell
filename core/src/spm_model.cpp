/// @file spm_model.cpp
/// @brief 单颗粒模型 (SPM) 实现：球坐标 FVM + Crank-Nicolson 扩散 + BV arcsinh 反解
///
/// 参考 PyBaMM pybamm/models/full_battery_models/lithium_ion/basic.py (SPM)：
///   ∂c_s/∂t = (1/r²) ∂/∂r (r² D_s ∂c_s/∂r)，r ∈ [0, R]
///   边界：∂c_s/∂r|_{r=0} = 0，-D_s ∂c_s/∂r|_{r=R} = j/F
///   η = (2RT/F) asinh(i / (2 i0))，V = U_p - U_n + η_p - η_n
#include <fluxcell/core/model/spm_model.hpp>

#include <fluxcell/core/discretization/discretization.hpp>
#include <fluxcell/core/physics/physics.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fluxcell::model {

namespace {
using fluxcell::physics::exchange_current_density_chen2020_neg;
using fluxcell::physics::exchange_current_density_chen2020_pos;
using fluxcell::physics::graphite_chen2020_ocp;
using fluxcell::physics::nmc_chen2020_ocp;
using fluxcell::physics::overpotential_from_bv;

constexpr real_t kFaraday = 96485.0; // C/mol
constexpr real_t kPi = 3.14159265358979323846;
} // namespace

// ══════════════════════════════════════════════════════════════
//  OCP 查表/内置公式双路径
// ══════════════════════════════════════════════════════════════

real_t SPMModel::Parameters::neg_ocp_at(real_t x) const {
    return neg_ocp_stoich.empty()
               ? physics::graphite_chen2020_ocp(x)
               : physics::ocp_from_lut(neg_ocp_stoich, neg_ocp_voltage, x);
}

real_t SPMModel::Parameters::pos_ocp_at(real_t y) const {
    return pos_ocp_stoich.empty()
               ? physics::nmc_chen2020_ocp(y)
               : physics::ocp_from_lut(pos_ocp_stoich, pos_ocp_voltage, y);
}

// ══════════════════════════════════════════════════════════════
//  Crank-Nicolson 预分解算子
// ══════════════════════════════════════════════════════════════

struct SPMModel::CNOperators {
    // 每电极：(I - dt/2·D) 的 LU 分解、(I + dt/2·D) 矩阵、表面控制体几何
    Eigen::SparseLU<SparseMatrix> solve_neg;
    Eigen::SparseLU<SparseMatrix> solve_pos;
    SparseMatrix plus_neg;
    SparseMatrix plus_pos;
    real_t surf_vol_neg = 0.0;   ///< 负极颗粒表面控制体体积 [m³]
    real_t surf_area_neg = 0.0;  ///< 负极颗粒外表面面积 [m²]
    real_t surf_vol_pos = 0.0;   ///< 正极颗粒表面控制体体积 [m³]
    real_t surf_area_pos = 0.0;  ///< 正极颗粒外表面面积 [m²]
};

SPMModel::SPMModel(const Parameters& params, real_t initial_current)
    : params_(params) {
    n_ = params.particle_cells;
    if (n_ < 3) {
        throw std::invalid_argument("SPMModel: particle_cells must be >= 3");
    }
    if (params.temperature <= 0.0) {
        throw std::invalid_argument("SPMModel: temperature must be positive");
    }
    if (params.initial_soc < 0.0 || params.initial_soc > 1.0) {
        throw std::invalid_argument("SPMModel: initial_soc must be in [0, 1]");
    }

    area_ = params.electrode_height * params.electrode_width;

    // 电极比表面积 a = 3·ε_s/R [m²/m³]
    a_neg_ = 3.0 * params.neg_active_fraction / params.neg_particle_radius;
    a_pos_ = 3.0 * params.pos_active_fraction / params.pos_particle_radius;

    // 球坐标网格 + 扩散算子（r=0 处零通量已隐含在 FVM 算子中）
    auto mesh_neg = discretization::SphericalMesh::create(params.neg_particle_radius, n_);
    auto mesh_pos = discretization::SphericalMesh::create(params.pos_particle_radius, n_);

    D_neg_ = discretization::assemble_spherical_diffusion(mesh_neg, params.neg_diffusion);
    D_pos_ = discretization::assemble_spherical_diffusion(mesh_pos, params.pos_diffusion);

    cn_ = std::make_unique<CNOperators>();
    cn_->surf_vol_neg = mesh_neg.volumes()[n_ - 1];
    cn_->surf_area_neg = mesh_neg.face_areas()[n_];
    cn_->surf_vol_pos = mesh_pos.volumes()[n_ - 1];
    cn_->surf_area_pos = mesh_pos.face_areas()[n_];

    vol_neg_.assign(mesh_neg.volumes().begin(), mesh_neg.volumes().end());
    vol_pos_.assign(mesh_pos.volumes().begin(), mesh_pos.volumes().end());

    cn_dt_ = -1.0; // 首次 step 时按实际 dt 构建/分解

    // 初始浓度（均匀）：由初始 SOC 线性映射化学计量端点
    // x_n(SOC) = x_0 + SOC·(x_1 − x_0)；y_p(SOC) = y_0 + SOC·(y_1 − y_0)
    const real_t x0 = params.neg_stoich_at_soc0
                    + params.initial_soc
                          * (params.neg_stoich_at_soc1 - params.neg_stoich_at_soc0);
    const real_t y0 = params.pos_stoich_at_soc0
                    + params.initial_soc
                          * (params.pos_stoich_at_soc1 - params.pos_stoich_at_soc0);
    c_s_neg_.setConstant(n_, x0 * params.neg_max_conc);
    c_s_pos_.setConstant(n_, y0 * params.pos_max_conc);

    dae_current_ = initial_current;
    update_voltage(initial_current);
}

SPMModel::~SPMModel() = default;

real_t SPMModel::neg_surface_conc() const {
    // FVM 边界浓度：均匀网格下 c_edge ≈ 1.5·c_last - 0.5·c_prev（二阶外推）。
    // PyBaMM 的颗粒变量定义在网格节点（edge）上，其表面浓度为边界值；
    // 直接用最外控制体平均会使表面浓度产生 O(dx) 系统偏差（放电时偏高）。
    const real_t c_last = c_s_neg_(n_ - 1);
    const real_t c_prev = c_s_neg_(n_ - 2);
    return std::clamp(1.5 * c_last - 0.5 * c_prev, 0.0, params_.neg_max_conc);
}

real_t SPMModel::pos_surface_conc() const {
    const real_t c_last = c_s_pos_(n_ - 1);
    const real_t c_prev = c_s_pos_(n_ - 2);
    return std::clamp(1.5 * c_last - 0.5 * c_prev, 0.0, params_.pos_max_conc);
}

void SPMModel::ensure_cn(real_t dt) {
    if (cn_dt_ == dt) {
        return;
    }

    const real_t half_dt = 0.5 * dt;

    // lhs = I - (dt/2)·D，plus = I + (dt/2)·D
    auto build = [this, half_dt](const SparseMatrix& D, SparseMatrix& plus,
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
            throw std::runtime_error("SPMModel: CN matrix factorization failed");
        }
    };

    build(D_neg_, cn_->plus_neg, cn_->solve_neg);
    build(D_pos_, cn_->plus_pos, cn_->solve_pos);

    cn_dt_ = dt;
}

void SPMModel::step(real_t current, real_t dt) {
    if (!(dt > 0.0)) {
        throw std::invalid_argument("SPMModel::step: dt must be positive");
    }
    ensure_cn(dt);

    // 表面摩尔通量 N = i/F [mol/(m²·s)]（正 = 锂流出颗粒，即放电方向）
    // SPM 中反应电流密度由外加电流唯一确定（锂守恒），BV 仅用于反解过电位
    const real_t N_neg = current / (area_ * a_neg_ * params_.neg_thickness) / kFaraday;
    const real_t N_pos = -current / (area_ * a_pos_ * params_.pos_thickness) / kFaraday;

    // Crank-Nicolson：c_new = (I - dt/2·D)⁻¹ [(I + dt/2·D)·c_old + dt·s]
    // 源项仅作用于表面控制体：s = -A_surf·N / V_surf
    Vector b_neg = cn_->plus_neg * c_s_neg_;
    b_neg(n_ - 1) -= dt * cn_->surf_area_neg / cn_->surf_vol_neg * N_neg;

    Vector b_pos = cn_->plus_pos * c_s_pos_;
    b_pos(n_ - 1) -= dt * cn_->surf_area_pos / cn_->surf_vol_pos * N_pos;

    c_s_neg_ = cn_->solve_neg.solve(b_neg);
    c_s_pos_ = cn_->solve_pos.solve(b_pos);

    // 限幅保证 OCP/BV 良定义（仅在接近满充/满放时触发，轻微破坏锂守恒）
    c_s_neg_ = c_s_neg_.cwiseMax(0.0).cwiseMin(params_.neg_max_conc);
    c_s_pos_ = c_s_pos_.cwiseMax(0.0).cwiseMin(params_.pos_max_conc);

    time_ += dt;
    update_voltage(current);
}

void SPMModel::update_voltage(real_t current) {
    const real_t T = params_.temperature;

    // 颗粒表面电流密度 [A/m²]（正 = 锂流出颗粒，即放电方向）
    const real_t i_neg = current / (area_ * a_neg_ * params_.neg_thickness);
    const real_t i_pos = -current / (area_ * a_pos_ * params_.pos_thickness);

    const real_t i0_neg = exchange_current_density_chen2020_neg(
        params_.electrolyte_conc, neg_surface_conc(), params_.neg_max_conc, T);
    const real_t i0_pos = exchange_current_density_chen2020_pos(
        params_.electrolyte_conc, pos_surface_conc(), params_.pos_max_conc, T);

    eta_neg_ = overpotential_from_bv(i_neg, i0_neg, T);
    eta_pos_ = overpotential_from_bv(i_pos, i0_pos, T);

    u_neg_ = params_.neg_ocp_at(neg_stoichiometry());
    u_pos_ = params_.pos_ocp_at(pos_stoichiometry());

    // SPM 端电压：无固相/液相欧姆损失（已用 PyBaMM SPM 数值验证）
    voltage_ = u_pos_ - u_neg_ + eta_pos_ - eta_neg_;
}

real_t SPMModel::total_lithium_mol() const {
    // 电极颗粒数：N = A·L·ε_s / ((4/3)πR³)
    auto electrode_mol = [&](const Vector& c, const std::vector<real_t>& volumes,
                             real_t L, real_t eps_s, real_t R) {
        const real_t particle_volume = (4.0 / 3.0) * kPi * R * R * R;
        const real_t n_particles = area_ * L * eps_s / particle_volume;
        real_t mol_per_particle = 0.0;
        for (int i = 0; i < n_; ++i) {
            mol_per_particle += volumes[i] * c(i);
        }
        return n_particles * mol_per_particle;
    };

    return electrode_mol(c_s_neg_, vol_neg_, params_.neg_thickness,
                         params_.neg_active_fraction, params_.neg_particle_radius)
         + electrode_mol(c_s_pos_, vol_pos_, params_.pos_thickness,
                         params_.pos_active_fraction, params_.pos_particle_radius);
}

SolutionSnapshot SPMModel::snapshot(real_t current) const {
    SolutionSnapshot snap;
    snap.time = time_;
    snap.voltage = voltage_;
    snap.current = current;
    snap.temperature = params_.temperature;
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
//  DAE 视图（SUNDIALS IDA）：纯 ODE，状态 = [c_s_n, c_s_p]
//  ODE 右端：c' = D·c − e_surf·N（表面控制体通量注入，同 CN 源项）
// ══════════════════════════════════════════════════════════════

int SPMModel::dae_size() const {
    return 2 * n_;
}

void SPMModel::dae_mask(Vector& id) const {
    id.setOnes(2 * n_); // 全部为微分变量
}

void SPMModel::dae_atol(Vector& atol) const {
    // 浓度尺度 ~3e4/6e4 mol/m³：0.5 mol/m³ ≈ 1e-5 相对精度
    atol.resize(2 * n_);
    atol.setConstant(0.5);
}

void SPMModel::dae_pack_state(Vector& y, Vector& yp) const {
    y.resize(2 * n_);
    yp.resize(2 * n_);
    y.segment(0, n_) = c_s_neg_;
    y.segment(n_, n_) = c_s_pos_;

    const real_t N_neg = dae_current_
                       / (area_ * a_neg_ * params_.neg_thickness) / kFaraday;
    const real_t N_pos = -dae_current_
                        / (area_ * a_pos_ * params_.pos_thickness) / kFaraday;
    yp.segment(0, n_) = D_neg_ * c_s_neg_;
    yp(n_ - 1) -= cn_->surf_area_neg / cn_->surf_vol_neg * N_neg;
    yp.segment(n_, n_) = D_pos_ * c_s_pos_;
    yp(2 * n_ - 1) -= cn_->surf_area_pos / cn_->surf_vol_pos * N_pos;
}

void SPMModel::dae_residual(real_t t, const Eigen::Ref<const Vector>& y,
                                const Eigen::Ref<const Vector>& yp,
                                Eigen::Ref<Vector> res) {
    (void)t;
    const real_t N_neg = dae_current_
                       / (area_ * a_neg_ * params_.neg_thickness) / kFaraday;
    const real_t N_pos = -dae_current_
                        / (area_ * a_pos_ * params_.pos_thickness) / kFaraday;

    Eigen::Map<const Vector> c_n(y.data(), n_);
    Eigen::Map<const Vector> c_p(y.data() + n_, n_);
    Eigen::Map<const Vector> d_n(yp.data(), n_);
    Eigen::Map<const Vector> d_p(yp.data() + n_, n_);
    Eigen::Map<Vector> r_n(res.data(), n_);
    Eigen::Map<Vector> r_p(res.data() + n_, n_);

    // 微分行：y' − (D·c − e_surf·N)
    r_n = d_n - (D_neg_ * c_n);
    r_n(n_ - 1) += cn_->surf_area_neg / cn_->surf_vol_neg * N_neg;
    r_p = d_p - (D_pos_ * c_p);
    r_p(n_ - 1) += cn_->surf_area_pos / cn_->surf_vol_pos * N_pos;
}

void SPMModel::dae_load_state(real_t t, const Eigen::Ref<const Vector>& y) {
    Eigen::Map<const Vector> c_n(y.data(), n_);
    Eigen::Map<const Vector> c_p(y.data() + n_, n_);
    c_s_neg_ = c_n.cwiseMax(0.0).cwiseMin(params_.neg_max_conc);
    c_s_pos_ = c_p.cwiseMax(0.0).cwiseMin(params_.pos_max_conc);
    time_ = t;
    update_voltage(dae_current_);
}

} // namespace fluxcell::model
