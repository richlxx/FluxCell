/// @file physics.cpp
/// @brief 电化学本构关系实现
#include <fluxcell/core/physics/physics.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

namespace fluxcell::physics {

// ══════════════════════════════════════════════════════════════
//  常量
// ══════════════════════════════════════════════════════════════
static constexpr real_t R_GAS = 8.314;     // J/(mol·K)
static constexpr real_t FARADAY = 96485.0; // C/mol

// ══════════════════════════════════════════════════════════════
//  Butler-Volmer 动力学
// ══════════════════════════════════════════════════════════════

real_t butler_volmer(real_t eta, real_t c_e, real_t c_s_surf, real_t c_s_max,
                     real_t k, real_t alpha_a, real_t alpha_c, real_t T) {
    // 交换电流密度：i₀ = k * c_e^α_a * (c_s_max - c_s_surf)^α_a * c_s_surf^α_c
    real_t c_s_surf_clamped = std::clamp(c_s_surf, 1.0, c_s_max - 1.0);
    real_t i0 = k * std::pow(c_e, alpha_a)
                  * std::pow(c_s_max - c_s_surf_clamped, alpha_a)
                  * std::pow(c_s_surf_clamped, alpha_c);

    real_t FRT = FARADAY / (R_GAS * T);

    // j = i₀ * [exp(α_a * F/(RT) * η) - exp(-α_c * F/(RT) * η)]
    real_t j = i0 * (std::exp(alpha_a * FRT * eta) - std::exp(-alpha_c * FRT * eta));
    return j;
}

real_t butler_volmer_deta(real_t eta, real_t c_e, real_t c_s_surf, real_t c_s_max,
                          real_t k, real_t alpha_a, real_t alpha_c, real_t T) {
    real_t c_s_surf_clamped = std::clamp(c_s_surf, 1.0, c_s_max - 1.0);
    real_t i0 = k * std::pow(c_e, alpha_a)
                  * std::pow(c_s_max - c_s_surf_clamped, alpha_a)
                  * std::pow(c_s_surf_clamped, alpha_c);

    real_t FRT = FARADAY / (R_GAS * T);

    return i0 * FRT * (alpha_a * std::exp(alpha_a * FRT * eta)
                      + alpha_c * std::exp(-alpha_c * FRT * eta));
}

// ══════════════════════════════════════════════════════════════
//  热力学：开路电位 (OCP)
// ══════════════════════════════════════════════════════════════

real_t graphite_ocp(real_t x) {
    // 石墨 OCP — Safari & Delacourt (2009) 简化拟合
    // x = stoichiometry = c_s / c_s_max
    x = std::clamp(x, 0.001, 0.999);

    return 0.194 + 1.5 * std::exp(-120.0 * x)
         + 0.0351 * std::tanh((x - 0.286) / 0.083)
         - 0.0045 * std::tanh((x - 0.9) / 0.119)
         - 0.035 * std::tanh((x - 0.99) / 0.05)
         - 0.0147 * std::tanh((x - 0.5) / 0.034)
         - 0.102 * std::tanh((x - 0.194) / 0.142)
         - 0.022 * std::tanh((x - 0.9) / 0.0164)
         - 0.011 * std::tanh((x - 0.124) / 0.0226)
         + 0.0155 * std::tanh((x - 0.105) / 0.029);
}

real_t nmc111_ocp(real_t y) {
    // NMC111 正极 OCP — Lu et al. (2013) 拟合
    y = std::clamp(y, 0.001, 0.999);

    return -0.8090 * y + 4.4875
         - 0.0428 * std::tanh(18.5138 * (y - 0.5542))
         - 17.7326 * std::tanh(15.7890 * (y - 0.3117))
         + 17.5842 * std::tanh(15.9308 * (y - 0.3120));
}

real_t lfp_ocp(real_t y) {
    // LFP 正极 OCP — 典型两相平台
    y = std::clamp(y, 0.001, 0.999);

    return 3.42
         - 0.05 * std::tanh(40.0 * (y - 0.05))
         + 0.05 * std::tanh(40.0 * (y - 0.95))
         + 0.1  * (1.0 / (1.0 + std::exp(-40.0 * (y - 0.1))) - 1.0 / (1.0 + std::exp(-40.0 * (y - 0.9))));
}

// ══════════════════════════════════════════════════════════════
//  Chen2020 (LG M50) 热力学与动力学
// ══════════════════════════════════════════════════════════════

real_t graphite_chen2020_ocp(real_t x) {
    // 石墨负极 OCP — Chen et al. (2020), PyBaMM graphite_LGM50_ocp_Chen2020
    x = std::clamp(x, 1e-6, 1.0 - 1e-6);

    return 1.9793 * std::exp(-39.3631 * x) + 0.2482
         - 0.0909 * std::tanh(29.8538 * (x - 0.1234))
         - 0.04478 * std::tanh(14.9159 * (x - 0.2769))
         - 0.0205 * std::tanh(30.4444 * (x - 0.6103));
}

real_t nmc_chen2020_ocp(real_t y) {
    // NMC 正极 OCP — Chen et al. (2020), PyBaMM nmc_LGM50_ocp_Chen2020
    y = std::clamp(y, 1e-6, 1.0 - 1e-6);

    return -0.8090 * y + 4.4875
         - 0.0428 * std::tanh(18.5138 * (y - 0.5542))
         - 17.7326 * std::tanh(15.7890 * (y - 0.3117))
         + 17.5842 * std::tanh(15.9308 * (y - 0.3120));
}

real_t ocp_from_lut(const std::vector<real_t>& stoich,
                    const std::vector<real_t>& voltage, real_t x) {
    return lut_interpolate(stoich, voltage, x);
}

real_t exchange_current_density_chen2020_neg(real_t c_e, real_t c_s_surf,
                                             real_t c_s_max, real_t T) {
    // i0 = m_ref · exp(E_r/R · (1/298.15 - 1/T)) · √(c_e · c_s · (c_s_max - c_s))
    static constexpr real_t m_ref = 6.48e-7; // (A/m²)(m³/mol)^1.5
    static constexpr real_t E_r = 35000.0;   // J/mol

    // 两端各留 1 mol/m³ 余量（与 butler_volmer 的 clamp 一致）：饱和端
    // c_surf 精确贴 c_max 时 i0 严格为 0，反解过电位 asinh(i/2i0) → ∞，
    // 深过充 CV→静置→放电切换的 Gummel 重估会直接 sinh 溢出；
    // 保留余量后 i0 ≥ k·√(c_e·1) 量级，反解有限且连续
    real_t c_s_clamped = std::clamp(c_s_surf, 0.0, c_s_max - 1.0);
    real_t prod = std::max(c_e * c_s_clamped * (c_s_max - c_s_clamped), 0.0);
    return m_ref * std::exp(E_r / R_GAS * (1.0 / 298.15 - 1.0 / T)) * std::sqrt(prod);
}

real_t exchange_current_density_chen2020_pos(real_t c_e, real_t c_s_surf,
                                             real_t c_s_max, real_t T) {
    static constexpr real_t m_ref = 3.42e-6; // (A/m²)(m³/mol)^1.5
    static constexpr real_t E_r = 17800.0;   // J/mol

    // 同负极：饱和端留 1 mol/m³ 余量，保证 i0 > 0（见 neg 版注释）
    real_t c_s_clamped = std::clamp(c_s_surf, 0.0, c_s_max - 1.0);
    real_t prod = std::max(c_e * c_s_clamped * (c_s_max - c_s_clamped), 0.0);
    return m_ref * std::exp(E_r / R_GAS * (1.0 / 298.15 - 1.0 / T)) * std::sqrt(prod);
}

real_t overpotential_from_bv(real_t i, real_t i0, real_t T) {
    // 对称 BV（α_a = α_c = 0.5）解析反解：η = (2RT/F) · asinh(i / (2·i0))
    if (i0 <= 0.0) {
        return std::numeric_limits<real_t>::quiet_NaN();
    }
    return (2.0 * R_GAS * T / FARADAY) * std::asinh(i / (2.0 * i0));
}

// ══════════════════════════════════════════════════════════════
//  通用查表插值
// ══════════════════════════════════════════════════════════════

real_t lut_interpolate(const std::vector<real_t>& nodes,
                       const std::vector<real_t>& values, real_t x) {
    if (nodes.empty() || values.size() != nodes.size()) {
        return std::numeric_limits<real_t>::quiet_NaN();
    }
    if (x <= nodes.front()) {
        return values.front(); // 端点外钳位
    }
    if (x >= nodes.back()) {
        return values.back();
    }
    const auto it = std::lower_bound(nodes.begin(), nodes.end(), x);
    const int i = static_cast<int>(it - nodes.begin());
    const real_t x0 = nodes[static_cast<std::size_t>(i) - 1];
    const real_t x1 = nodes[static_cast<std::size_t>(i)];
    if (x1 <= x0) {
        return values[static_cast<std::size_t>(i)]; // 重复采样点防护
    }
    const real_t v0 = values[static_cast<std::size_t>(i) - 1];
    const real_t v1 = values[static_cast<std::size_t>(i)];
    return v0 + (v1 - v0) * (x - x0) / (x1 - x0);
}

// ══════════════════════════════════════════════════════════════
//  Nyman2008 电解液输运（PyBaMM Chen2020 默认）
// ══════════════════════════════════════════════════════════════

real_t nyman2008_conductivity(real_t c_e) {
    const real_t m = c_e / 1000.0;
    return 0.1297 * m * m * m - 2.51 * std::pow(m, 1.5) + 3.329 * m;
}

real_t nyman2008_diffusivity(real_t c_e) {
    const real_t m = c_e / 1000.0;
    return 8.794e-11 * m * m - 3.972e-10 * m + 4.862e-10;
}

// ══════════════════════════════════════════════════════════════
//  电解液传输性质
// ══════════════════════════════════════════════════════════════

real_t electrolyte_conductivity(real_t c_e, real_t T) {
    // 简化模型：κ = κ_ref * (c_e / c_ref)^α * exp(E_a/R * (1/T_ref - 1/T))
    // 第一阶段使用常数近似
    (void)c_e;
    (void)T;
    return 1.0; // [S/m] TODO: 替换为浓度依赖模型
}

real_t electrolyte_diffusivity(real_t c_e, real_t T) {
    (void)c_e;
    (void)T;
    return 7.5e-10; // [m²/s] TODO: 替换为浓度/温度依赖模型
}

real_t thermodynamic_factor(real_t c_e, real_t T) {
    (void)c_e;
    (void)T;
    return 1.0; // 理想溶液近似
}

} // namespace fluxcell::physics
