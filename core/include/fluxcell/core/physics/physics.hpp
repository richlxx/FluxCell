/// @file physics.hpp
/// @brief 电化学本构关系：Butler-Volmer 动力学、热力学（OCP）、传输性质
#pragma once

#include <fluxcell/core/types.hpp>

#include <cmath>
#include <functional>
#include <vector>

namespace fluxcell::physics {

// ══════════════════════════════════════════════════════════════
//  Butler-Volmer 动力学
// ══════════════════════════════════════════════════════════════

/// Butler-Volmer 方程：计算局部反应电流密度
/// j = k * c_e^α_a * (c_s_max - c_s_surf)^α_a * c_s_surf^α_c
///     * [exp(α_a * F / (R*T) * η) - exp(-α_c * F / (R*T) * η)]
///
/// @param eta        过电位 η = φ_s - φ_e - U(c_s_surf) [V]
/// @param c_e        电解液浓度 [mol/m³]
/// @param c_s_surf   固相表面浓度 [mol/m³]
/// @param c_s_max    固相最大浓度 [mol/m³]
/// @param k          反应速率常数 [m/(mol^0.5·s)]
/// @param alpha_a    阳极传递系数 [-]
/// @param alpha_c    阴极传递系数 [-]
/// @param T          温度 [K]
/// @return 局部反应电流密度 j [A/m²]
real_t butler_volmer(real_t eta, real_t c_e, real_t c_s_surf, real_t c_s_max,
                     real_t k, real_t alpha_a, real_t alpha_c, real_t T);

/// Butler-Volmer 对过电位的偏导数 dj/dη
real_t butler_volmer_deta(real_t eta, real_t c_e, real_t c_s_surf, real_t c_s_max,
                          real_t k, real_t alpha_a, real_t alpha_c, real_t T);

// ══════════════════════════════════════════════════════════════
//  热力学：开路电位 (OCP)
// ══════════════════════════════════════════════════════════════

/// OCP 函数类型：U(stoichiometry) → [V]
using OCPFunction = std::function<real_t(real_t stoich)>;

/// OCP 查表线性插值：stoich/voltage 等长采样序列（升序），端点外钳位。
/// @return 插值电压 [V]；表为空或两序列长度不一致时返回 NaN
real_t ocp_from_lut(const std::vector<real_t>& stoich,
                    const std::vector<real_t>& voltage, real_t x);

/// 石墨负极 OCP（基于 stoichiometry x = c_s / c_s_max）
/// 参考：Safari et al., J. Electrochem. Soc. 156, A61 (2009)
real_t graphite_ocp(real_t x);

/// NMC111 正极 OCP
/// 参考：Lu et al., J. Electrochem. Soc. 160, A610 (2013)
real_t nmc111_ocp(real_t y);

/// LFP 正极 OCP（平坦平台 + 缓变）
real_t lfp_ocp(real_t y);

/// 石墨负极 OCP — Chen2020 (LG M50)
/// 参考：Chen et al., J. Electrochem. Soc. 167, 080534 (2020)
real_t graphite_chen2020_ocp(real_t x);

/// NMC 正极 OCP — Chen2020 (LG M50)
/// 参考：Chen et al., J. Electrochem. Soc. 167, 080534 (2020)
real_t nmc_chen2020_ocp(real_t y);

/// Chen2020 负极交换电流密度（石墨/LiPF6 EC:DMC）
/// i0 = m_ref * exp(E_r/R (1/298.15 - 1/T)) * sqrt(c_e · c_s_surf · (c_s_max - c_s_surf))
/// m_ref = 6.48e-7 (A/m²)(m³/mol)^1.5, E_r = 35000 J/mol
real_t exchange_current_density_chen2020_neg(real_t c_e, real_t c_s_surf,
                                             real_t c_s_max, real_t T);

/// Chen2020 正极交换电流密度（NMC/LiPF6 EC:DMC）
/// m_ref = 3.42e-6 (A/m²)(m³/mol)^1.5, E_r = 17800 J/mol
real_t exchange_current_density_chen2020_pos(real_t c_e, real_t c_s_surf,
                                             real_t c_s_max, real_t T);

/// 对称 Butler-Volmer 反解过电位（α_a = α_c = 0.5）
/// η = (2RT/F) · asinh(i / (2·i0))
/// @param i   局部电流密度 [A/m²]（正=阳极/氧化）
/// @param i0  交换电流密度 [A/m²]
/// @param T   温度 [K]
/// @return 过电位 η [V]
real_t overpotential_from_bv(real_t i, real_t i0, real_t T);

// ══════════════════════════════════════════════════════════════
//  通用查表插值
// ══════════════════════════════════════════════════════════════

/// 查表线性插值：升序节点、端点钳位；空表或长度不一致返回 NaN
/// （OCP 曲线与电解液输运表共用）
real_t lut_interpolate(const std::vector<real_t>& nodes,
                       const std::vector<real_t>& values, real_t x);

// ══════════════════════════════════════════════════════════════
//  Nyman2008 电解液输运（PyBaMM Chen2020 默认，LiPF6 EC:DMC）
//  用于标定浓度表采样与测试对照；模型内部经 YAML 表 + 线性插值使用
// ══════════════════════════════════════════════════════════════

/// 离子电导率 κ(c_e) [S/m]：κ = 0.1297·m³ − 2.51·m^1.5 + 3.329·m，m = c_e/1000
real_t nyman2008_conductivity(real_t c_e);

/// 扩散系数 D(c_e) [m²/s]：D = 8.794e-11·m² − 3.972e-10·m + 4.862e-10
real_t nyman2008_diffusivity(real_t c_e);

// ══════════════════════════════════════════════════════════════
//  电解液传输性质
// ══════════════════════════════════════════════════════════════

/// 电解液离子电导率 κ(c_e, T) [S/m]
/// 默认模型：简化常数模型，可替换为浓度/温度依赖函数
real_t electrolyte_conductivity(real_t c_e, real_t T);

/// 电解液扩散系数 D_e(c_e, T) [m²/s]
real_t electrolyte_diffusivity(real_t c_e, real_t T);

/// 电解液热力学因子 (1 + d ln f / d ln c_e)
real_t thermodynamic_factor(real_t c_e, real_t T);

} // namespace fluxcell::physics
