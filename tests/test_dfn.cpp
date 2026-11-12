/// @file test_dfn.cpp
/// @brief 全阶 DFN (P2D) 模型单元测试：静置极限 / 电流守恒 / 锂守恒 /
///        反应分布 / 与 SPMe 的电压对比 / 快照字段 / 引擎接入 / 参数 IO /
///        SEI 降解（机理寿命）/ 循环工况
#include <fluxcell/core/model/dfn_model.hpp>
#include <fluxcell/core/model/spme_model.hpp>
#include <fluxcell/core/physics/physics.hpp>
#include <fluxcell/core/simulation/simulation.hpp>
#include <fluxcell/core/io/io.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <stdexcept>

namespace {

using fluxcell::real_t;
using fluxcell::SolutionSnapshot;
using fluxcell::SimulationEngine;
using fluxcell::SimulationStatus;
using fluxcell::CellModelType;
using fluxcell::CCStep;
using fluxcell::MultiStepProtocol;
using fluxcell::model::DFNModel;
using fluxcell::model::SPMeModel;
using fluxcell::model::ThermalType;

/// Chen2020 默认参数（DFN 形式）
DFNModel::Parameters chen2020_dfn() {
    DFNModel::Parameters p; // 默认值即 Chen2020
    return p;
}

TEST(DFNModel, InitialVoltageEqualsPureOCPDifference) {
    auto params = chen2020_dfn();
    const auto& solid = params.solid;
    const real_t x0 = solid.neg_stoich_at_soc0
                    + solid.initial_soc
                          * (solid.neg_stoich_at_soc1 - solid.neg_stoich_at_soc0);
    const real_t y0 = solid.pos_stoich_at_soc0
                    + solid.initial_soc
                          * (solid.pos_stoich_at_soc1 - solid.pos_stoich_at_soc0);
    const real_t ocp_diff =
        fluxcell::physics::nmc_chen2020_ocp(y0)
        - fluxcell::physics::graphite_chen2020_ocp(x0);

    // 零电流构造：均匀浓度 + 电流守恒 → 逐格 η = 0，无任何损失项
    DFNModel model(params, 0.0);
    EXPECT_NEAR(model.voltage(), ocp_diff, 1e-9);

    // 逐格过电位全为零（静置极限）
    for (int k = 0; k < model.neg_reaction_current().size(); ++k) {
        EXPECT_NEAR(model.neg_reaction_current()(k), 0.0, 1e-12);
    }
    for (int k = 0; k < model.pos_reaction_current().size(); ++k) {
        EXPECT_NEAR(model.pos_reaction_current()(k), 0.0, 1e-12);
    }
}

TEST(DFNModel, RestPreservesStateAndVoltage) {
    DFNModel model(chen2020_dfn(), 0.0);
    const real_t v0 = model.voltage();

    for (int i = 0; i < 5; ++i) {
        model.step(0.0, 60.0);
    }

    // 静置无源：c_e 保持均匀初值，φ_e 保持 0，电压不变
    const auto& c_e = model.electrolyte_conc();
    ASSERT_EQ(c_e.size(), 100);
    for (int i = 0; i < c_e.size(); ++i) {
        EXPECT_NEAR(c_e(i), 1000.0, 1e-9);
        EXPECT_NEAR(model.electrolyte_potential()(i), 0.0, 1e-12);
    }
    EXPECT_NEAR(model.electrolyte_overpotential(), 0.0, 1e-12);
    EXPECT_NEAR(model.voltage(), v0, 1e-12);
}

TEST(DFNModel, CurrentConservedOnBothElectrodes) {
    DFNModel model(chen2020_dfn(), 5.0);
    for (int i = 0; i < 5; ++i) {
        model.step(5.0, 10.0);
    }

    // Σ a·dx·i_surf = i_app（Newton 电流守恒残差 ≤ tolerance × ref）
    const auto& i_n = model.neg_reaction_current();
    const auto& i_p = model.pos_reaction_current();
    real_t sum_n = 0.0;
    real_t sum_p = 0.0;
    for (int k = 0; k < i_n.size(); ++k) {
        sum_n += i_n(k);
    }
    for (int k = 0; k < i_p.size(); ++k) {
        sum_p += i_p(k);
    }
    const auto& p = model.parameters();
    const real_t dx = (p.solid.neg_thickness + p.solid.sep_thickness
                       + p.solid.pos_thickness)
                    / static_cast<real_t>(p.electrolyte_cells);
    const real_t a_n = 3.0 * p.solid.neg_active_fraction / p.solid.neg_particle_radius;
    const real_t a_p = 3.0 * p.solid.pos_active_fraction / p.solid.pos_particle_radius;

    // Σ a·dx·i_surf = i_app = I/A（面电流密度语义）
    const real_t area = p.solid.electrode_height * p.solid.electrode_width;
    const real_t i_app = 5.0 / area;
    EXPECT_NEAR(sum_n * a_n * dx, i_app, 1e-6);
    EXPECT_NEAR(sum_p * a_p * dx, -i_app, 1e-6);
}

TEST(DFNModel, TotalLithiumConservedUnderCycling) {
    DFNModel model(chen2020_dfn(), 5.0);
    const real_t total0 = model.total_lithium_mol();
    EXPECT_GT(total0, 0.1); // 5 Ah 电池固相 + 液相锂约 0.2 mol

    // 放电 + 充电 + 静置循环（液相源逐格但总量守恒，固相零通量）
    for (int i = 0; i < 10; ++i) {
        model.step(5.0, 10.0);
    }
    for (int i = 0; i < 5; ++i) {
        model.step(-5.0, 10.0);
    }
    model.step(0.0, 60.0);

    EXPECT_NEAR(model.total_lithium_mol(), total0, 1e-9 * total0);
}

TEST(DFNModel, DischargeBuildsElectrolyteGradient) {
    DFNModel model(chen2020_dfn(), 5.0);
    for (int i = 0; i < 5; ++i) {
        model.step(5.0, 10.0);
    }

    // 与 SPMe 相同的液相趋势：负极区浓度升高、正极区降低
    const auto& c_e = model.electrolyte_conc();
    EXPECT_GT(c_e(0), 1000.0);
    EXPECT_LT(c_e(c_e.size() - 1), 1000.0);
}

TEST(DFNModel, ReactionDistributionVariesAlongElectrode) {
    DFNModel model(chen2020_dfn(), 5.0);
    model.step(5.0, 10.0);

    // 电位损失沿厚度变化 → η 分布不均 → 孔壁电流非均匀 → 各列表面浓度分化
    const auto& c_surf = model.neg_surface_conc_profile();
    EXPECT_GT(c_surf.maxCoeff() - c_surf.minCoeff(), 1.0); // > 1 mol/m³ 差异

    // 放电耗锂：平均表面浓度低于初始值（100% SOC 石墨 x ≈ 0.9）
    const auto& p = model.parameters().solid;
    EXPECT_LT(c_surf.mean(), p.initial_soc * p.neg_max_conc);

    // 放电方向：负极阳极极化 η > 0，正极阴极极化 η < 0（快照平均）
    const SolutionSnapshot snap = model.snapshot(5.0);
    EXPECT_GT(snap.eta_neg, 0.0);
    EXPECT_LT(snap.eta_pos, 0.0);
}

TEST(DFNModel, SolidOhmicLossReducesVoltage) {
    // 对照同一 DFN 仅改固相电导：分布效应相消，差异 ≈ 纯固相欧姆。
    // 正极 σ = 0.18 S/m（Chen2020）@1C → Δφ_s ≈ 19 mV。
    auto p_low = chen2020_dfn();
    auto p_high = chen2020_dfn();
    p_high.neg_conductivity = 1.0e6;
    p_high.pos_conductivity = 1.0e6;

    DFNModel m_low(p_low, 5.0);
    DFNModel m_high(p_high, 5.0);
    for (int i = 0; i < 10; ++i) {
        m_low.step(5.0, 10.0);
        m_high.step(5.0, 10.0);
    }

    const real_t dv = m_high.voltage() - m_low.voltage();
    EXPECT_GT(dv, 0.010);  // 固相欧姆显著
    EXPECT_LT(dv, 0.030);  // 但在几十 mV 量级
}

TEST(DFNModel, SmallCurrentApproachesOCP) {
    auto params = chen2020_dfn();
    DFNModel loaded(params, 5.0);
    DFNModel gentle(params, 0.25);

    for (int i = 0; i < 10; ++i) {
        loaded.step(5.0, 10.0);
        gentle.step(0.25, 10.0);
    }

    const real_t v_ocp = DFNModel(params, 0.0).voltage();
    const real_t loss_1c = std::abs(v_ocp - loaded.voltage());
    const real_t loss_c20 = std::abs(v_ocp - gentle.voltage());
    EXPECT_GT(loss_1c, 0.0);
    // C/20 活化 ~9 mV + 欧姆 ~2 mV，1C 活化 ~100 mV + 欧姆 ~40 mV
    EXPECT_LT(loss_c20, 0.3 * loss_1c);
}

TEST(DFNModel, RestRelaxesVoltageTowardOCP) {
    DFNModel model(chen2020_dfn(), 5.0);
    model.step(5.0, 10.0);
    const real_t v_loaded = model.voltage();

    for (int i = 0; i < 5; ++i) {
        model.step(0.0, 60.0);
    }

    // 静置弛豫：欧姆/活化损失消失，电压向 OCP 回升
    EXPECT_GT(model.voltage(), v_loaded);
}

TEST(DFNModel, SnapshotCarriesPotentialFields) {
    DFNModel model(chen2020_dfn(), 5.0);
    model.step(5.0, 10.0);

    const SolutionSnapshot snap = model.snapshot(5.0);
    const int m = static_cast<int>(model.electrolyte_potential().size());
    EXPECT_EQ(snap.c_e.size(), m);
    EXPECT_EQ(snap.phi_e.size(), m);
    ASSERT_EQ(snap.phi_s.size(), m); // DFN 提供全厚度固相电位
    EXPECT_EQ(snap.c_s_neg.size(), 20); // 颗粒径向网格（平均颗粒）
    EXPECT_EQ(snap.c_s_pos.size(), 20);

    const auto& p = model.parameters().solid;
    const real_t dx = (p.neg_thickness + p.sep_thickness + p.pos_thickness)
                    / static_cast<real_t>(m);
    const int sep_start = static_cast<int>(std::lround(p.neg_thickness / dx));
    const int pos_start = static_cast<int>(std::lround(
        (p.neg_thickness + p.sep_thickness) / dx));
    // 隔膜区固相电位为负极末端延伸（常数）
    for (int i = sep_start; i < pos_start; ++i) {
        EXPECT_NEAR(snap.phi_s(i), snap.phi_s(sep_start - 1), 1e-12);
    }
    // 全部场变量有限
    for (int i = 0; i < m; ++i) {
        EXPECT_TRUE(std::isfinite(snap.phi_e(i)));
        EXPECT_TRUE(std::isfinite(snap.phi_s(i)));
        EXPECT_TRUE(std::isfinite(snap.c_e(i)));
    }
}

// ══════════════════════════════════════════════════════════════
//  耦合热模型（lumped，参考 PyBaMM options {"thermal": "lumped"}）
// ══════════════════════════════════════════════════════════════

TEST(DFNModel, LumpedRestKeepsAmbientTemperature) {
    auto p = chen2020_dfn();
    p.thermal_type = ThermalType::Lumped;
    p.thermal_params.t_ambient = 303.15;
    p.thermal_params.t_initial = 303.15;

    DFNModel model(p, 0.0);
    for (int i = 0; i < 10; ++i) {
        model.step(0.0, 60.0);
    }

    // 静置：无产热、无温差 → 温度精确不变，产热为 0
    EXPECT_NEAR(model.temperature(), 303.15, 1e-9);
    EXPECT_NEAR(model.total_heating(), 0.0, 1e-12);
}

TEST(DFNModel, IsothermalTemperatureStaysConstant) {
    auto p = chen2020_dfn(); // 默认等温
    p.solid.temperature = 305.0;
    DFNModel model(p, 5.0);
    for (int i = 0; i < 10; ++i) {
        model.step(5.0, 10.0);
    }
    // 等温模式温度恒为参数温度（与热耦合前行为一致）
    EXPECT_NEAR(model.temperature(), 305.0, 1e-12);
}

TEST(DFNModel, LumpedAdiabaticHeatingConservesEnergy) {
    auto p = chen2020_dfn();
    p.thermal_type = ThermalType::Lumped;
    p.thermal_params.h_total = 0.0;    // 绝热：无散热
    p.thermal_params.cooling_area = 0.0;

    // 有效热容（与模型 rho_c_p_eff 同口径：各层按厚度加权，含集流体）
    const auto& th = p.thermal_params;
    const auto& s = p.solid;
    const real_t l_total = s.neg_thickness + s.sep_thickness + s.pos_thickness
                         + th.cc_thickness_neg + th.cc_thickness_pos;
    const real_t rho_cp_eff =
        (th.rho_c_p_cc_neg * th.cc_thickness_neg
         + th.rho_c_p_neg * s.neg_thickness
         + th.rho_c_p_sep * s.sep_thickness
         + th.rho_c_p_pos * s.pos_thickness
         + th.rho_c_p_cc_pos * th.cc_thickness_pos) / l_total;
    const real_t c_th = rho_cp_eff * th.cell_volume;

    DFNModel model(p, 5.0); // 1C 恒流放电
    const real_t t0 = model.temperature();

    // 对产热做梯形积分（CN 温度更新同样用 (q_old+q_new)/2，口径一致）
    real_t q_prev = model.total_heating();
    real_t heat_integral = 0.0;
    for (int i = 0; i < 60; ++i) {
        model.step(5.0, 10.0);
        const real_t q_now = model.total_heating();
        heat_integral += 0.5 * (q_prev + q_now) * 10.0;
        q_prev = q_now;
    }

    // 能量守恒：C_th·ΔT = ∫Q dt（二阶 CN + 梯形，残差为高阶项）
    const real_t expected_dt = heat_integral / c_th;
    const real_t actual_dt = model.temperature() - t0;
    EXPECT_GT(actual_dt, 1.0); // 1C 放电 10 min 绝热升温明显
    EXPECT_NEAR(actual_dt, expected_dt, 0.02 * std::abs(expected_dt));
}

TEST(DFNModel, LumpedDischargeHeatsAndNewtonCoolingLimitsRise) {
    auto p_cool = chen2020_dfn();
    p_cool.thermal_type = ThermalType::Lumped;
    p_cool.thermal_params.h_total = 10.0;
    p_cool.thermal_params.cooling_area = 0.00531;

    auto p_adia = p_cool;
    p_adia.thermal_params.h_total = 0.0;
    p_adia.thermal_params.cooling_area = 0.0;

    DFNModel m_cool(p_cool, 5.0);
    DFNModel m_adia(p_adia, 5.0);
    for (int i = 0; i < 60; ++i) {
        m_cool.step(5.0, 10.0);
        m_adia.step(5.0, 10.0);
    }

    // 1C 放电：产热为正、温度上升；散热使温升低于绝热
    EXPECT_GT(m_cool.total_heating(), 0.0);
    EXPECT_GT(m_cool.temperature(), 298.15);
    EXPECT_LT(m_cool.temperature(), m_adia.temperature());
    // 产热三分量自洽：Q_total = Q_ohm + Q_rxn + Q_rev
    EXPECT_GE(m_cool.ohmic_heating(), 0.0);
    EXPECT_GE(m_cool.reaction_heating(), 0.0);
    EXPECT_NEAR(m_cool.total_heating(),
                m_cool.ohmic_heating() + m_cool.reaction_heating()
                    + m_cool.reversible_heating(),
                1e-12);
}

TEST(EntropyLut, PiecewiseLinearInterpolation) {
    using fluxcell::model::EntropyLut;
    EntropyLut lut;
    lut.sto = {0.2, 0.6, 0.8};
    lut.dudt = {1.0e-4, 3.0e-4, -2.0e-4};

    // 区间内分段线性（中点取两端均值）
    EXPECT_NEAR(lut.eval(0.4), 2.0e-4, 1e-12);
    EXPECT_NEAR(lut.eval(0.7), 0.5e-4, 1e-12);
    EXPECT_NEAR(lut.eval(0.6), 3.0e-4, 1e-12); // 节点处取节点值
    // 范围外 clamp 端点（与 numpy.interp 一致）
    EXPECT_NEAR(lut.eval(0.0), 1.0e-4, 1e-12);
    EXPECT_NEAR(lut.eval(1.0), -2.0e-4, 1e-12);
    // 空表恒 0（未配置熵变系数时的语义）
    EXPECT_EQ(EntropyLut{}.eval(0.5), 0.0);
}

TEST(DFNModel, EntropicLutProducesAnalyticReversibleHeat) {
    // 恒定 dU/dT LUT 下可逆熵热有解析式：q_rev = I·T·(dU_n/dT − dU_p/dT)，
    // （放电时负极 j_n > 0、正极 j_p < 0，符号相反；电流守恒
    // Σ a·dx·i_surf,n = I/A 精确成立，因此期望值不含数值近似）
    auto p = chen2020_dfn();
    p.thermal_type = ThermalType::Lumped;
    constexpr real_t kDudtN = 2.0e-4;
    constexpr real_t kDudtP = 1.0e-4;
    p.thermal_params.entropic_lut_neg = {{0.0, 1.0}, {kDudtN, kDudtN}};
    p.thermal_params.entropic_lut_pos = {{0.0, 1.0}, {kDudtP, kDudtP}};
    DFNModel model(p);
    const real_t t0 = model.temperature();
    model.step(5.0, 0.5); // 1C 放电
    const real_t expected = 5.0 * t0 * (kDudtN - kDudtP);
    EXPECT_NEAR(model.reversible_heating(), expected,
                1e-3 * std::abs(expected));
}

TEST(DFNModel, InvalidThermalParametersRejected) {
    auto p = chen2020_dfn();
    p.thermal_type = ThermalType::Lumped;
    p.thermal_params.cell_volume = 0.0;
    EXPECT_THROW(DFNModel{p}, std::invalid_argument);

    p = chen2020_dfn();
    p.thermal_type = ThermalType::Lumped;
    p.thermal_params.t_initial = -1.0;
    EXPECT_THROW(DFNModel{p}, std::invalid_argument);

    p = chen2020_dfn();
    p.thermal_type = ThermalType::Lumped;
    p.thermal_params.rho_c_p_neg = 0.0;
    EXPECT_THROW(DFNModel{p}, std::invalid_argument);
}

TEST(DFNModel, InvalidParametersRejected) {
    auto p = chen2020_dfn();
    p.neg_conductivity = 0.0;
    EXPECT_THROW(DFNModel{p}, std::invalid_argument);

    p = chen2020_dfn();
    p.pos_conductivity = -1.0;
    EXPECT_THROW(DFNModel{p}, std::invalid_argument);

    p = chen2020_dfn();
    p.max_iterations = 0;
    EXPECT_THROW(DFNModel{p}, std::invalid_argument);

    p = chen2020_dfn();
    p.tolerance = 0.0;
    EXPECT_THROW(DFNModel{p}, std::invalid_argument);

    p = chen2020_dfn();
    p.electrolyte_cells = 5;
    EXPECT_THROW(DFNModel{p}, std::invalid_argument);

    p = chen2020_dfn();
    p.solid.particle_cells = 2;
    EXPECT_THROW(DFNModel{p}, std::invalid_argument);
}

TEST(SimulationEngine, DfnSingleStepProducesPotentialData) {
    SimulationEngine engine;
    MultiStepProtocol proto;
    proto.model = CellModelType::DFN;
    proto.time_step = 1.0;
    CCStep step;
    step.mode = CCStep::Mode::Discharge;
    step.current = 5.0;
    step.cutoff_time = 3.0;
    proto.steps.push_back(step);

    engine.single_step({}, proto); // 初始化 + t = 0 快照
    for (int i = 0; i < 3; ++i) {
        engine.single_step({}, proto);
    }
    EXPECT_EQ(engine.status(), SimulationStatus::Completed);

    const auto snaps = engine.collect_snapshots();
    ASSERT_EQ(snaps.size(), 4u); // t = 0..3
    // DFN 快照带固相电位分布；放电段电压下降且液相有损失
    for (const auto& snap : snaps) {
        ASSERT_EQ(snap.phi_s.size(), snap.phi_e.size());
        EXPECT_GT(snap.phi_s.size(), 0);
    }
    EXPECT_LT(snaps.back().voltage, snaps.front().voltage);
    EXPECT_LT(snaps.back().eta_e, 0.0);
}

TEST(Io, LoadDfnParametersReadsConductivity) {
    const std::filesystem::path yml = std::filesystem::path(FLUXCELL_SOURCE_DIR)
                                    / "params" / "dfn_chen2020.yml";
    const auto p = fluxcell::io::load_dfn_parameters(yml);

    // DFN 特有键
    EXPECT_NEAR(p.neg_conductivity, 215.0, 1e-9);
    EXPECT_NEAR(p.pos_conductivity, 0.18, 1e-9);

    // SPMe/固相部分复用既有读取
    EXPECT_NEAR(p.solid.capacity, 5.0, 1e-9);
    EXPECT_NEAR(p.solid.initial_soc, 1.0, 1e-9);
    EXPECT_NEAR(p.solid.neg_thickness, 85.2e-6, 1e-12);
    EXPECT_NEAR(p.electrolyte_diffusivity, 1.7694e-10, 1e-16);
    EXPECT_EQ(p.electrolyte_cells, 100);

    // 浓度输运表随 SPMe 复用读取一并加载（76 个 Nyman2008 节点）
    ASSERT_EQ(p.electrolyte_conc_lut.size(), 76u);
    EXPECT_NEAR(p.conductivity_at(1000.0), 0.9487, 1e-12);
    EXPECT_NEAR(p.diffusivity_at(1000.0), 1.7694e-10, 1e-16);
}

TEST(Io, LoadDfnParametersDefaultsWithoutConductivityKeys) {
    // SPMe 参数集无 conductivity 键 → DFN 特有参数取类内默认值
    const std::filesystem::path yml = std::filesystem::path(FLUXCELL_SOURCE_DIR)
                                    / "params" / "spme_chen2020.yml";
    const auto p = fluxcell::io::load_dfn_parameters(yml);

    EXPECT_NEAR(p.neg_conductivity, 215.0, 1e-12);
    EXPECT_NEAR(p.pos_conductivity, 0.18, 1e-12);
    // 无 degradation 段 → 缺省关闭
    EXPECT_FALSE(p.degradation.enabled);
}

TEST(Io, LoadDfnParametersReadsDegradationSection) {
    // dfn_chen2020.yml 的 degradation 段（enabled: false + OKane2022
    // SEI 参数）应完整读入且与结构默认值一致
    const std::filesystem::path yml = std::filesystem::path(FLUXCELL_SOURCE_DIR)
                                    / "params" / "dfn_chen2020.yml";
    const auto p = fluxcell::io::load_dfn_parameters(yml);

    EXPECT_FALSE(p.degradation.enabled);
    EXPECT_NEAR(p.degradation.l_sei_0, 5.0e-9, 1e-24);
    EXPECT_NEAR(p.degradation.v_bar_sei, 9.585e-05, 1e-14);
    EXPECT_NEAR(p.degradation.z_sei, 1.0, 1e-15);
    EXPECT_NEAR(p.degradation.u_sei, 0.4, 1e-15);
    EXPECT_NEAR(p.degradation.r_sei, 2.0e5, 1e-6);
    EXPECT_NEAR(p.degradation.k_sei, 1.0e-12, 1e-28);
    EXPECT_NEAR(p.degradation.c_ec_0, 4541.0, 1e-9);
    EXPECT_NEAR(p.degradation.d_ec, 2.0e-18, 1e-34);
    EXPECT_NEAR(p.degradation.e_sei, 38000.0, 1e-6);
    EXPECT_NEAR(p.degradation.t_ref, 298.15, 1e-12);
    EXPECT_NEAR(p.degradation.alpha_sei, 0.5, 1e-15);
}

// ══════════════════════════════════════════════════════════════
//  SEI 降解（机理寿命仿真；PyBaMM options {"SEI": "ec reaction limited",
//  "SEI film resistance": "distributed"} 口径，参数取 OKane2022）
// ══════════════════════════════════════════════════════════════

TEST(DFNModel, DegradationDisabledMatchesTinyKineticsBaseline) {
    // 降解关闭（默认）与"启用但速率常数/电阻率取数值零"应给出相同轨迹：
    // r_sei = 0 使膜电阻修正（对主反应 j 生效，与 k_sei 无关）一并关闭；
    // 关闭侧诊断恒 0，开启侧 SEI 恒为初始厚度、锂损失为 0
    auto p_off = chen2020_dfn();
    auto p_tiny = chen2020_dfn();
    p_tiny.degradation.enabled = true;
    p_tiny.degradation.k_sei = 1e-30; // 构造校验要求 > 0；对轨迹影响可忽略
    p_tiny.degradation.r_sei = 0.0;

    DFNModel m_off(p_off, 5.0);
    DFNModel m_tiny(p_tiny, 5.0);
    for (int i = 0; i < 20; ++i) {
        m_off.step(5.0, 10.0);
        m_tiny.step(5.0, 10.0);
    }
    EXPECT_NEAR(m_off.voltage(), m_tiny.voltage(), 1e-12);
    EXPECT_NEAR(m_tiny.sei_thickness_avg(), p_tiny.degradation.l_sei_0, 1e-18);
    EXPECT_NEAR(m_tiny.lithium_loss_ah(), 0.0, 1e-15);
    EXPECT_NEAR(m_tiny.degraded_capacity_ah(), 5.0, 1e-12);
    // 关闭时诊断恒 0（sei_thickness_avg/lithium_loss_ah 直接短路）
    EXPECT_FALSE(m_off.degradation_enabled());
    EXPECT_EQ(m_off.sei_thickness_avg(), 0.0);
    EXPECT_EQ(m_off.lithium_loss_ah(), 0.0);
}

TEST(DFNModel, SeiGrowsChargedNotDischarged) {
    // η_SEI = (φ_s−φ_e) − U_sei 驱动：满电 U_n(x≈0.91) ≈ 0.09 V < U_sei
    // → 还原反应进行，SEI 生长；放空 U_n ≈ 1.12 V > U_sei → 生长停滞
    auto p = chen2020_dfn();
    p.degradation.enabled = true;
    p.degradation.r_sei = 0.0; // 排除膜电阻对 φ 的反馈，聚焦生长方向性

    DFNModel full(p, 0.0);
    for (int i = 0; i < 60; ++i) {
        full.step(0.0, 60.0); // 满电静置 1 h（自放电量级生长）
    }

    auto p_empty = p;
    p_empty.solid.initial_soc = 0.0; // SOC = 0：负极几乎完全脱嵌
    DFNModel empty(p_empty, 0.0);
    for (int i = 0; i < 60; ++i) {
        empty.step(0.0, 60.0);
    }

    const real_t l0 = p.degradation.l_sei_0;
    EXPECT_GT(full.sei_thickness_avg(), l0 + 1e-8);   // 净生长 > 10 nm
    EXPECT_LT(empty.sei_thickness_avg(), l0 + 1e-15); // 生长可忽略
}

TEST(DFNModel, SeiSingleStepMatchesImplicitEulerAnalytic) {
    // 单步小 dt 解析对照：dL/dt = g/(1+qL) 的隐式 Euler 一元二次根。
    // g/q 用模型步末快照的 φ_s−φ_e 复算（advance_sei_cn 冻结同一套
    // 自洽解）；等温 T = t_ref → Arrhenius 因子 = 1。r_sei = 0 使
    // η_SEI 与 i_surf 解耦，公式链（Arrhenius 位置/单位/符号）完全独立
    auto p = chen2020_dfn();
    p.degradation.enabled = true;
    p.degradation.r_sei = 0.0;

    DFNModel model(p, 0.0);
    const real_t dt = 1.0;
    model.step(0.0, dt);

    const auto snap = model.snapshot(0.0);
    const auto& deg = p.degradation;
    const real_t T = model.temperature(); // 等温 = t_ref → arrh = 1
    constexpr real_t kF = 96485.0;        // 与 dfn_model.cpp 一致
    constexpr real_t kR = 8.314;

    const auto& par = p.solid;
    const real_t dx = (par.neg_thickness + par.sep_thickness
                       + par.pos_thickness)
                    / static_cast<real_t>(model.electrolyte_conc().size());
    const int m_neg = static_cast<int>(std::lround(par.neg_thickness / dx));

    real_t l_sum = 0.0;
    for (int k = 0; k < m_neg; ++k) {
        const real_t eta_sei = snap.phi_s(k) - snap.phi_e(k) - deg.u_sei;
        const real_t K = deg.k_sei
            * std::exp(-deg.alpha_sei * kF / (kR * T) * eta_sei);
        const real_t g = deg.v_bar_sei * deg.c_ec_0 * K / deg.z_sei;
        const real_t q = K / deg.d_ec;
        const real_t L = deg.l_sei_0;
        const real_t b = 1.0 - q * L;
        l_sum += (-b + std::sqrt(b * b + 4.0 * q * (L + g * dt))) / (2.0 * q);
    }
    EXPECT_NEAR(model.sei_thickness_avg(), l_sum / m_neg, 1e-9 * l_sum);
}

TEST(DFNModel, LithiumLossMatchesSeiAreaIntegral) {
    // 锂损失 = SEI 增厚的面积分（PyBaMM Q_sei = Δc_sei·V_k·z 口径）：
    // loss = z/V̄·a_n·dx·m_neg·A·(L_avg − L0)·F/3600（按离散网格复算，
    // dx·m_neg ≠ L_n：负极厚度未对齐网格，差 ~0.6%）
    auto p = chen2020_dfn();
    p.degradation.enabled = true;
    DFNModel model(p, 5.0);
    for (int i = 0; i < 50; ++i) {
        model.step(5.0, 10.0);
    }

    const auto& par = model.parameters();
    const real_t a_n = 3.0 * par.solid.neg_active_fraction
                            / par.solid.neg_particle_radius;
    const real_t area = par.solid.electrode_height * par.solid.electrode_width;
    const real_t dx = (par.solid.neg_thickness + par.solid.sep_thickness
                       + par.solid.pos_thickness)
                    / static_cast<real_t>(model.electrolyte_conc().size());
    const int m_neg = static_cast<int>(std::lround(par.solid.neg_thickness / dx));
    const real_t expected = par.degradation.z_sei / par.degradation.v_bar_sei
                          * a_n * dx * m_neg * area
                          * (model.sei_thickness_avg()
                             - par.degradation.l_sei_0)
                          * 96485.0 / 3600.0;

    EXPECT_GT(model.lithium_loss_ah(), 0.0);
    EXPECT_NEAR(model.lithium_loss_ah(), expected, 1e-9 * expected);
    EXPECT_NEAR(model.degraded_capacity_ah(),
                par.solid.capacity - model.lithium_loss_ah(), 1e-12);
}

TEST(DFNModel, SeiFilmResistanceReducesDischargeVoltage) {
    // 分布式膜电阻：BV η 修正 −(j + j_sei)·L·R_sei；L0·R_sei = 1e-3 Ω·m²，
    // 1C 孔壁电流 ~1.5 A/m² → 初始压降 ~1.5 mV，随 SEI 生长近似线性放大
    auto p_ref = chen2020_dfn();
    auto p_sei = chen2020_dfn();
    p_sei.degradation.enabled = true; // L0 = 5e-9 m、R_sei = 2e5 Ω·m

    DFNModel m_ref(p_ref, 5.0);
    DFNModel m_sei(p_sei, 5.0);
    m_ref.step(5.0, 10.0);
    m_sei.step(5.0, 10.0);
    const real_t dv1 = m_ref.voltage() - m_sei.voltage();
    EXPECT_GT(dv1, 5e-4);  // L0·R_sei·i_surf ≈ 1.5 mV
    EXPECT_LT(dv1, 5e-3);

    for (int i = 0; i < 29; ++i) {
        m_ref.step(5.0, 10.0);
        m_sei.step(5.0, 10.0);
    }
    const real_t dv30 = m_ref.voltage() - m_sei.voltage();
    // SEI 生长（300 s 放电 ~2.6 倍）使膜压降同步放大
    EXPECT_GT(dv30, 1.5 * dv1);
    EXPECT_LT(dv30, 0.030);
}

TEST(DFNModel, InvalidDegradationParametersRejected) {
    // 构造期校验（enabled 时生效）：正性/范围约束
    auto p = chen2020_dfn();
    p.degradation.enabled = true;
    p.degradation.l_sei_0 = -1e-9;
    EXPECT_THROW(DFNModel{p}, std::invalid_argument);

    p = chen2020_dfn();
    p.degradation.enabled = true;
    p.degradation.v_bar_sei = 0.0;
    EXPECT_THROW(DFNModel{p}, std::invalid_argument);

    p = chen2020_dfn();
    p.degradation.enabled = true;
    p.degradation.k_sei = 0.0;
    EXPECT_THROW(DFNModel{p}, std::invalid_argument);

    p = chen2020_dfn();
    p.degradation.enabled = true;
    p.degradation.d_ec = -1.0;
    EXPECT_THROW(DFNModel{p}, std::invalid_argument);

    p = chen2020_dfn();
    p.degradation.enabled = true;
    p.degradation.alpha_sei = 1.5;
    EXPECT_THROW(DFNModel{p}, std::invalid_argument);

    // disabled 时不校验（关闭 = 全路径零开销，诊断参数不抛错）
    p = chen2020_dfn();
    p.degradation.k_sei = -1.0;
    EXPECT_NO_THROW(DFNModel{p});
}

TEST(SimulationEngine, DfnCycleProtocolRepeatsStepsAndMarksCycleEnd) {
    // 循环工况：cycles = 2 → 3 s 段重复两圈（6 步 + t = 0 快照）；
    // cycle_number 跨圈递增，cycle_end 仅在每圈最后一段的截止快照为真
    SimulationEngine engine;
    MultiStepProtocol proto;
    proto.model = CellModelType::DFN;
    proto.time_step = 1.0;
    proto.cycles = 2;
    CCStep step;
    step.mode = CCStep::Mode::Discharge;
    step.current = 5.0;
    step.cutoff_time = 3.0;
    proto.steps.push_back(step);

    engine.single_step({}, proto); // 初始化 + t = 0 快照（圈 1）
    for (int i = 0; i < 6; ++i) {
        engine.single_step({}, proto);
    }
    EXPECT_EQ(engine.status(), SimulationStatus::Completed);

    const auto snaps = engine.collect_snapshots();
    ASSERT_EQ(snaps.size(), 7u); // t = 0..6
    for (std::size_t i = 0; i < snaps.size(); ++i) {
        EXPECT_EQ(snaps[i].cycle_number, i < 4 ? 1 : 2);
        EXPECT_EQ(snaps[i].cycle_end, i == 3 || i == 6);
    }
}

TEST(SimulationEngine, DfnDegradationFadesCapacityAcrossCycles) {
    // 循环 + 降解：满电静置圈（SEI 生长 = 自放电）下锂损失逐圈累积、
    // 可用容量单调衰减、SEI 厚度单调增长（寿命曲线数据源）
    SimulationEngine engine;
    MultiStepProtocol proto;
    proto.model = CellModelType::DFN;
    proto.time_step = 2.0;
    proto.cycles = 2;
    proto.dfn_params.degradation.enabled = true;
    CCStep rest;
    rest.mode = CCStep::Mode::Rest;
    rest.cutoff_time = 60.0; // 每圈静置 60 s（30 步）
    proto.steps.push_back(rest);

    engine.single_step({}, proto);
    for (int i = 0; i < 60; ++i) {
        engine.single_step({}, proto);
    }
    EXPECT_EQ(engine.status(), SimulationStatus::Completed);

    const auto snaps = engine.collect_snapshots();
    ASSERT_EQ(snaps.size(), 61u); // t = 0..120
    EXPECT_NEAR(snaps.front().li_loss_ah, 0.0, 1e-15);
    EXPECT_GT(snaps.back().li_loss_ah, 1e-4);

    // 每圈末（cycle_end）诊断：锂损失累积、容量衰减、SEI 增厚
    const SolutionSnapshot& end1 = snaps[30];
    ASSERT_TRUE(end1.cycle_end);
    EXPECT_EQ(end1.cycle_number, 1);
    const SolutionSnapshot& end2 = snaps[60];
    ASSERT_TRUE(end2.cycle_end);
    EXPECT_EQ(end2.cycle_number, 2);
    EXPECT_GT(end2.li_loss_ah, end1.li_loss_ah);
    EXPECT_NEAR(end2.degraded_capacity_ah, 5.0 - end2.li_loss_ah, 1e-12);
    EXPECT_GT(end2.sei_thickness, end1.sei_thickness);
    EXPECT_GT(end2.sei_thickness, 5e-9);
}

} // namespace
