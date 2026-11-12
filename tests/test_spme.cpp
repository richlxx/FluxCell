/// @file test_spme.cpp
/// @brief SPMe 模型单元测试：初始条件 / 电解液传输 / 液相过电位 / 守恒性 /
///        与 SPM 的退化一致性 / 引擎多模型接入 / 参数 IO
#include <fluxcell/core/model/spme_model.hpp>
#include <fluxcell/core/model/spm_model.hpp>
#include <fluxcell/core/simulation/simulation.hpp>
#include <fluxcell/core/io/io.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
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
using fluxcell::model::SPMModel;
using fluxcell::model::SPMeModel;

/// Chen2020 默认参数（SPMe 形式）
SPMeModel::Parameters chen2020_spme() {
    SPMeModel::Parameters p; // 默认值即 Chen2020
    return p;
}

/// 单段恒流放电 SPMe 协议
MultiStepProtocol spme_discharge(real_t current, real_t dt, real_t cutoff_time) {
    MultiStepProtocol proto;
    proto.model = CellModelType::SPMe;
    proto.time_step = dt;
    CCStep step;
    step.mode = CCStep::Mode::Discharge;
    step.current = current;
    step.cutoff_time = cutoff_time;
    proto.steps.push_back(step);
    return proto;
}

TEST(SPMeModel, InitialSocMapsStoichiometryLinearly) {
    auto params = chen2020_spme();
    params.solid.initial_soc = 0.5;
    SPMeModel model(params);

    const real_t x_expected = 0.5 * (params.solid.neg_stoich_at_soc0
                                   + params.solid.neg_stoich_at_soc1);
    const real_t y_expected = 0.5 * (params.solid.pos_stoich_at_soc0
                                   + params.solid.pos_stoich_at_soc1);
    EXPECT_NEAR(model.neg_stoichiometry(), x_expected, 1e-12);
    EXPECT_NEAR(model.pos_stoichiometry(), y_expected, 1e-12);
}

TEST(SPMeModel, RestKeepsElectrolyteUniform) {
    SPMeModel model(chen2020_spme(), 0.0);
    const real_t v0 = model.voltage();

    for (int i = 0; i < 5; ++i) {
        model.step(0.0, 60.0);
    }

    // 静置无源：c_e 保持均匀初值，phi_e 保持 0，电压不变
    const auto& c_e = model.electrolyte_conc();
    ASSERT_EQ(c_e.size(), 100);
    for (int i = 0; i < c_e.size(); ++i) {
        EXPECT_NEAR(c_e(i), 1000.0, 1e-9);
        EXPECT_NEAR(model.electrolyte_potential()(i), 0.0, 1e-12);
    }
    EXPECT_NEAR(model.electrolyte_overpotential(), 0.0, 1e-12);
    EXPECT_NEAR(model.voltage(), v0, 1e-12);
}

TEST(SPMeModel, SolidPhaseDegeneratesToSPM) {
    // SPMe 的固相方程与 SPM 相同（颗粒表面通量由总电流决定，与动力学解耦）：
    // 静置电压应严格一致；放电后固相嵌入率应一致（BV 的 i0 因液相浓度而异，
    // 只影响过电位，不影响固相浓度场）。
    SPMModel::Parameters sp; // 默认 Chen2020
    SPMModel spm(sp, 0.0);

    auto pp = chen2020_spme();
    pp.solid = sp;
    SPMeModel spme(pp, 0.0);

    EXPECT_NEAR(spm.voltage(), spme.voltage(), 1e-12);

    for (int i = 0; i < 10; ++i) {
        spm.step(5.0, 10.0);
        spme.step(5.0, 10.0);
    }
    EXPECT_NEAR(spm.neg_stoichiometry(), spme.neg_stoichiometry(), 1e-9);
    EXPECT_NEAR(spm.pos_stoichiometry(), spme.pos_stoichiometry(), 1e-9);
}

TEST(SPMeModel, DischargeBuildsElectrolyteGradient) {
    SPMeModel model(chen2020_spme(), 5.0);
    for (int i = 0; i < 10; ++i) {
        model.step(5.0, 10.0); // 1C 放电 100 s
    }

    // 放电：锂经电解液从负极流向正极，负极区浓度升高、正极区降低
    const auto& c_e = model.electrolyte_conc();
    EXPECT_GT(c_e(0), 1000.0);
    EXPECT_LT(c_e(c_e.size() - 1), 1000.0);
    EXPECT_GT(c_e(0), c_e(c_e.size() - 1));

    // 液相过电位 = 欧姆 + 浓差，放电方向为负（拉低端电压）
    EXPECT_LT(model.electrolyte_overpotential(), 0.0);
    EXPECT_LT(model.voltage(), model.pos_ocp() - model.neg_ocp()
                                  + model.pos_overpotential()
                                  - model.neg_overpotential());
}

TEST(SPMeModel, ElectrolyteOverpotentialMatchesOhmicAnalytic) {
    // 恒 κ 近似的解析欧姆积分（负极/正极线性电流分布 → L/(2κ) 系数）：
    // η_e = −j·(L_n/(2κ_n) + L_s/κ_s + L_p/(2κ_p))，均匀初浓度下无浓差项。
    // 离散求和与解析式存在 O(dx) 差（网格边界舍入 + 矩形积分），留 5% 容差。
    auto params = chen2020_spme();
    const real_t area = params.solid.electrode_height * params.solid.electrode_width;
    const real_t j = 5.0 / area; // 液相面电流密度 [A/m²]

    const real_t kappa_n = params.electrolyte_conductivity
                         * std::pow(params.porosity_neg, params.bruggeman_exponent);
    const real_t kappa_s = params.electrolyte_conductivity
                         * std::pow(params.porosity_sep, params.bruggeman_exponent);
    const real_t kappa_p = params.electrolyte_conductivity
                         * std::pow(params.porosity_pos, params.bruggeman_exponent);
    const real_t analytic =
        -j * (params.solid.neg_thickness / (2.0 * kappa_n)
              + params.solid.sep_thickness / kappa_s
              + params.solid.pos_thickness / (2.0 * kappa_p));

    SPMeModel model(params, 5.0); // 构造即计算 t = 0 过电位
    const real_t eta_e = model.electrolyte_overpotential();

    EXPECT_LT(eta_e, 0.0);
    EXPECT_NEAR(eta_e, analytic, 0.05 * std::abs(analytic));
}

TEST(SPMeModel, TotalLithiumConservedUnderDischarge) {
    SPMeModel model(chen2020_spme(), 5.0);
    const real_t total0 = model.total_lithium_mol();
    EXPECT_GT(total0, 0.1); // 5 Ah 电池固相锂约 0.19 mol

    for (int i = 0; i < 20; ++i) {
        model.step(5.0, 10.0);
    }
    EXPECT_NEAR(model.total_lithium_mol(), total0, 1e-9 * total0);
}

TEST(SPMeModel, SnapshotCarriesElectrolyteFields) {
    SPMeModel model(chen2020_spme(), 5.0);
    model.step(5.0, 10.0);

    const SolutionSnapshot snap = model.snapshot(5.0);
    EXPECT_EQ(snap.c_e.size(), 100);
    EXPECT_EQ(snap.phi_e.size(), 100);
    ASSERT_EQ(snap.c_s_surf_neg.size(), 1);
    ASSERT_EQ(snap.c_s_surf_pos.size(), 1);

    // 电压分解闭合：V = U_p − U_n + η_p − η_n + η_e
    EXPECT_NEAR(snap.u_pos - snap.u_neg + snap.eta_pos - snap.eta_neg + snap.eta_e,
                snap.voltage, 1e-12);
    EXPECT_NEAR(snap.eta_e, model.electrolyte_overpotential(), 1e-15);
}

TEST(SPMeModel, HigherTemperatureReducesActivationLoss) {
    auto cold = chen2020_spme();
    cold.solid.temperature = 298.15;
    auto hot = cold;
    hot.solid.temperature = 318.15;

    SPMeModel m_cold(cold, 5.0);
    SPMeModel m_hot(hot, 5.0);
    m_cold.step(5.0, 10.0);
    m_hot.step(5.0, 10.0);

    // Arrhenius：高温交换电流密度更大，活化过电位绝对值更小
    EXPECT_GT(m_hot.neg_overpotential(), 0.0);
    EXPECT_LT(m_hot.neg_overpotential(), m_cold.neg_overpotential());
}

TEST(SPMeModel, InvalidParametersRejected) {
    auto p = chen2020_spme();
    p.electrolyte_cells = 5;
    EXPECT_THROW(SPMeModel{p}, std::invalid_argument);

    p = chen2020_spme();
    p.porosity_neg = 0.0;
    EXPECT_THROW(SPMeModel{p}, std::invalid_argument);

    p = chen2020_spme();
    p.electrolyte_diffusivity = 0.0;
    EXPECT_THROW(SPMeModel{p}, std::invalid_argument);

    p = chen2020_spme();
    p.solid.temperature = 0.0;
    EXPECT_THROW(SPMeModel{p}, std::invalid_argument);

    p = chen2020_spme();
    p.solid.initial_soc = 1.5;
    EXPECT_THROW(SPMeModel{p}, std::invalid_argument);
}

TEST(SimulationEngine, SpmeSingleStepProducesElectrolyteData) {
    SimulationEngine engine;
    const MultiStepProtocol proto = spme_discharge(5.0, 1.0, 3.0);

    engine.single_step({}, proto); // 初始化 + t = 0 快照
    for (int i = 0; i < 3; ++i) {
        engine.single_step({}, proto);
    }
    EXPECT_EQ(engine.status(), SimulationStatus::Completed);

    const auto snaps = engine.collect_snapshots();
    ASSERT_EQ(snaps.size(), 4u); // t = 0..3
    EXPECT_EQ(snaps.front().c_e.size(), 100);
    EXPECT_EQ(snaps.front().phi_e.size(), 100);
    // 放电段每一步都有液相欧姆损失
    for (std::size_t i = 1; i < snaps.size(); ++i) {
        EXPECT_LT(snaps[i].eta_e, 0.0);
        EXPECT_NEAR(snaps[i].u_pos - snaps[i].u_neg + snaps[i].eta_pos
                        - snaps[i].eta_neg + snaps[i].eta_e,
                    snaps[i].voltage, 1e-12);
    }
}

TEST(Io, LoadSpmeParametersReadsElectrolyteSection) {
    const std::filesystem::path yml = std::filesystem::path(FLUXCELL_SOURCE_DIR)
                                    / "params" / "spme_chen2020.yml";
    const auto p = fluxcell::io::load_spme_parameters(yml);

    EXPECT_NEAR(p.electrolyte_diffusivity, 1.7694e-10, 1e-16);
    EXPECT_NEAR(p.transference_number, 0.2594, 1e-9);
    EXPECT_NEAR(p.porosity_neg, 0.25, 1e-9);
    EXPECT_NEAR(p.porosity_sep, 0.47, 1e-9);
    EXPECT_NEAR(p.porosity_pos, 0.335, 1e-9);
    EXPECT_NEAR(p.bruggeman_exponent, 1.5, 1e-9);
    EXPECT_NEAR(p.c_e_init, 1000.0, 1e-6);
    EXPECT_EQ(p.electrolyte_cells, 100);

    // 固相部分复用 SPM 读取
    EXPECT_NEAR(p.solid.capacity, 5.0, 1e-9);
    EXPECT_NEAR(p.solid.initial_soc, 1.0, 1e-9);
    EXPECT_NEAR(p.solid.neg_thickness, 85.2e-6, 1e-12);
}

TEST(Io, LoadSpmeParametersReadsTransportTable) {
    const std::filesystem::path yml = std::filesystem::path(FLUXCELL_SOURCE_DIR)
                                    / "params" / "spme_chen2020.yml";
    const auto p = fluxcell::io::load_spme_parameters(yml);

    // Nyman2008 标定表：76 节点，端点值与曲线一致
    ASSERT_EQ(p.electrolyte_conc_lut.size(), 76u);
    ASSERT_EQ(p.conductivity_lut.size(), 76u);
    ASSERT_EQ(p.diffusivity_lut.size(), 76u);
    EXPECT_NEAR(p.electrolyte_conc_lut.front(), 0.0, 1e-12);
    EXPECT_NEAR(p.electrolyte_conc_lut.back(), 5000.0, 1e-9);
    EXPECT_NEAR(p.conductivity_lut.front(), 0.0, 1e-12);
    EXPECT_NEAR(p.conductivity_lut.back(), 4.79485, 1e-6);
    EXPECT_NEAR(p.diffusivity_lut.front(), 4.862e-10, 1e-17);
    EXPECT_NEAR(p.diffusivity_lut.back(), 6.987e-10, 1e-16);

    // 有表时按浓度线性插值：c=1000 节点精确命中标定值
    EXPECT_NEAR(p.conductivity_at(1000.0), 0.9487, 1e-12);
    EXPECT_NEAR(p.diffusivity_at(1000.0), 1.7694e-10, 1e-16);
    // 节点值与端点外钳位
    EXPECT_NEAR(p.conductivity_at(50.0), 0.138404, 1e-7);
    EXPECT_NEAR(p.conductivity_at(-10.0), 0.0, 1e-12);
    EXPECT_NEAR(p.diffusivity_at(9000.0), 6.987e-10, 1e-16);
}

TEST(Io, TransportTableInconsistentLengthRejected) {
    // 长度不一致的表在模型构造时校验拒绝
    SPMeModel::Parameters p;
    p.electrolyte_conc_lut = {0.0, 1000.0, 5000.0};
    p.conductivity_lut = {0.0, 0.9487}; // 长度不匹配
    EXPECT_THROW(SPMeModel model(p), std::invalid_argument);
}

TEST(TransportTable, FallbackToConstantWhenEmpty) {
    // 无表时回退单点常量（既有行为不变）
    SPMeModel::Parameters p;
    EXPECT_NEAR(p.conductivity_at(1000.0), p.electrolyte_conductivity, 1e-15);
    EXPECT_NEAR(p.diffusivity_at(500.0), p.electrolyte_diffusivity, 1e-20);
}

// ══════════════════════════════════════════════════════════════
//  参数保存 round-trip（save → load 逐字段校验）
// ══════════════════════════════════════════════════════════════

namespace {

/// 固相参数逐字段位级相等（save_spm/spme/dfn_parameters round-trip 辅助）
void expect_solid_roundtrip(const SPMModel::Parameters& a,
                            const SPMModel::Parameters& b) {
    EXPECT_DOUBLE_EQ(a.capacity, b.capacity);
    EXPECT_DOUBLE_EQ(a.initial_soc, b.initial_soc);
    EXPECT_DOUBLE_EQ(a.temperature, b.temperature);
    EXPECT_EQ(a.particle_cells, b.particle_cells);
    EXPECT_DOUBLE_EQ(a.neg_thickness, b.neg_thickness);
    EXPECT_DOUBLE_EQ(a.neg_particle_radius, b.neg_particle_radius);
    EXPECT_DOUBLE_EQ(a.neg_active_fraction, b.neg_active_fraction);
    EXPECT_DOUBLE_EQ(a.neg_diffusion, b.neg_diffusion);
    EXPECT_DOUBLE_EQ(a.neg_max_conc, b.neg_max_conc);
    EXPECT_DOUBLE_EQ(a.neg_m_ref, b.neg_m_ref);
    EXPECT_DOUBLE_EQ(a.neg_activation_energy, b.neg_activation_energy);
    EXPECT_DOUBLE_EQ(a.pos_thickness, b.pos_thickness);
    EXPECT_DOUBLE_EQ(a.pos_particle_radius, b.pos_particle_radius);
    EXPECT_DOUBLE_EQ(a.pos_active_fraction, b.pos_active_fraction);
    EXPECT_DOUBLE_EQ(a.pos_diffusion, b.pos_diffusion);
    EXPECT_DOUBLE_EQ(a.pos_max_conc, b.pos_max_conc);
    EXPECT_DOUBLE_EQ(a.pos_m_ref, b.pos_m_ref);
    EXPECT_DOUBLE_EQ(a.pos_activation_energy, b.pos_activation_energy);
    EXPECT_DOUBLE_EQ(a.sep_thickness, b.sep_thickness);
    EXPECT_DOUBLE_EQ(a.electrolyte_conc, b.electrolyte_conc);
    EXPECT_DOUBLE_EQ(a.electrode_height, b.electrode_height);
    EXPECT_DOUBLE_EQ(a.electrode_width, b.electrode_width);
    EXPECT_DOUBLE_EQ(a.neg_stoich_at_soc0, b.neg_stoich_at_soc0);
    EXPECT_DOUBLE_EQ(a.neg_stoich_at_soc1, b.neg_stoich_at_soc1);
    EXPECT_DOUBLE_EQ(a.pos_stoich_at_soc0, b.pos_stoich_at_soc0);
    EXPECT_DOUBLE_EQ(a.pos_stoich_at_soc1, b.pos_stoich_at_soc1);
}

} // namespace

TEST(Io, SaveLoadRoundTripSpmParameters) {
    const std::filesystem::path src = std::filesystem::path(FLUXCELL_SOURCE_DIR)
                                    / "params" / "chen2020.yml";
    const auto p = fluxcell::io::load_spm_parameters(src);

    const auto tmp = std::filesystem::temp_directory_path()
                   / "fluxcell_test_spm_roundtrip.yml";
    fluxcell::io::save_spm_parameters(tmp, p);
    const auto q = fluxcell::io::load_spm_parameters(tmp);
    std::error_code ec;
    std::filesystem::remove(tmp, ec);

    expect_solid_roundtrip(p, q);
    // OCP 表（若源文件提供）逐点还原（相对容差兼容长尾小数）
    ASSERT_EQ(q.neg_ocp_stoich.size(), p.neg_ocp_stoich.size());
    ASSERT_EQ(q.pos_ocp_stoich.size(), p.pos_ocp_stoich.size());
    for (std::size_t i = 0; i < p.neg_ocp_stoich.size(); ++i) {
        EXPECT_NEAR(p.neg_ocp_stoich[i], q.neg_ocp_stoich[i],
                    1e-12 * std::abs(p.neg_ocp_stoich[i]));
        EXPECT_NEAR(p.neg_ocp_voltage[i], q.neg_ocp_voltage[i],
                    1e-12 * std::abs(p.neg_ocp_voltage[i]));
    }
    for (std::size_t i = 0; i < p.pos_ocp_stoich.size(); ++i) {
        EXPECT_NEAR(p.pos_ocp_stoich[i], q.pos_ocp_stoich[i],
                    1e-12 * std::abs(p.pos_ocp_stoich[i]));
        EXPECT_NEAR(p.pos_ocp_voltage[i], q.pos_ocp_voltage[i],
                    1e-12 * std::abs(p.pos_ocp_voltage[i]));
    }
}

TEST(Io, SaveLoadRoundTripSpmeParameters) {
    const std::filesystem::path src = std::filesystem::path(FLUXCELL_SOURCE_DIR)
                                    / "params" / "spme_chen2020.yml";
    const auto p = fluxcell::io::load_spme_parameters(src);

    const auto tmp = std::filesystem::temp_directory_path()
                   / "fluxcell_test_spme_roundtrip.yml";
    fluxcell::io::save_spme_parameters(tmp, p);
    const auto q = fluxcell::io::load_spme_parameters(tmp);
    std::error_code ec;
    std::filesystem::remove(tmp, ec);

    expect_solid_roundtrip(p.solid, q.solid);
    EXPECT_DOUBLE_EQ(p.electrolyte_diffusivity, q.electrolyte_diffusivity);
    EXPECT_DOUBLE_EQ(p.transference_number, q.transference_number);
    EXPECT_DOUBLE_EQ(p.porosity_neg, q.porosity_neg);
    EXPECT_DOUBLE_EQ(p.porosity_sep, q.porosity_sep);
    EXPECT_DOUBLE_EQ(p.porosity_pos, q.porosity_pos);
    EXPECT_DOUBLE_EQ(p.bruggeman_exponent, q.bruggeman_exponent);
    EXPECT_DOUBLE_EQ(p.electrolyte_conductivity, q.electrolyte_conductivity);
    EXPECT_DOUBLE_EQ(p.c_e_init, q.c_e_init);
    EXPECT_EQ(p.electrolyte_cells, q.electrolyte_cells);
    // 浓度相关输运表（76 点 Nyman2008）逐点还原
    ASSERT_EQ(q.electrolyte_conc_lut.size(), p.electrolyte_conc_lut.size());
    ASSERT_EQ(q.conductivity_lut.size(), p.conductivity_lut.size());
    ASSERT_EQ(q.diffusivity_lut.size(), p.diffusivity_lut.size());
    for (std::size_t i = 0; i < p.electrolyte_conc_lut.size(); ++i) {
        EXPECT_NEAR(p.electrolyte_conc_lut[i], q.electrolyte_conc_lut[i],
                    1e-12 * std::abs(p.electrolyte_conc_lut[i]));
        EXPECT_NEAR(p.conductivity_lut[i], q.conductivity_lut[i],
                    1e-12 * std::abs(p.conductivity_lut[i]));
        EXPECT_NEAR(p.diffusivity_lut[i], q.diffusivity_lut[i],
                    1e-12 * std::abs(p.diffusivity_lut[i]));
    }
}

TEST(Io, SaveLoadRoundTripDfnParameters) {
    const std::filesystem::path src = std::filesystem::path(FLUXCELL_SOURCE_DIR)
                                    / "params" / "dfn_chen2020.yml";
    const auto p = fluxcell::io::load_dfn_parameters(src);

    const auto tmp = std::filesystem::temp_directory_path()
                   / "fluxcell_test_dfn_roundtrip.yml";
    fluxcell::io::save_dfn_parameters(tmp, p);
    const auto q = fluxcell::io::load_dfn_parameters(tmp);
    std::error_code ec;
    std::filesystem::remove(tmp, ec);

    expect_solid_roundtrip(p.solid, q.solid);
    EXPECT_DOUBLE_EQ(p.neg_conductivity, q.neg_conductivity);
    EXPECT_DOUBLE_EQ(p.pos_conductivity, q.pos_conductivity);
    // 热模型段 round-trip（DFN 专属，与读取/写出段对称）
    EXPECT_EQ(p.thermal_type, q.thermal_type);
    EXPECT_DOUBLE_EQ(p.thermal_params.t_ambient, q.thermal_params.t_ambient);
    EXPECT_DOUBLE_EQ(p.thermal_params.t_initial, q.thermal_params.t_initial);
    EXPECT_DOUBLE_EQ(p.thermal_params.h_total, q.thermal_params.h_total);
    EXPECT_DOUBLE_EQ(p.thermal_params.cooling_area, q.thermal_params.cooling_area);
    EXPECT_DOUBLE_EQ(p.thermal_params.cell_volume, q.thermal_params.cell_volume);
    EXPECT_DOUBLE_EQ(p.thermal_params.rho_c_p_neg, q.thermal_params.rho_c_p_neg);
    EXPECT_DOUBLE_EQ(p.thermal_params.rho_c_p_sep, q.thermal_params.rho_c_p_sep);
    EXPECT_DOUBLE_EQ(p.thermal_params.rho_c_p_pos, q.thermal_params.rho_c_p_pos);
    EXPECT_DOUBLE_EQ(p.thermal_params.rho_c_p_cc_neg,
                     q.thermal_params.rho_c_p_cc_neg);
    EXPECT_DOUBLE_EQ(p.thermal_params.rho_c_p_cc_pos,
                     q.thermal_params.rho_c_p_cc_pos);
    EXPECT_DOUBLE_EQ(p.thermal_params.cc_thickness_neg,
                     q.thermal_params.cc_thickness_neg);
    EXPECT_DOUBLE_EQ(p.thermal_params.cc_thickness_pos,
                     q.thermal_params.cc_thickness_pos);
    EXPECT_DOUBLE_EQ(p.thermal_params.entropic_change_neg,
                     q.thermal_params.entropic_change_neg);
    EXPECT_DOUBLE_EQ(p.thermal_params.entropic_change_pos,
                     q.thermal_params.entropic_change_pos);
    // 熵变 LUT：源文件用 CSV 路径引用加载，内联映射写出后读回应逐点一致
    ASSERT_EQ(p.thermal_params.entropic_lut_neg.size(),
              q.thermal_params.entropic_lut_neg.size());
    ASSERT_GT(p.thermal_params.entropic_lut_neg.size(), std::size_t{2});
    ASSERT_EQ(p.thermal_params.entropic_lut_pos.size(),
              q.thermal_params.entropic_lut_pos.size());
    ASSERT_GT(p.thermal_params.entropic_lut_pos.size(), std::size_t{2});
    for (std::size_t i = 0; i < p.thermal_params.entropic_lut_neg.size(); ++i) {
        EXPECT_NEAR(p.thermal_params.entropic_lut_neg.sto[i],
                    q.thermal_params.entropic_lut_neg.sto[i],
                    1e-12 * std::abs(p.thermal_params.entropic_lut_neg.sto[i]));
        EXPECT_NEAR(p.thermal_params.entropic_lut_neg.dudt[i],
                    q.thermal_params.entropic_lut_neg.dudt[i],
                    1e-12 * std::abs(p.thermal_params.entropic_lut_neg.dudt[i])
                        + 1e-15);
    }
    for (std::size_t i = 0; i < p.thermal_params.entropic_lut_pos.size(); ++i) {
        EXPECT_NEAR(p.thermal_params.entropic_lut_pos.sto[i],
                    q.thermal_params.entropic_lut_pos.sto[i],
                    1e-12 * std::abs(p.thermal_params.entropic_lut_pos.sto[i]));
        EXPECT_NEAR(p.thermal_params.entropic_lut_pos.dudt[i],
                    q.thermal_params.entropic_lut_pos.dudt[i],
                    1e-12 * std::abs(p.thermal_params.entropic_lut_pos.dudt[i])
                        + 1e-15);
    }
    // 电解液扩展结构与输运表在 SPMe round-trip 已逐点覆盖，此处抽查
    EXPECT_DOUBLE_EQ(p.electrolyte_diffusivity, q.electrolyte_diffusivity);
    ASSERT_EQ(q.electrolyte_conc_lut.size(), p.electrolyte_conc_lut.size());
    EXPECT_NEAR(q.conductivity_at(1000.0), 0.9487, 1e-12);
    EXPECT_NEAR(q.diffusivity_at(1000.0), 1.7694e-10, 1e-16);
}

} // namespace
