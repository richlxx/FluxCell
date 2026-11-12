/// @file test_spm.cpp
/// @brief 单颗粒模型 (SPM) 单元测试
#include <fluxcell/core/model/spm_model.hpp>
#include <fluxcell/core/io/io.hpp>
#include <fluxcell/core/physics/physics.hpp>
#include <fluxcell/core/simulation/simulation.hpp>

#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iterator>
#include <stdexcept>
#include <thread>
#include <vector>

#include <algorithm>

using namespace fluxcell;
using fluxcell::model::SPMModel;

// ══════════════════════════════════════════════════════════════
//  Chen2020 本构关系
// ══════════════════════════════════════════════════════════════

TEST(BVInverse, RoundTripWithButlerVolmer) {
    // η = (2RT/F)·asinh(i/2i0) 反解后应使对称 BV 精确还原 i
    const real_t T = 298.15;
    const real_t i0 = 3.5;
    for (real_t i : {0.1, 1.0, 5.0, 20.0, -8.0, -100.0}) {
        const real_t eta = physics::overpotential_from_bv(i, i0, T);
        // 对称 BV：j = i0·[exp(Fη/2RT) - exp(-Fη/2RT)] = 2·i0·sinh(Fη/2RT)
        const real_t j = i0 * (std::exp(0.5 * 96485.0 / (8.314 * T) * eta)
                             - std::exp(-0.5 * 96485.0 / (8.314 * T) * eta));
        EXPECT_NEAR(j, i, std::abs(i) * 1e-9) << "at i=" << i;
    }
}

TEST(BVInverse, SignConvention) {
    const real_t T = 298.15;
    EXPECT_GT(physics::overpotential_from_bv(1.0, 1.0, T), 0.0);
    EXPECT_LT(physics::overpotential_from_bv(-1.0, 1.0, T), 0.0);
    EXPECT_NEAR(physics::overpotential_from_bv(0.0, 1.0, T), 0.0, 1e-15);
}

TEST(Chen2020OCP, GraphiteRangeAndTrend) {
    // 石墨 OCP：x 增大 → 电位降低，工作区间 [0.02, 0.95] 内 0~1.5V
    for (real_t x = 0.02; x <= 0.95; x += 0.05) {
        const real_t u = physics::graphite_chen2020_ocp(x);
        EXPECT_GE(u, -0.05) << "at x=" << x;
        EXPECT_LE(u, 1.5) << "at x=" << x;
    }
    EXPECT_GT(physics::graphite_chen2020_ocp(0.1), physics::graphite_chen2020_ocp(0.9));
}

TEST(Chen2020OCP, NMCAgreesWithNMC111Fit) {
    // Chen2020 NMC OCP 与既有 NMC111 拟合公式一致（PyBaMM 中同源）
    for (real_t y = 0.3; y <= 0.95; y += 0.05) {
        EXPECT_NEAR(physics::nmc_chen2020_ocp(y), physics::nmc111_ocp(y), 1e-12)
            << "at y=" << y;
    }
}

TEST(SPMOcpLut, YamlLutMatchesBuiltinFormula) {
    // params/chen2020.yml 的 ocp: 段 = 内置 Chen2020 公式的采样点；
    // 恒流放电全程对比两条路径的端电压（201 点 LUT 线性插值误差 < 1mV）
    const auto lut_params = io::load_spm_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params" / "chen2020.yml");
    ASSERT_FALSE(lut_params.neg_ocp_stoich.empty());
    ASSERT_FALSE(lut_params.pos_ocp_stoich.empty());

    SPMModel m_lut(lut_params);
    SPMModel m_formula(SPMModel::Parameters{}); // 默认值 = 内置 Chen2020 公式

    for (int k = 0; k < 300; ++k) {
        m_lut.step(2.0, 1.0);
        m_formula.step(2.0, 1.0);
        ASSERT_NEAR(m_lut.voltage(), m_formula.voltage(), 2e-3)
            << "t=" << m_lut.time();
    }
}

TEST(Chen2020ExchangeCurrent, PositiveAndTemperatureDependent) {
    const real_t i0_ref = physics::exchange_current_density_chen2020_neg(
        1000.0, 15000.0, 33133.0, 298.15);
    EXPECT_GT(i0_ref, 0.0);

    // Arrhenius：升温 → i0 增大（正活化能）
    const real_t i0_hot = physics::exchange_current_density_chen2020_neg(
        1000.0, 15000.0, 33133.0, 318.15);
    EXPECT_GT(i0_hot, i0_ref);
}

// ══════════════════════════════════════════════════════════════
//  SPMModel
// ══════════════════════════════════════════════════════════════

namespace {
SPMModel::Parameters chen2020_params() {
    SPMModel::Parameters p; // 默认值即 Chen2020
    return p;
}

/// 单段恒流放电协议（引擎测试的便捷构造）
MultiStepProtocol cc_protocol(real_t current, real_t time_step, real_t cutoff_time,
                              real_t cutoff_voltage = 0.0) {
    MultiStepProtocol proto;
    proto.time_step = time_step;
    CCStep step;
    step.mode = CCStep::Mode::Discharge;
    step.current = current;
    step.cutoff_time = cutoff_time;
    step.cutoff_voltage = cutoff_voltage;
    proto.steps.push_back(step);
    return proto;
}
} // namespace // namespace

TEST(SPMModel, InitialVoltageMatchesOCV) {
    SPMModel model(chen2020_params());
    // 100% SOC：V0 = U_p(0.2700) - U_n(0.9014) ≈ 4.19V（PyBaMM SPM 同参数）
    EXPECT_NEAR(model.voltage(), 4.19, 0.05);
    EXPECT_NEAR(model.time(), 0.0, 1e-12);
}

TEST(SPMModel, ZeroCurrentKeepsState) {
    SPMModel model(chen2020_params());
    const real_t v0 = model.voltage();
    for (int i = 0; i < 100; ++i) {
        model.step(0.0, 1.0);
    }
    EXPECT_NEAR(model.voltage(), v0, 1e-9);
    EXPECT_NEAR(model.time(), 100.0, 1e-9);
}

TEST(SPMModel, LithiumConservationUnderCC) {
    SPMModel model(chen2020_params());
    const real_t m0 = model.total_lithium_mol();

    // 1C 放电 600s，远离截止，不触发限幅
    for (int i = 0; i < 600; ++i) {
        model.step(5.0, 1.0);
    }
    const real_t m1 = model.total_lithium_mol();
    EXPECT_NEAR(m1, m0, std::abs(m0) * 1e-6);
}

TEST(SPMModel, DischargeDepletesNegChargesPos) {
    SPMModel model(chen2020_params());
    const real_t x0 = model.neg_stoichiometry();
    const real_t y0 = model.pos_stoichiometry();

    for (int i = 0; i < 300; ++i) {
        model.step(5.0, 1.0);
    }
    EXPECT_LT(model.neg_stoichiometry(), x0); // 负极脱锂
    EXPECT_GT(model.pos_stoichiometry(), y0); // 正极嵌锂
}

TEST(SPMModel, DischargeVoltageDecreases) {
    SPMModel model(chen2020_params());
    real_t v_prev = model.voltage();
    for (int i = 0; i < 900; ++i) {
        model.step(5.0, 1.0);
        EXPECT_LE(model.voltage(), v_prev + 1e-6) << "step " << i;
        v_prev = model.voltage();
    }
}

TEST(SPMModel, PyBaMMReferencePoint1C1800s) {
    // PyBaMM SPM（Chen2020, I=5A, dt=1s）t=1800s: V = 3.5684V
    SPMModel model(chen2020_params());
    for (int i = 0; i < 1800; ++i) {
        model.step(5.0, 1.0);
    }
    EXPECT_NEAR(model.voltage(), 3.5684, 0.01);
}

TEST(SPMModel, TimeStepConvergence) {
    // dt = 1s 与 dt = 0.5s 结果应接近（二阶 CN）
    SPMModel m1(chen2020_params());
    SPMModel m05(chen2020_params());
    for (int i = 0; i < 1800; ++i) m1.step(5.0, 1.0);
    for (int i = 0; i < 3600; ++i) m05.step(5.0, 0.5);
    EXPECT_NEAR(m1.voltage(), m05.voltage(), 0.005);
}

TEST(SPMModel, CRateScales) {
    // 2C 放电同一时刻电压应低于 1C（更极端状态）
    SPMModel m1c(chen2020_params());
    SPMModel m2c(chen2020_params());
    for (int i = 0; i < 600; ++i) {
        m1c.step(5.0, 1.0);
        m2c.step(10.0, 1.0);
    }
    EXPECT_LT(m2c.voltage(), m1c.voltage());
}

TEST(SPMModel, SnapshotVoltageDecompositionConsistent) {
    // 快照分量满足 V = U_p - U_n + η_p - η_n，且化学计量数在 (0,1) 内
    SPMModel model(chen2020_params(), 5.0);
    for (int i = 0; i < 120; ++i) model.step(5.0, 1.0);

    const auto snap = model.snapshot(5.0);
    EXPECT_EQ(snap.c_s_neg.size(), chen2020_params().particle_cells);
    EXPECT_EQ(snap.c_s_pos.size(), chen2020_params().particle_cells);
    EXPECT_NEAR(snap.u_pos - snap.u_neg + snap.eta_pos - snap.eta_neg,
                snap.voltage, 1e-12);
    EXPECT_GT(snap.x_surf, 0.0);
    EXPECT_LT(snap.x_surf, 1.0);
    EXPECT_GT(snap.y_surf, 0.0);
    EXPECT_LT(snap.y_surf, 1.0);
    // 放电：石墨脱锂（x 下降）、NMC 嵌锂（y 上升）；
    // i_n > 0（锂流出负极）→ η_n > 0，i_p < 0 → η_p < 0，故 η_p - η_n < 0（V 低于开路）
    const auto& p = chen2020_params();
    EXPECT_LT(snap.x_surf, p.neg_stoich_at_soc1); // 初始 SOC = 100%
    EXPECT_GT(snap.y_surf, p.pos_stoich_at_soc1);
    EXPECT_GT(snap.eta_neg, 0.0);
    EXPECT_LT(snap.eta_pos, 0.0);
}

// ══════════════════════════════════════════════════════════
//  初始 SOC 与环境温度
// ══════════════════════════════════════════════════════════

TEST(SPMModel, InitialSocMapsStoichiometryLinearly) {
    auto params = chen2020_params();
    params.initial_soc = 0.5;
    SPMModel model(params);

    // 初始均匀浓度下表面外推 = 本体：SOC 线性映射端点中点
    const real_t x_expected = 0.5 * (params.neg_stoich_at_soc0 + params.neg_stoich_at_soc1);
    const real_t y_expected = 0.5 * (params.pos_stoich_at_soc0 + params.pos_stoich_at_soc1);
    EXPECT_NEAR(model.neg_stoichiometry(), x_expected, 1e-12);
    EXPECT_NEAR(model.pos_stoichiometry(), y_expected, 1e-12);

    // OCV 随 SOC 单调上升：V(100%) > V(50%) > V(0%)
    auto soc0_params = chen2020_params();
    soc0_params.initial_soc = 0.0;
    SPMModel soc1(chen2020_params());
    SPMModel soc0(soc0_params);
    EXPECT_GT(soc1.voltage(), model.voltage());
    EXPECT_GT(model.voltage(), soc0.voltage());
}

TEST(SPMModel, AmbientTemperatureReducesOverpotential) {
    auto cold_params = chen2020_params();
    cold_params.temperature = 298.15;
    auto hot_params = cold_params;
    hot_params.temperature = 318.15;

    SPMModel cold(cold_params);
    SPMModel hot(hot_params);
    cold.step(5.0, 1.0);
    hot.step(5.0, 1.0);

    // 升温 → 交换电流密度增大（Arrhenius）→ 活化过电位绝对值减小
    EXPECT_LT(std::abs(hot.neg_overpotential()), std::abs(cold.neg_overpotential()));
    EXPECT_LT(std::abs(hot.pos_overpotential()), std::abs(cold.pos_overpotential()));
    // 放电时端电压更高
    EXPECT_GT(hot.voltage(), cold.voltage());
}

TEST(SPMModel, InvalidTemperatureAndSocRejected) {
    auto params = chen2020_params();

    params.temperature = 0.0;
    EXPECT_THROW(SPMModel{params}, std::invalid_argument);
    params.temperature = -5.0;
    EXPECT_THROW(SPMModel{params}, std::invalid_argument);

    params.temperature = 298.15;
    params.initial_soc = -0.1;
    EXPECT_THROW(SPMModel{params}, std::invalid_argument);
    params.initial_soc = 1.5;
    EXPECT_THROW(SPMModel{params}, std::invalid_argument);
}

// ══════════════════════════════════════════════════════════════
//  SimulationEngine
// ══════════════════════════════════════════════════════════════

TEST(SimulationEngine, CompletesShortRun) {
    SimulationEngine engine;
    MultiStepProtocol proto = cc_protocol(5.0, 1.0, 10.0);
    proto.real_time_factor = 0.0; // 全速

    engine.start(chen2020_params(), proto);

    // 等待完成（轮询上限 5s）
    int waited = 0;
    while (engine.status() != SimulationStatus::Completed && waited < 500) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        ++waited;
    }
    ASSERT_EQ(engine.status(), SimulationStatus::Completed);

    auto snaps = engine.collect_snapshots();
    EXPECT_EQ(snaps.size(), 11u); // t = 0..10
    EXPECT_NEAR(snaps.front().time, 0.0, 1e-9);
    EXPECT_NEAR(snaps.back().time, 10.0, 1e-9);
    // 放电电压单调下降
    for (size_t i = 1; i < snaps.size(); ++i) {
        EXPECT_LE(snaps[i].voltage, snaps[i - 1].voltage + 1e-6);
    }
}

TEST(SimulationEngine, StopWorksMidRun) {
    SimulationEngine engine;
    MultiStepProtocol proto = cc_protocol(5.0, 1.0, 3600.0);
    proto.real_time_factor = 0.0;

    engine.start(chen2020_params(), proto);
    std::this_thread::sleep_for(std::chrono::milliseconds(50));
    engine.stop();

    EXPECT_EQ(engine.status(), SimulationStatus::Idle);
    // 引擎可重启
    engine.start(chen2020_params(), proto);
    engine.stop();
    EXPECT_EQ(engine.status(), SimulationStatus::Idle);
}

TEST(SimulationEngine, PauseThenResumeKeepsRunning) {
    // 回归：连续运行中暂停后必须可恢复（pause() 需同步 status_，
    // 否则 resume() 因状态非 Paused 而永不唤醒线程）
    SimulationEngine engine;
    MultiStepProtocol proto = cc_protocol(5.0, 1.0, 3600.0);
    proto.real_time_factor = 5.0; // 限速：每仿真秒 ≈ 200ms 墙钟

    engine.start(chen2020_params(), proto);
    std::this_thread::sleep_for(std::chrono::milliseconds(300));
    ASSERT_EQ(engine.status(), SimulationStatus::Running);
    engine.collect_snapshots();

    // 暂停：状态立即变为 Paused
    engine.pause();
    EXPECT_EQ(engine.status(), SimulationStatus::Paused);

    // 等 in-flight 步收尾后清空队列，静默期内不应再有新输出
    std::this_thread::sleep_for(std::chrono::milliseconds(300));
    engine.collect_snapshots();
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    EXPECT_EQ(engine.collect_snapshots().size(), 0u);

    // 恢复：状态回到 Running，曲线继续推进
    engine.resume();
    EXPECT_EQ(engine.status(), SimulationStatus::Running);
    std::this_thread::sleep_for(std::chrono::milliseconds(400));
    EXPECT_GT(engine.collect_snapshots().size(), 0u);

    engine.stop();
    EXPECT_EQ(engine.status(), SimulationStatus::Idle);
}

// ══════════════════════════════════════════════════════════
//  SimulationEngine 恒压 / 恒功率工况
// ══════════════════════════════════════════════════════════

TEST(SimulationEngine, CVStepTracksVoltageAndCutsOffByCurrent) {
    // CC 充电 → CV：端电压跟踪 4.2V，尾电流衰减到阈值时段结束。
    // 初始 SOC 0.5：默认满电起步充电会立刻过充发散
    SimulationEngine engine;
    MultiStepProtocol proto = cc_protocol(5.0, 1.0, 900.0);
    proto.steps.front().mode = CCStep::Mode::Charge;

    CCStep cv;
    cv.mode = CCStep::Mode::Charge;
    cv.control = CCStep::Control::Voltage;
    cv.voltage = 4.2;
    cv.cutoff_time = 0.0;
    cv.cutoff_current = 4.0;
    proto.steps.push_back(cv);

    auto params = chen2020_params();
    params.initial_soc = 0.5;
    engine.start(params, proto);
    int waited = 0;
    while (engine.status() != SimulationStatus::Completed && waited < 5000) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        ++waited;
    }
    ASSERT_EQ(engine.status(), SimulationStatus::Completed);

    auto snaps = engine.collect_snapshots();
    ASSERT_GE(snaps.size(), 3u);
    // CV 段（t >= 900）：端电压跟踪目标（容差 0.1V 吸收步级波动）
    for (const auto& s : snaps) {
        if (s.time >= 900.0) {
            EXPECT_NEAR(s.voltage, 4.2, 0.1) << "t = " << s.time;
        }
    }
    // 尾电流截止：末步 |I| <= 阈值（10% 容差）
    EXPECT_LE(std::abs(snaps.back().current), 4.4);
}

TEST(SimulationEngine, CPStepMaintainsPower) {
    // 恒功率放电：每步末 U·I ≈ P（步初电压换算 → 步末偏差为步级误差）
    SimulationEngine engine;
    MultiStepProtocol proto;
    proto.time_step = 2.0;
    CCStep cp;
    cp.mode = CCStep::Mode::Discharge;
    cp.control = CCStep::Control::Power;
    cp.power = 20.0;
    cp.cutoff_time = 60.0;
    proto.steps.push_back(cp);

    engine.start(chen2020_params(), proto);
    int waited = 0;
    while (engine.status() != SimulationStatus::Completed && waited < 1000) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        ++waited;
    }
    ASSERT_EQ(engine.status(), SimulationStatus::Completed);

    auto snaps = engine.collect_snapshots();
    ASSERT_GE(snaps.size(), 2u);
    for (std::size_t i = 1; i < snaps.size(); ++i) {
        const real_t p = snaps[i].voltage * std::abs(snaps[i].current);
        EXPECT_NEAR(p, 20.0, 2.0) << "t = " << snaps[i].time;
    }
}

TEST(SimulationEngine, MixedControlProtocolCompletes) {
    // 放电 CC → 静置 → CP 充电：状态连续传递，CP 充电电流为负。
    // 初始 SOC 0.5：默认满电起步充电会立刻过充
    SimulationEngine engine;
    MultiStepProtocol proto;
    proto.time_step = 2.0;

    CCStep cc = cc_protocol(5.0, 2.0, 60.0).steps.front();
    CCStep rest;
    rest.mode = CCStep::Mode::Rest;
    rest.cutoff_time = 60.0;
    CCStep cp;
    cp.mode = CCStep::Mode::Charge;
    cp.control = CCStep::Control::Power;
    cp.power = 15.0;
    cp.cutoff_time = 60.0;
    proto.steps = {cc, rest, cp};

    auto params = chen2020_params();
    params.initial_soc = 0.5;
    engine.start(params, proto);
    int waited = 0;
    while (engine.status() != SimulationStatus::Completed && waited < 1000) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        ++waited;
    }
    ASSERT_EQ(engine.status(), SimulationStatus::Completed);

    auto snaps = engine.collect_snapshots();
    ASSERT_GE(snaps.size(), 4u);
    // CP 充电段（t >= 120）：电流为负，且 U·I ≈ -P
    for (const auto& s : snaps) {
        if (s.time > 120.0) {
            EXPECT_LT(s.current, 0.0) << "t = " << s.time;
            EXPECT_NEAR(-s.voltage * s.current, 15.0, 1.5) << "t = " << s.time;
        }
    }
}

TEST(MultiStepProtocol, ValidateRejectsBadControlSteps) {
    MultiStepProtocol proto;
    proto.time_step = 1.0;

    // CV 目标电压必须为正
    CCStep cv;
    cv.control = CCStep::Control::Voltage;
    cv.voltage = 0.0;
    cv.cutoff_time = 600.0;
    proto.steps = {cv};
    EXPECT_THROW(proto.validate(), std::invalid_argument);

    // CV 缺少时间/尾电流截止
    cv.voltage = 4.2;
    cv.cutoff_time = 0.0;
    proto.steps = {cv};
    EXPECT_THROW(proto.validate(), std::invalid_argument);

    // 尾电流截止合法
    cv.cutoff_current = 0.5;
    proto.steps = {cv};
    EXPECT_NO_THROW(proto.validate());

    // 静置段不支持恒压控制
    CCStep rest_cv;
    rest_cv.mode = CCStep::Mode::Rest;
    rest_cv.control = CCStep::Control::Voltage;
    rest_cv.voltage = 4.2;
    rest_cv.cutoff_time = 600.0;
    proto.steps = {rest_cv};
    EXPECT_THROW(proto.validate(), std::invalid_argument);

    // CP 功率必须为正
    CCStep cp_bad;
    cp_bad.control = CCStep::Control::Power;
    cp_bad.power = 0.0;
    cp_bad.cutoff_time = 600.0;
    proto.steps = {cp_bad};
    EXPECT_THROW(proto.validate(), std::invalid_argument);
}

// ══════════════════════════════════════════════════════════
//  SimulationEngine 单步模式
// ══════════════════════════════════════════════════════════

TEST(SimulationEngine, SingleStepAdvancesAndPauses) {
    SimulationEngine engine;
    MultiStepProtocol proto = cc_protocol(5.0, 1.0, 3600.0);

    engine.single_step(chen2020_params(), proto);
    engine.single_step(chen2020_params(), proto);
    engine.single_step(chen2020_params(), proto);

    EXPECT_EQ(engine.status(), SimulationStatus::Paused);

    auto snaps = engine.collect_snapshots();
    ASSERT_EQ(snaps.size(), 4u); // t = 0..3
    EXPECT_NEAR(snaps.front().time, 0.0, 1e-9);
    EXPECT_NEAR(snaps.back().time, 3.0, 1e-9);
    // 单步电压与直接模型步进一致
    SPMModel ref(chen2020_params(), 5.0);
    for (int i = 0; i < 3; ++i) ref.step(5.0, 1.0);
    EXPECT_NEAR(snaps.back().voltage, ref.voltage(), 1e-12);

    engine.stop();
}

TEST(SimulationEngine, SingleStepThenResumeCompletes) {
    SimulationEngine engine;
    MultiStepProtocol proto = cc_protocol(5.0, 1.0, 10.0);

    engine.single_step(chen2020_params(), proto);
    engine.single_step(chen2020_params(), proto);
    ASSERT_EQ(engine.status(), SimulationStatus::Paused);

    // 恢复连续运行直至完成
    engine.resume();
    int waited = 0;
    while (engine.status() != SimulationStatus::Completed && waited < 500) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        ++waited;
    }
    ASSERT_EQ(engine.status(), SimulationStatus::Completed);

    auto snaps = engine.collect_snapshots();
    EXPECT_EQ(snaps.size(), 11u); // t = 0..10，单步的两帧连续衔接
    EXPECT_NEAR(snaps.front().time, 0.0, 1e-9);
    EXPECT_NEAR(snaps.back().time, 10.0, 1e-9);
}

TEST(SimulationEngine, SingleStepAfterStopStartsFresh) {
    SimulationEngine engine;
    MultiStepProtocol proto = cc_protocol(5.0, 1.0, 3600.0);

    engine.single_step(chen2020_params(), proto);
    engine.single_step(chen2020_params(), proto);
    engine.stop();
    ASSERT_EQ(engine.status(), SimulationStatus::Idle);

    // stop 后单步应从 t = 0 重新开始（旧进度作废）
    engine.single_step(chen2020_params(), proto);
    EXPECT_EQ(engine.status(), SimulationStatus::Paused);
    auto snaps = engine.collect_snapshots();
    ASSERT_EQ(snaps.size(), 2u);
    EXPECT_NEAR(snaps.front().time, 0.0, 1e-9);
    EXPECT_NEAR(snaps.back().time, 1.0, 1e-9);
    engine.stop();
}

TEST(SimulationEngine, SingleStepDuringRunKeepsProgress) {
    SimulationEngine engine;
    MultiStepProtocol proto = cc_protocol(5.0, 1.0, 3600.0);
    proto.real_time_factor = 10.0; // 限速回放：保证 50ms 窗口内仍处于运行中

    engine.start(chen2020_params(), proto);
    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    // 运行中单步：停线程续走当前进度，不重置
    engine.single_step(chen2020_params(), proto);
    EXPECT_EQ(engine.status(), SimulationStatus::Paused);

    auto snaps = engine.collect_snapshots();
    EXPECT_GE(snaps.size(), 2u);
    // 时间单调且无重叠（续接而非重启）
    for (size_t i = 1; i < snaps.size(); ++i) {
        EXPECT_GT(snaps[i].time, snaps[i - 1].time);
    }
    engine.stop();
}

// ══════════════════════════════════════════════════════════
//  SimulationEngine 多段工况
// ══════════════════════════════════════════════════════════

TEST(SimulationEngine, SocCoulombCounting) {
    SimulationEngine engine;
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    proto.real_time_factor = 0.0;
    // 放电 600 s → 充电 600 s（1C@5Ah 对称充放）
    CCStep discharge;
    discharge.mode = CCStep::Mode::Discharge;
    discharge.current = 5.0;
    discharge.cutoff_time = 600.0;
    CCStep charge;
    charge.mode = CCStep::Mode::Charge;
    charge.current = 5.0;
    charge.cutoff_time = 600.0;
    proto.steps.push_back(discharge);
    proto.steps.push_back(charge);

    const auto params = chen2020_params(); // initial_soc = 1.0，capacity = 5.0 Ah
    engine.start(params, proto);
    int waited = 0;
    while (engine.status() != SimulationStatus::Completed && waited < 3000) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        ++waited;
    }
    ASSERT_EQ(engine.status(), SimulationStatus::Completed);

    const auto snaps = engine.collect_snapshots();
    ASSERT_EQ(snaps.size(), 1201u); // t = 0..1200

    // 初始点 SOC = 初始 SOC（此前 bug：初始快照未经 push_snapshot，SOC 恒为 0）
    EXPECT_NEAR(snaps.front().soc, params.initial_soc, 1e-12);

    // 放电段单调下降；600 s 后 SOC = 1 − (5A·600s/3600)/5Ah ≈ 0.833
    //（此前 bug：吞吐量双重除以 3600，SOC 几乎不变）
    for (std::size_t i = 1; i <= 600; ++i) {
        EXPECT_LT(snaps[i].soc, snaps[i - 1].soc) << "at t=" << snaps[i].time;
    }
    const real_t expected_after_discharge =
        params.initial_soc - (5.0 * 600.0 / 3600.0) / params.capacity;
    EXPECT_NEAR(snaps[600].soc, expected_after_discharge, 1e-9);

    // 充电段单调回升（此前 bug：abs 使充电段继续下降）；对称充放后回到初始 SOC
    for (std::size_t i = 601; i < snaps.size(); ++i) {
        EXPECT_GT(snaps[i].soc, snaps[i - 1].soc) << "at t=" << snaps[i].time;
    }
    EXPECT_NEAR(snaps.back().soc, params.initial_soc, 1e-9);
    engine.stop();
}

TEST(SimulationEngine, MultiStepDischargeRestCharge) {
    SimulationEngine engine;
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    // 放电 10s → 静置 10s → 充电 10s，均按时间截止
    CCStep discharge;
    discharge.mode = CCStep::Mode::Discharge;
    discharge.current = 5.0;
    discharge.cutoff_time = 10.0;
    CCStep rest;
    rest.mode = CCStep::Mode::Rest;
    rest.cutoff_time = 10.0;
    CCStep charge;
    charge.mode = CCStep::Mode::Charge;
    charge.current = 5.0;
    charge.cutoff_time = 10.0;
    proto.steps = {discharge, rest, charge};

    engine.start(chen2020_params(), proto);
    int waited = 0;
    while (engine.status() != SimulationStatus::Completed && waited < 500) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        ++waited;
    }
    ASSERT_EQ(engine.status(), SimulationStatus::Completed);

    auto snaps = engine.collect_snapshots();
    ASSERT_EQ(snaps.size(), 31u); // t = 0..30
    EXPECT_NEAR(snaps.back().time, 30.0, 1e-9);

    // 段电流符号：放电 > 0，静置 = 0，充电 < 0
    EXPECT_GT(snaps[5].current, 0.0);
    EXPECT_NEAR(snaps[15].current, 0.0, 1e-12);
    EXPECT_LT(snaps[25].current, 0.0);

    // 静置段过电位弛豫：放电段末电压 < 静置段末电压
    EXPECT_LT(snaps[10].voltage, snaps[20].voltage);
    // 充电段电压回升
    EXPECT_GT(snaps[30].voltage, snaps[20].voltage);
}

TEST(SimulationEngine, MultiStepVoltageCutoffSwitchesToRest) {
    SimulationEngine engine;
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    // 5C 放电至 4.0V 截止（早于时间上限），随后静置 5s
    CCStep discharge;
    discharge.mode = CCStep::Mode::Discharge;
    discharge.current = 25.0;
    discharge.cutoff_time = 3600.0;
    discharge.cutoff_voltage = 4.0;
    CCStep rest;
    rest.mode = CCStep::Mode::Rest;
    rest.cutoff_time = 5.0;
    proto.steps = {discharge, rest};

    engine.start(chen2020_params(), proto);
    int waited = 0;
    while (engine.status() != SimulationStatus::Completed && waited < 3000) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        ++waited;
    }
    ASSERT_EQ(engine.status(), SimulationStatus::Completed);

    auto snaps = engine.collect_snapshots();
    // 电压先降（放电）后升（静置弛豫）→ 内部存在电压极小点
    const auto min_it = std::min_element(
        snaps.begin(), snaps.end(),
        [](const SolutionSnapshot& a, const SolutionSnapshot& b) {
            return a.voltage < b.voltage;
        });
    ASSERT_NE(min_it, snaps.begin());          // 确有下降
    ASSERT_NE(min_it, std::prev(snaps.end())); // 极小点不是最后一个
    EXPECT_GT(snaps.back().voltage, min_it->voltage);
    // 总时长不超过放电时间上限 + 静置段（电压截止先于时间上限触发）
    EXPECT_LE(snaps.back().time, 3600.0 + 5.0 + 1e-9);
}

TEST(SimulationEngine, MultiStepSingleStepCrossesSegments) {
    SimulationEngine engine;
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    CCStep discharge;
    discharge.mode = CCStep::Mode::Discharge;
    discharge.current = 5.0;
    discharge.cutoff_time = 3.0;
    CCStep rest;
    rest.mode = CCStep::Mode::Rest;
    rest.cutoff_time = 2.0;
    proto.steps = {discharge, rest};

    // 步 1-3：放电段，第 3 步后段截止切换到静置段
    for (int i = 0; i < 3; ++i) {
        engine.single_step(chen2020_params(), proto);
    }
    ASSERT_EQ(engine.status(), SimulationStatus::Paused);

    // 步 4-5：静置段，第 5 步后全部段完成
    engine.single_step(chen2020_params(), proto);
    EXPECT_EQ(engine.status(), SimulationStatus::Paused);
    engine.single_step(chen2020_params(), proto);
    EXPECT_EQ(engine.status(), SimulationStatus::Completed);

    auto snaps = engine.collect_snapshots();
    ASSERT_EQ(snaps.size(), 6u); // t = 0..5
    EXPECT_NEAR(snaps[3].current, 5.0, 1e-12); // t = 3，放电段最后一步
    EXPECT_NEAR(snaps[4].current, 0.0, 1e-12); // t = 4，静置段首步
    EXPECT_GT(snaps[5].voltage, snaps[2].voltage); // 静置弛豫回升
}

TEST(SimulationEngine, MultiStepValidation) {
    MultiStepProtocol proto;

    // 空 steps
    EXPECT_THROW(proto.validate(), std::invalid_argument);

    // 段缺少截止条件
    CCStep no_cutoff;
    no_cutoff.mode = CCStep::Mode::Discharge;
    no_cutoff.current = 5.0;
    no_cutoff.cutoff_time = 0.0;
    no_cutoff.cutoff_voltage = 0.0;
    proto.steps.push_back(no_cutoff);
    EXPECT_THROW(proto.validate(), std::invalid_argument);

    // 电流非正（非静置段）
    CCStep bad_current;
    bad_current.mode = CCStep::Mode::Discharge;
    bad_current.current = 0.0;
    bad_current.cutoff_time = 10.0;
    proto.steps[0] = bad_current;
    EXPECT_THROW(proto.validate(), std::invalid_argument);

    // 合法：时间截止
    proto.steps[0] = CCStep{}; // 默认放电 5A 3600s
    EXPECT_NO_THROW(proto.validate());

    // 合法：纯电压截止
    CCStep voltage_only;
    voltage_only.mode = CCStep::Mode::Charge;
    voltage_only.current = 2.0;
    voltage_only.cutoff_time = 0.0;
    voltage_only.cutoff_voltage = 4.2;
    proto.steps[0] = voltage_only;
    EXPECT_NO_THROW(proto.validate());
}
