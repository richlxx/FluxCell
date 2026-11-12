/// @file test_ida_solver.cpp
/// @brief SUNDIALS IDA 自适应求解器单元测试：与固定步长 Crank-Nicolson 的
///        电压一致性（SPM/SPMe/DFN）/ 自适应步长诊断 / 段电流切换 reinit
#include <fluxcell/core/io/io.hpp>
#include <fluxcell/core/model/dfn_model.hpp>
#include <fluxcell/core/model/spme_model.hpp>
#include <fluxcell/core/model/spm_model.hpp>
#include <fluxcell/core/simulation/simulation.hpp>
#include <fluxcell/core/solver/dae_problem.hpp>
#include <fluxcell/core/solver/ida_solver.hpp>

#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace fluxcell;
using fluxcell::model::DFNModel;
using fluxcell::model::SPMModel;
using fluxcell::model::SPMeModel;
using fluxcell::solver::DaeProblem;
using fluxcell::solver::IdaSolver;

namespace {

/// 用 IDA 自适应推进恒流段到 t_end（终态经 dae_load_state 写回模型）
template <typename ModelT>
real_t ida_final_voltage(ModelT& model, real_t current, real_t t_end,
                         real_t rel_tol) {
    DaeProblem* problem = model.dae_problem();
    EXPECT_NE(problem, nullptr);
    if (problem == nullptr) {
        return 0.0;
    }
    problem->dae_set_current(current);
    IdaSolver ida(*problem, 0.0, rel_tol);
    ida.advance_to(t_end);
    return model.voltage();
}

/// 固定步长 Crank-Nicolson 推进恒流段（DFN 内部按需细化子步）
template <typename ModelT>
real_t cn_final_voltage(ModelT& model, real_t current, real_t dt,
                        real_t t_end) {
    const int steps = static_cast<int>(t_end / dt + 0.5);
    for (int i = 0; i < steps; ++i) {
        model.step(current, dt);
    }
    return model.voltage();
}

} // namespace

// ══════════════════════════════════════════════════════════════
//  与 Crank-Nicolson 的一致性
// ══════════════════════════════════════════════════════════════

TEST(IdaSolverSpm, MatchesCrankNicolson) {
    SPMModel::Parameters params; // 默认值即 Chen2020
    const real_t i0 = 5.0;
    const real_t t_end = 600.0;

    SPMModel ida_model(params, i0);
    const real_t v_ida = ida_final_voltage(ida_model, i0, t_end, 1e-8);

    SPMModel cn_model(params, i0);
    const real_t v_cn = cn_final_voltage(cn_model, i0, 1.0, t_end);

    EXPECT_TRUE(std::isfinite(v_ida));
    EXPECT_NEAR(v_ida, v_cn, 5e-4);
}

TEST(IdaSolverSpme, MatchesCrankNicolson) {
    SPMeModel::Parameters params;
    const real_t i0 = 5.0;
    const real_t t_end = 600.0;

    SPMeModel ida_model(params, i0);
    const real_t v_ida = ida_final_voltage(ida_model, i0, t_end, 1e-8);

    SPMeModel cn_model(params, i0);
    const real_t v_cn = cn_final_voltage(cn_model, i0, 1.0, t_end);

    EXPECT_TRUE(std::isfinite(v_ida));
    EXPECT_NEAR(v_ida, v_cn, 5e-4);
}

TEST(IdaSolverDfn, MatchesCrankNicolsonShortRun) {
    DFNModel::Parameters params;
    const real_t i0 = 5.0;
    const real_t t_end = 300.0;

    DFNModel ida_model(params, i0);
    const real_t v_ida = ida_final_voltage(ida_model, i0, t_end, 1e-6);

    DFNModel cn_model(params, i0);
    const real_t v_cn = cn_final_voltage(cn_model, i0, 0.5, t_end);

    EXPECT_TRUE(std::isfinite(v_ida));
    EXPECT_NEAR(v_ida, v_cn, 3e-3);
}

TEST(IdaSolverDfn, SeiDegradationMatchesCrankNicolson) {
    // CN 固定步长 vs IDA 自适应：含 SEI 降解（生长 + 膜电阻 + 分流）
    // 的电压与 SEI 厚度一致性（DAE 微分状态段 off_sei 的回归护栏）
    DFNModel::Parameters params;
    params.degradation.enabled = true;
    const real_t i0 = 5.0;
    const real_t t_end = 300.0;

    DFNModel cn_model(params, i0);
    for (int i = 0; i < 30; ++i) {
        cn_model.step(i0, 10.0);
    }

    DFNModel ida_model(params, i0);
    const real_t v_ida = ida_final_voltage(ida_model, i0, t_end, 1e-8);

    EXPECT_NEAR(v_ida, cn_model.voltage(), 5e-4);
    const real_t growth_cn = cn_model.sei_thickness_avg()
                           - params.degradation.l_sei_0;
    EXPECT_GT(growth_cn, 1e-9);
    EXPECT_NEAR(ida_model.sei_thickness_avg() - params.degradation.l_sei_0,
                growth_cn, 0.05 * growth_cn);
}

// ══════════════════════════════════════════════════════════════
//  自适应步长诊断与段切换
// ══════════════════════════════════════════════════════════════

TEST(IdaSolver, AdaptiveStepCountAndDiagnostics) {
    SPMModel model(SPMModel::Parameters{}, 5.0);
    DaeProblem* problem = model.dae_problem();
    ASSERT_NE(problem, nullptr);

    problem->dae_set_current(5.0);
    IdaSolver ida(*problem, 0.0, 1e-8);
    ida.advance_to(600.0);

    // 平滑扩散问题：自适应内部步应远少于 600 个 1 s 输出间隔
    EXPECT_GT(ida.last_step(), 0.0);
    EXPECT_NEAR(ida.time(), 600.0, 1e-9);
    EXPECT_LT(ida.total_steps(), 300);
    EXPECT_TRUE(std::isfinite(model.voltage()));
}

TEST(IdaSolver, ReinitAfterCurrentChange) {
    SPMModel model(SPMModel::Parameters{}, 5.0);
    DaeProblem* problem = model.dae_problem();
    ASSERT_NE(problem, nullptr);

    // 恒流放电 300 s（极化建立，端电压低于开路）
    problem->dae_set_current(5.0);
    IdaSolver ida(*problem, 0.0, 1e-8);
    ida.advance_to(300.0);
    const real_t v_discharge = model.voltage();

    // 段切换为静置：reinit 丢弃历史后继续推进，极化弛豫、电压回升
    problem->dae_set_current(0.0);
    ida.reinit(model.time(), 1.0);
    ida.advance_to(600.0);
    const real_t v_rest = model.voltage();

    EXPECT_TRUE(std::isfinite(v_rest));
    EXPECT_GT(v_rest, v_discharge);
}

// ══════════════════════════════════════════════════════════════
//  引擎多段工况（复现 GUI 场景）
// ══════════════════════════════════════════════════════════════

namespace {

/// 运行引擎直至 Completed/Failed（限时），返回最终状态
SimulationStatus run_until_done(SimulationEngine& engine, int timeout_ms) {
    SimulationStatus st = engine.status();
    int waited = 0;
    while ((st == SimulationStatus::Running || st == SimulationStatus::Paused)
           && waited < timeout_ms) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        ++waited;
        st = engine.status();
    }
    return st;
}

/// 复现 GUI 三段工况：1C 放电至 2.5V → 静置 600 s → 2C 充电（IDA 自适应），
/// 返回最终状态；失败时通过 ADD_FAILURE 输出引擎错误信息
SimulationStatus run_cutoff_scenario(CellModelType model,
                                     SimulationEngine& engine) {
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    proto.real_time_factor = 0.0;
    proto.solver = SolverType::SundialsIDA;
    proto.model = model;

    CCStep discharge;
    discharge.mode = CCStep::Mode::Discharge;
    discharge.current = 5.0;
    discharge.cutoff_time = 3600.0;
    discharge.cutoff_voltage = 2.5;
    CCStep rest;
    rest.mode = CCStep::Mode::Rest;
    rest.current = 0.0;
    rest.cutoff_time = 600.0;
    CCStep charge;
    charge.mode = CCStep::Mode::Charge;
    charge.current = 10.0;
    charge.cutoff_time = 600.0;
    proto.steps = {discharge, rest, charge};

    // start 的参数仅用于 SPM 模型路径；SPMe/DFN 参数经协议携带
    engine.start(SPMModel::Parameters{}, proto);
    const SimulationStatus st = run_until_done(engine, 120000);
    if (st == SimulationStatus::Failed) {
        ADD_FAILURE() << "engine failed: " << engine.error_message();
    }
    return st;
}

} // namespace

TEST(IdaEngine, SPMDischargeCutoffRestThen2CCharge) {
    SimulationEngine engine;
    ASSERT_EQ(run_cutoff_scenario(CellModelType::SPM, engine),
              SimulationStatus::Completed);

    const auto snaps = engine.collect_snapshots();
    EXPECT_GE(snaps.size(), 100u);
    // 充电段电压应明显高于放电截止电压
    EXPECT_GT(snaps.back().voltage, 2.5);
    engine.stop();
}

TEST(IdaEngine, SPMeDischargeCutoffRestThen2CCharge) {
    SimulationEngine engine;
    ASSERT_EQ(run_cutoff_scenario(CellModelType::SPMe, engine),
              SimulationStatus::Completed);

    const auto snaps = engine.collect_snapshots();
    EXPECT_GE(snaps.size(), 100u);
    EXPECT_GT(snaps.back().voltage, 2.5);
    engine.stop();
}

TEST(IdaEngine, DFNDischargeCutoffRestThen2CCharge) {
    SimulationEngine engine;
    ASSERT_EQ(run_cutoff_scenario(CellModelType::DFN, engine),
              SimulationStatus::Completed);

    const auto snaps = engine.collect_snapshots();
    EXPECT_GE(snaps.size(), 100u);
    EXPECT_GT(snaps.back().voltage, 2.5);
    engine.stop();
}

// ══════════════════════════════════════════════════════════════
//  CC-CV 充电复现（GUI 场景）：chen2020 DFN(IDA)，从 0 SOC 完整充电
//  1C 恒流充电至 4.2V → 4.2V 恒压保持；
//  曾报告：SOC 冲到 ~110% 而端电压停在 ~4.18V（未锁定目标电压）
// ══════════════════════════════════════════════════════════════

TEST(IdaEngine, DFNChargeCCToCV) {
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    proto.real_time_factor = 0.0;
    proto.solver = SolverType::SundialsIDA;
    proto.model = CellModelType::DFN;
    proto.dfn_params = io::load_dfn_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params"
        / "dfn_chen2020.yml");
    proto.dfn_params.solid.initial_soc = 0.0; // 用户场景：从 0 SOC 开始充电
    // 本测试为等温基准（参数文件默认 lumped，显式覆盖保持历史行为）
    proto.dfn_params.thermal_type = model::ThermalType::Isothermal;

    CCStep cc;
    cc.mode = CCStep::Mode::Charge;
    cc.current = 5.0;             // 1C（5 Ah）
    cc.cutoff_voltage = 4.2;      // 充电截止上限
    cc.cutoff_time = 7200.0;
    CCStep cv;
    cv.mode = CCStep::Mode::Charge;
    cv.control = CCStep::Control::Voltage;
    cv.voltage = 4.2;
    cv.cutoff_time = 3600.0;      // 恒压保持 1 h
    proto.steps = {cc, cv};

    SimulationEngine engine;
    engine.start(SPMModel::Parameters{}, proto);
    const SimulationStatus st = run_until_done(engine, 600000);
    ASSERT_EQ(st, SimulationStatus::Completed)
        << "engine failed: " << engine.error_message();

    const auto snaps = engine.collect_snapshots();
    // 诊断输出：全部快照每 120 s 的状态（soc/voltage 自洽性）
    std::printf("\n---- 全程诊断（每 120 s）----\n");
    for (std::size_t i = 0; i < snaps.size(); i += 120) {
        const auto& s = snaps[i];
        std::printf("t=%7.0f s  V=%.4f  I=%+7.3f A  SOC=%6.2f%%\n",
                    s.time, s.voltage, s.current, s.soc * 100.0);
    }
    std::printf("末端: V=%.4f  I=%+7.4f  SOC=%.2f%%\n", snaps.back().voltage,
                snaps.back().current, snaps.back().soc * 100.0);

    // CV 段末：端电压应锁定在目标 ±5 mV
    EXPECT_NEAR(snaps.back().voltage, 4.2, 5e-3);
    // 从 0 SOC 充电：CC 段至 4.2V 约 86~90%，CV 段尾电流衰减后收敛在
    // OCV(4.2V) 对应值附近（<100%）；冲到 110% 即库仑计数/反解异常
    EXPECT_LT(snaps.back().soc, 1.05);
    engine.stop();
}

// ══════════════════════════════════════════════════════════════
//  4.4V 高电压恒压复现（GUI 场景）：报告 IDA+KLU 不收敛。
//  诊断测试：失败时也输出轨迹与引擎错误，定位 CV 段崩溃时刻与形态
// ══════════════════════════════════════════════════════════════

TEST(IdaEngine, DFNCCCV44VDiagnostic) {
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    proto.real_time_factor = 0.0;
    proto.solver = SolverType::SundialsIDA;
    proto.model = CellModelType::DFN;
    proto.dfn_params = io::load_dfn_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params"
        / "dfn_chen2020.yml");
    proto.dfn_params.solid.initial_soc = 0.0; // 用户场景：从空电完整充电
    // 保留参数文件默认 lumped 热模型（等温下 4.4V CV 已验证收敛）

    CCStep cc;
    cc.mode = CCStep::Mode::Charge;
    cc.current = 5.0;             // 1C（GUI 一键 CCCV 典型倍率）
    cc.cutoff_voltage = 4.4;
    cc.cutoff_time = 7200.0;
    CCStep cv;
    cv.mode = CCStep::Mode::Charge;
    cv.control = CCStep::Control::Voltage;
    cv.voltage = 4.4;
    cv.cutoff_time = 3600.0;
    cv.cutoff_current = 0.1;
    proto.steps = {cc, cv};

    SimulationEngine engine;
    engine.start(SPMModel::Parameters{}, proto);
    const SimulationStatus st = run_until_done(engine, 900000);

    // 无论成败都输出轨迹（失败时 dae_dump_state 已打印末状态摘要）
    const auto snaps = engine.collect_snapshots();
    std::printf("\n---- 4.4V CCCV 诊断（每 60 s）状态=%d ----\n",
                static_cast<int>(st));
    for (std::size_t i = 0; i < snaps.size(); i += 60) {
        const auto& s = snaps[i];
        std::printf("t=%7.0f s  V=%.4f  I=%+7.3f A  SOC=%6.2f%%\n",
                    s.time, s.voltage, s.current, s.soc * 100.0);
    }
    if (!snaps.empty()) {
        std::printf("末端: t=%.0f V=%.4f  I=%+7.4f  SOC=%.2f%%\n",
                    snaps.back().time, snaps.back().voltage,
                    snaps.back().current, snaps.back().soc * 100.0);
    }
    if (st == SimulationStatus::Failed) {
        ADD_FAILURE() << "engine failed: " << engine.error_message();
    }
    engine.stop();
}

/// 大偏差 CV 切换复现（GUI 报告「CCCV 4.4V 还是不收敛」的候选形态）：
/// CC 段时间截止（3000 s，V≈4.25）先于电压截止触发 → 恒压段初值
/// 离目标 ≈0.15 V，增广端口 CalcIC 需把 I_cell 从 −5 A 一口推到 ≈−12 A。
TEST(IdaEngine, DFNCCCV44VBigJumpSwitch) {
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    proto.real_time_factor = 0.0;
    proto.solver = SolverType::SundialsIDA;
    proto.model = CellModelType::DFN;
    proto.dfn_params = io::load_dfn_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params"
        / "dfn_chen2020.yml");
    proto.dfn_params.solid.initial_soc = 0.0;

    CCStep cc;
    cc.mode = CCStep::Mode::Charge;
    cc.current = 5.0;
    cc.cutoff_voltage = 4.4;
    cc.cutoff_time = 3000.0; // 时间截止先触发：V≈4.25 时切 CV
    CCStep cv;
    cv.mode = CCStep::Mode::Charge;
    cv.control = CCStep::Control::Voltage;
    cv.voltage = 4.4;
    cv.cutoff_time = 3600.0;
    cv.cutoff_current = 0.1;
    proto.steps = {cc, cv};

    SimulationEngine engine;
    engine.start(SPMModel::Parameters{}, proto);
    const SimulationStatus st = run_until_done(engine, 900000);

    const auto snaps = engine.collect_snapshots();
    std::printf("\n---- 4.4V CCCV 大偏差切换（每 60 s）状态=%d ----\n",
                static_cast<int>(st));
    for (std::size_t i = 0; i < snaps.size(); i += 60) {
        const auto& s = snaps[i];
        std::printf("t=%7.0f s  V=%.4f  I=%+7.3f A  SOC=%6.2f%%\n",
                    s.time, s.voltage, s.current, s.soc * 100.0);
    }
    if (!snaps.empty()) {
        std::printf("末端: t=%.0f V=%.4f  I=%+7.4f  SOC=%.2f%%\n",
                    snaps.back().time, snaps.back().voltage,
                    snaps.back().current, snaps.back().soc * 100.0);
    }
    if (st == SimulationStatus::Failed) {
        ADD_FAILURE() << "engine failed: " << engine.error_message();
    }
    engine.stop();
}

/// GUI 一键 CCCV 默认形态：StepUI 默认值（cutoff_time=3600 s、
/// 尾电流 0.05 A）展开的 CC+CV 两段，1C 0 SOC 充至 4.4 V。
/// 尾电流减半推进到更深的过充区（SOC ≈112.7%+），负极表面嵌锂率逼近 1。
TEST(IdaEngine, DFNCCCV44VGuiDefaults) {
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    proto.real_time_factor = 0.0;
    proto.solver = SolverType::SundialsIDA;
    proto.model = CellModelType::DFN;
    proto.dfn_params = io::load_dfn_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params"
        / "dfn_chen2020.yml");
    proto.dfn_params.solid.initial_soc = 0.0;

    CCStep cc;
    cc.mode = CCStep::Mode::Charge;
    cc.current = 5.0;
    cc.cutoff_voltage = 4.4;
    cc.cutoff_time = 3600.0; // StepUI 默认兑底时间
    CCStep cv;
    cv.mode = CCStep::Mode::Charge;
    cv.control = CCStep::Control::Voltage;
    cv.voltage = 4.4;
    cv.cutoff_time = 14400.0;
    cv.cutoff_current = 0.05; // StepUI 默认尾电流；时间兑底放开 → 走到
                              // |I|→0.05 A 的最深过充区（SOC >112.8%，
                              // 负极表面嵌锂率→1，BV 行 i0→0 极限）
    proto.steps = {cc, cv};

    SimulationEngine engine;
    engine.start(SPMModel::Parameters{}, proto);
    const SimulationStatus st = run_until_done(engine, 900000);

    const auto snaps = engine.collect_snapshots();
    std::printf("\n---- 4.4V CCCV GUI 默认形态深过充（每 60 s）状态=%d ----\n",
                static_cast<int>(st));
    for (std::size_t i = 0; i < snaps.size(); i += 120) {
        const auto& s = snaps[i];
        std::printf("t=%7.0f s  V=%.4f  I=%+7.4f A  SOC=%6.2f%%\n",
                    s.time, s.voltage, s.current, s.soc * 100.0);
    }
    if (!snaps.empty()) {
        std::printf("末端: t=%.0f V=%.4f  I=%+7.4f  SOC=%.2f%%\n",
                    snaps.back().time, snaps.back().voltage,
                    snaps.back().current, snaps.back().soc * 100.0);
    }
    if (st == SimulationStatus::Failed) {
        ADD_FAILURE() << "engine failed: " << engine.error_message();
    }
    engine.stop();
}

/// 用户实际三段工况复现：CCCV 1C 充至 4.4V → 静置 600s → 放电 1C 至 2.5V。
/// CV 段末负极表面嵌锂率贴 1（c_surf = c_max，exchange i0 → 0），
/// 静置→放电切换的 Gummel 重估从饱和表面出 5A → BV sinh 溢出。
TEST(IdaEngine, DFNTriStepCCCVRestDischarge) {
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    proto.real_time_factor = 0.0;
    proto.solver = SolverType::SundialsIDA;
    proto.model = CellModelType::DFN;
    proto.dfn_params = io::load_dfn_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params"
        / "dfn_chen2020.yml");
    proto.dfn_params.solid.initial_soc = 0.0;

    CCStep cc;
    cc.mode = CCStep::Mode::Charge;
    cc.current = 5.0;
    cc.cutoff_voltage = 4.4;
    cc.cutoff_time = 3600.0;
    CCStep cv;
    cv.mode = CCStep::Mode::Charge;
    cv.control = CCStep::Control::Voltage;
    cv.voltage = 4.4;
    cv.cutoff_time = 3600.0;
    cv.cutoff_current = 0.05;
    CCStep rest;
    rest.mode = CCStep::Mode::Rest;
    rest.cutoff_time = 600.0;
    CCStep dc;
    dc.mode = CCStep::Mode::Discharge;
    dc.current = 5.0;
    dc.cutoff_voltage = 2.5;
    dc.cutoff_time = 7200.0;
    proto.steps = {cc, cv, rest, dc};

    SimulationEngine engine;
    engine.start(SPMModel::Parameters{}, proto);
    const SimulationStatus st = run_until_done(engine, 900000);

    const auto snaps = engine.collect_snapshots();
    std::printf("\n---- CCCV+静置+放电 三段（每 120 s）状态=%d ----\n",
                static_cast<int>(st));
    for (std::size_t i = 0; i < snaps.size(); i += 120) {
        const auto& s = snaps[i];
        std::printf("t=%7.0f s  V=%.4f  I=%+7.4f A  SOC=%6.2f%%\n",
                    s.time, s.voltage, s.current, s.soc * 100.0);
    }
    if (!snaps.empty()) {
        std::printf("末端: t=%.0f V=%.4f  I=%+7.4f  SOC=%.2f%%\n",
                    snaps.back().time, snaps.back().voltage,
                    snaps.back().current, snaps.back().soc * 100.0);
    }
    if (st == SimulationStatus::Failed) {
        ADD_FAILURE() << "engine failed: " << engine.error_message();
    }
    engine.stop();
}

/// 用户完整循环工况：4.3V CCCV → 静置 600s → 1C 放电 2.5V × 100 圈
/// （GUI 全默认段参数，SOC 0 起）。此前用户在 BV 修复前的版本上第 3 圈
/// CV 段报 IDA Newton flag -4（饱和态试探步 sinh 溢出污染残差）；
/// 修复后本测试作为圈间稳定性回归基线
TEST(IdaEngine, DFNCycle100CCCV43VRestDischarge) {
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    proto.real_time_factor = 0.0;
    proto.solver = SolverType::SundialsIDA;
    proto.model = CellModelType::DFN;
    proto.dfn_params = io::load_dfn_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params"
        / "dfn_chen2020.yml");
    proto.dfn_params.solid.initial_soc = 0.0;

    CCStep cc;
    cc.mode = CCStep::Mode::Charge;
    cc.current = 5.0;
    cc.cutoff_voltage = 4.3;
    cc.cutoff_time = 3600.0;
    CCStep cv;
    cv.mode = CCStep::Mode::Charge;
    cv.control = CCStep::Control::Voltage;
    cv.voltage = 4.3;
    cv.cutoff_time = 3600.0;
    cv.cutoff_current = 0.05;
    CCStep rest;
    rest.mode = CCStep::Mode::Rest;
    rest.cutoff_time = 600.0;
    CCStep dc;
    dc.mode = CCStep::Mode::Discharge;
    dc.current = 5.0;
    dc.cutoff_voltage = 2.5;
    dc.cutoff_time = 3600.0; // GUI 默认兜底：1C 从 ~108% SOC 放到 2.5V
                             // 需 ~3800 s，3600 s 时间截止先触发（SOC ~8%）
    proto.steps = {cc, cv, rest, dc};
    proto.cycles = 100; // 用户完整工况：100 圈验证圈间稳定性

    SimulationEngine engine;
    engine.start(SPMModel::Parameters{}, proto);
    const SimulationStatus st = run_until_done(engine, 900000);

    const auto snaps = engine.collect_snapshots();
    std::printf("\n---- 4.3V 循环 x3（每 240 s）状态=%d ----\n",
                static_cast<int>(st));
    for (std::size_t i = 0; i < snaps.size(); i += 240) {
        const auto& s = snaps[i];
        std::printf("t=%7.0f s  V=%.4f  I=%+7.4f A  SOC=%6.2f%%\n",
                    s.time, s.voltage, s.current, s.soc * 100.0);
    }
    if (!snaps.empty()) {
        std::printf("末端: t=%.0f V=%.4f  I=%+7.4f  SOC=%.2f%%\n",
                    snaps.back().time, snaps.back().voltage,
                    snaps.back().current, snaps.back().soc * 100.0);
    }
    if (st == SimulationStatus::Failed) {
        ADD_FAILURE() << "engine failed: " << engine.error_message();
    }
    engine.stop();
}

/// 用户完整工况（dist/condition1.yml，200 圈）的前 5 圈：
/// CC 5A→4.25V、CV 4.25V→尾流 0.25A、静置 1600s、DC 5A→2.5V、静置 1800s；
/// 集总热 + SEI 降解开启。用户实测 t=27822 s（第 2~3 圈）报 IDA flag -4
///（V=4.25、I=-4.83 特征指向 CV 段）；本测试先以 5 圈判定复现性。
namespace {

/// 4.25V CV 循环复现主体（sei/lumped 双开关归因矩阵，圈数可调）
void cycle5_cv425_body(bool sei, bool lumped, int cycles = 5) {
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    proto.real_time_factor = 0.0;
    proto.solver = SolverType::SundialsIDA;
    proto.model = CellModelType::DFN;
    proto.dfn_params = io::load_dfn_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params"
        / "dfn_chen2020.yml");
    proto.dfn_params.solid.initial_soc = 0.0;
    // SEI 降解 = GUI「机理寿命」；lumped = YAML thermal.type（显式声明口径）
    proto.dfn_params.thermal_type = lumped ? model::ThermalType::Lumped
                                           : model::ThermalType::Isothermal;
    proto.dfn_params.degradation.enabled = sei;

    CCStep cc;   // 充电 5A → 4.25V 截止
    cc.mode = CCStep::Mode::Charge;
    cc.current = 5.0;
    cc.cutoff_voltage = 4.25;
    cc.cutoff_time = 3600.0;
    CCStep cv;   // 恒压 4.25V → 尾电流 0.25A
    cv.mode = CCStep::Mode::Charge;
    cv.control = CCStep::Control::Voltage;
    cv.voltage = 4.25;
    cv.cutoff_time = 3600.0;
    cv.cutoff_current = 0.25;
    CCStep rest1; // 静置 1600 s
    rest1.mode = CCStep::Mode::Rest;
    rest1.cutoff_time = 1600.0;
    CCStep dc;   // 放电 5A → 2.5V
    dc.mode = CCStep::Mode::Discharge;
    dc.current = 5.0;
    dc.cutoff_voltage = 2.5;
    dc.cutoff_time = 3600.0;
    CCStep rest2; // 静置 1800 s
    rest2.mode = CCStep::Mode::Rest;
    rest2.cutoff_time = 1800.0;
    proto.steps = {cc, cv, rest1, dc, rest2};
    proto.cycles = cycles;

    SimulationEngine engine;
    engine.start(SPMModel::Parameters{}, proto);
    const SimulationStatus st = run_until_done(engine, 900000);

    const auto snaps = engine.collect_snapshots();
    std::printf("\n---- 4.25V CV + SEI x5 圈（每 600 s）状态=%d ----\n",
                static_cast<int>(st));
    for (std::size_t i = 0; i < snaps.size(); i += 600) {
        const auto& s = snaps[i];
        std::printf("t=%7.0f s  V=%.4f  I=%+7.4f A  SOC=%6.2f%%\n",
                    s.time, s.voltage, s.current, s.soc * 100.0);
    }
    if (!snaps.empty()) {
        std::printf("末端: t=%.0f V=%.4f  I=%+7.4f  SOC=%.2f%%\n",
                    snaps.back().time, snaps.back().voltage,
                    snaps.back().current, snaps.back().soc * 100.0);
    }
    if (st == SimulationStatus::Failed) {
        ADD_FAILURE() << "engine failed: " << engine.error_message();
    }
    engine.stop();
}

} // namespace

TEST(IdaEngine, DFNCycle5CV425LumpedSei) {     // 用户口径：lumped + SEI
    cycle5_cv425_body(true, true);
}
TEST(IdaEngine, DFNCycle5CV425LumpedNoSei) {   // 归因：关 SEI
    cycle5_cv425_body(false, true);
}
TEST(IdaEngine, DFNCycle5CV425IsothermalSei) { // 归因：关集总热
    cycle5_cv425_body(true, false);
}
TEST(IdaEngine, DFNCycle5CV425IsothermalNoSei) {
    cycle5_cv425_body(false, false);
}
TEST(IdaEngine, DFNCycle20CV425LumpedSei) {    // 长程：CV 反复进出的鲁棒性
    cycle5_cv425_body(true, true, 20);
}

/// DFN 1C 恒流放电（满电起步，与 PyBaMM chen2020_dfn_discharge_1C.csv
/// 同起点同工况）：诊断 CC 段电压是否与基准对齐（充电偏差的对称性判别）
TEST(IdaEngine, DFNDischarge1CDiagnostic) {
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    proto.real_time_factor = 0.0;
    proto.solver = SolverType::SundialsIDA;
    proto.model = CellModelType::DFN;
    proto.dfn_params = io::load_dfn_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params"
        / "dfn_chen2020.yml"); // initial_soc = 1.0（满电，与放电基准同起点）
    // 与 PyBaMM chen2020_dfn_discharge_1C.csv 对齐（等温基准）：显式关闭
    // 参数文件默认的 lumped 热模型
    proto.dfn_params.thermal_type = model::ThermalType::Isothermal;

    CCStep d;
    d.mode = CCStep::Mode::Discharge;
    d.current = 5.0;
    d.cutoff_voltage = 2.5;
    d.cutoff_time = 7200.0;
    proto.steps = {d};

    SimulationEngine engine;
    engine.start(SPMModel::Parameters{}, proto);
    const SimulationStatus st = run_until_done(engine, 600000);
    ASSERT_EQ(st, SimulationStatus::Completed)
        << "engine failed: " << engine.error_message();

    const auto snaps = engine.collect_snapshots();
    std::printf("\n---- DFN 1C 放电（每 60 s）----\n");
    for (std::size_t i = 0; i < snaps.size(); i += 60) {
        const auto& s = snaps[i];
        std::printf("t=%7.0f s  V=%.4f  I=%+6.2f A\n",
                    s.time, s.voltage, s.current);
    }
    std::printf("末端: t=%.0f s  V=%.4f\n", snaps.back().time,
                snaps.back().voltage);
    engine.stop();
}

/// DFN 耦合热（lumped）1C 恒流放电（满电起步）：与 PyBaMM
/// chen2020_dfn_thermal_entropic_discharge_1C.csv（options {"thermal":
/// "lumped"}，Chen2020 热参数 + 用户测量熵变 dU/dT LUT）对标电压与
/// 电芯温度。参数文件 thermal 段即 lumped（含 CSV 熵变查找表），
/// 此处不覆盖 —— 端到端覆盖 IO thermal 段读取（CSV → EntropyLut）
/// + DAE 热行（逐格 sto 插值熵热）+ 温度反馈电化学。
TEST(IdaEngine, DFNThermalLumpedDischargeMatchesPyBaMM) {
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    proto.real_time_factor = 0.0;
    proto.solver = SolverType::SundialsIDA;
    proto.model = CellModelType::DFN;
    proto.dfn_params = io::load_dfn_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params"
        / "dfn_chen2020.yml"); // thermal 段默认 lumped
    ASSERT_EQ(proto.dfn_params.thermal_type, model::ThermalType::Lumped);
    // thermal 段的 CSV 熵变查找表应已加载（mV/K → V/K 换算抽查）
    ASSERT_FALSE(proto.dfn_params.thermal_params.entropic_lut_neg.empty());
    ASSERT_FALSE(proto.dfn_params.thermal_params.entropic_lut_pos.empty());
    EXPECT_NEAR(proto.dfn_params.thermal_params.entropic_lut_neg.dudt.front(),
                2.5e-5, 1e-12); // 0.025 mV/K → 2.5e-5 V/K

    CCStep d;
    d.mode = CCStep::Mode::Discharge;
    d.current = 5.0;
    d.cutoff_voltage = 2.5;
    d.cutoff_time = 3600.0;
    proto.steps = {d};

    SimulationEngine engine;
    engine.start(SPMModel::Parameters{}, proto);
    const SimulationStatus st = run_until_done(engine, 900000);
    ASSERT_EQ(st, SimulationStatus::Completed)
        << "engine failed: " << engine.error_message();

    const auto snaps = engine.collect_snapshots();
    // PyBaMM 26.8 DFN+lumped+熵变 LUT（Chen2020，1C → 2.5V）参考值
    // （生成脚本 python/validation/generate_thermal_reference.py）
    struct RefPoint { real_t t; real_t v; real_t T; };
    constexpr RefPoint kRefs[] = {
        {600.0,  3.82725, 303.5041},
        {1200.0, 3.67889, 306.0700},
        {1800.0, 3.52976, 306.4512},
        {2400.0, 3.41104, 306.8667},
        {3000.0, 3.24730, 308.7782},
    };
    for (const auto& ref : kRefs) {
        // 快照按 dt=1s 采样，线性扫描取首个 time >= ref.t 的快照
        std::size_t k = 0;
        while (k < snaps.size() && snaps[k].time < ref.t) {
            ++k;
        }
        ASSERT_LT(k, snaps.size()) << "missing snapshot at t=" << ref.t;
        const auto& snap = snaps[k];
        std::printf("t=%5.0f  V=%.4f (ref %.4f)  T=%.3f K (ref %.3f)\n",
                    snap.time, snap.voltage, ref.v,
                    snap.temperature, ref.T);
        // 容差覆盖：CN/IDA 与 PyBaMM BDF 的离散差异 + 温度经 Arrhenius
        // 反馈放大的电化学偏差（电压 20 mV / 温度 1.5 K）
        EXPECT_NEAR(snap.voltage, ref.v, 0.02) << "t=" << ref.t;
        EXPECT_NEAR(snap.temperature, ref.T, 1.5) << "t=" << ref.t;
    }
    engine.stop();
}

/// DFN + SEI 降解 2 圈 1C 放电/充电循环（满电起步）：与 PyBaMM
/// okane2022_dfn_sei_cc_2cycles.csv 对标（options {"SEI": "ec reaction
/// limited", "SEI film resistance": "distributed", "thermal":
/// "isothermal"}，OKane2022 参数）。dfn_chen2020.yml 的电化学部分 =
/// Chen2020（OKane2022 沿用），degradation 段 SEI 参数 = OKane2022；
/// 等温时其 Arrhenius 输运退化为常数基准，固相 Bruggeman 1.5 与
/// σ_eff = σ·ε^1.5 一致。圈末（cycle_end）锂损失/SEI 厚度为寿命
/// 曲线采样点（生成脚本
/// python/validation/generate_degradation_reference.py）
TEST(IdaEngine, DFNSeiDegradationCyclesMatchPyBaMM) {
    MultiStepProtocol proto;
    // 采样步 1 s（参考 CSV 为 30 s 采样）：引擎的电压截止在采样步边界
    // 检测（step_finished），30 s 步长下每段截止最多晚一个步长，两圈
    // 4 段累积 87~142 s 错位经圈 2 充电末段 ~0.31 mV/s 的陡峭斜率放大
    // 成 30~41 mV 假性电压偏差（电压轨迹本身仅差 <10 mV）；1 s 步长
    // 把截止错位压到秒级，使对比反映真实求解差异
    proto.time_step = 1.0;
    proto.real_time_factor = 0.0;
    proto.solver = SolverType::SundialsIDA;
    proto.model = CellModelType::DFN;
    proto.dfn_params = io::load_dfn_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params"
        / "dfn_chen2020.yml");
    proto.dfn_params.thermal_type = model::ThermalType::Isothermal;
    proto.dfn_params.degradation.enabled = true;
    proto.cycles = 2;

    CCStep d;
    d.mode = CCStep::Mode::Discharge;
    d.current = 5.0;         // 1C
    d.cutoff_voltage = 2.5;
    d.cutoff_time = 7200.0;
    CCStep c;
    c.mode = CCStep::Mode::Charge;
    c.current = 5.0;         // 1C
    c.cutoff_voltage = 4.2;
    c.cutoff_time = 7200.0;
    proto.steps = {d, c};

    SimulationEngine engine;
    engine.start(SPMModel::Parameters{}, proto);
    const SimulationStatus st = run_until_done(engine, 1200000);
    ASSERT_EQ(st, SimulationStatus::Completed)
        << "engine failed: " << engine.error_message();

    const auto snaps = engine.collect_snapshots();

    // 读取参考 CSV（time_s,voltage_V,li_loss_ah,sei_thickness_m,cycle_end）
    struct RefRow {
        real_t t = 0.0;
        real_t v = 0.0;
        real_t loss = 0.0;
        real_t l_sei = 0.0;
        bool cycle_end = false;
    };
    std::vector<RefRow> refs;
    {
        std::ifstream csv(std::filesystem::path(FLUXCELL_SOURCE_DIR)
                          / "python" / "validation" / "reference_data"
                          / "okane2022_dfn_sei_cc_2cycles.csv");
        ASSERT_TRUE(csv.is_open())
            << "参考 CSV 缺失（先运行 python/validation/"
               "generate_degradation_reference.py）";
        std::string line;
        std::getline(csv, line); // header
        while (std::getline(csv, line)) {
            if (line.empty()) {
                continue;
            }
            std::istringstream ss(line);
            std::string cell;
            RefRow row;
            std::getline(ss, cell, ',');
            row.t = static_cast<real_t>(std::stod(cell));
            std::getline(ss, cell, ',');
            row.v = static_cast<real_t>(std::stod(cell));
            std::getline(ss, cell, ',');
            row.loss = static_cast<real_t>(std::stod(cell));
            std::getline(ss, cell, ',');
            row.l_sei = static_cast<real_t>(std::stod(cell));
            std::getline(ss, cell, ',');
            row.cycle_end = cell == "1";
            refs.push_back(row);
        }
    }
    ASSERT_FALSE(refs.empty());

    // 全程电压对比（每 20 行抽 1 点；容差覆盖 CN/IDA 与 Casadi BDF
    // 的离散差异 + 段截止时刻的秒级错位 × 电压斜率）。
    // 段切换邻域（±400 s）的对比点跳过：截止前陡坡区（负极表面
    // 近空）对固相 LUT 插值/发散保护处理敏感，截止时刻秒级错位
    // 在陡坡（~12 mV/s）内放大为几十至上百 mV；切段后界面极化
    // 翻转 + 浓差分布慢弛豫（负极固相扩散时标 ~2240 s）需 ~400 s
    // 收敛 —— 该邻域的固有求解器差异（数十 mV）不代表稳态电化学
    // 精度，稳态区电压差 <10 mV、圈末寿命量差 <1% 为主要验证目标
    constexpr real_t kSwitchTimes[] = {3515.6, 5823.45, 8103.63, 10358.18};
    auto near_switch = [&](real_t t) {
        for (const real_t ts : kSwitchTimes) {
            if (std::abs(t - ts) < 400.0) {
                return true;
            }
        }
        return false;
    };
    for (std::size_t r = 0; r < refs.size(); r += 20) {
        const auto& ref = refs[r];
        if (near_switch(ref.t)) {
            continue;
        }
        std::size_t k = 0;
        while (k < snaps.size() && snaps[k].time < ref.t) {
            ++k;
        }
        ASSERT_LT(k, snaps.size()) << "missing snapshot at t=" << ref.t;
        EXPECT_NEAR(snaps[k].voltage, ref.v, 0.03) << "t=" << ref.t;
    }

    // 圈末寿命量对比（cycle_end 快照 vs 参考圈末采样窗口行）：
    // 锂损失与 SEI 厚度。容差 10%（相对增长量）：SEI 动力学指数敏感
    // （eta_SEI 差 10 mV → j_sei 差 ~20%），两求解器圈末时刻的
    // 秒级错位即引入同量级偏差
    std::vector<const SolutionSnapshot*> ends;
    std::vector<const RefRow*> ref_ends;
    for (const auto& s : snaps) {
        if (s.cycle_end) {
            ends.push_back(&s);
        }
    }
    for (const auto& r : refs) {
        if (r.cycle_end) {
            ref_ends.push_back(&r);
        }
    }
    ASSERT_EQ(ends.size(), ref_ends.size());
    for (std::size_t i = 0; i < ends.size(); ++i) {
        const auto& snap = *ends[i];
        const auto& ref = *ref_ends[i];
        std::printf(
            "圈 %d 末: t=%.0f V=%.4f loss=%.5f L=%.4e | ref t=%.0f "
            "V=%.4f loss=%.5f L=%.4e\n",
            static_cast<int>(i) + 1, snap.time, snap.voltage,
            snap.li_loss_ah, snap.sei_thickness, ref.t, ref.v,
            ref.loss, ref.l_sei);
        EXPECT_NEAR(snap.voltage, ref.v, 0.03);
        EXPECT_NEAR(snap.li_loss_ah, ref.loss, 0.10 * ref.loss);
        EXPECT_NEAR(snap.sei_thickness, ref.l_sei,
                    0.10 * (ref.l_sei - 5e-9) + 2e-9);
    }
    engine.stop();
}

/// SPM 从 0 SOC 1C 恒流充电（与 PyBaMM SPM 同工况）：最简模型隔离层 ——
/// 同步输出表面嵌入率（与 PyBaMM surface stoichiometry 对比），
/// 隔离“扩散响应差异”与“OCP/BV 公式差异”
TEST(IdaEngine, SPMCharge1CDiagnostic) {
    auto params = io::load_spm_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params"
        / "chen2020.yml");
    params.initial_soc = 0.0;
    model::SPMModel model(params, -5.0); // 1C 充电（放电为正）

    std::printf("\n---- SPM 1C 充电（每 120 s）----\n");
    std::printf("t=%7d s  V=%.4f  x_surf=%.4f  y_surf=%.4f\n", 0,
                model.voltage(), model.neg_stoichiometry(),
                model.pos_stoichiometry());
    for (int s = 1; s <= 3300; ++s) {
        model.step(-5.0, 1.0);
        if (s % 120 == 0 || model.voltage() >= 4.2) {
            std::printf("t=%7d s  V=%.4f  x_surf=%.4f  y_surf=%.4f\n", s,
                        model.voltage(), model.neg_stoichiometry(),
                        model.pos_stoichiometry());
        }
        if (model.voltage() >= 4.2) {
            break;
        }
    }
}

/// SPMe 从 0 SOC 1C 恒流充电（与 PyBaMM SPMe 同工况）：隔离 DFN 特有项
/// （反应分布/固相电位场）与共有层面（固相扩散/电解液/BV/OCV）
TEST(IdaEngine, SPMeCharge1CDiagnostic) {
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    proto.real_time_factor = 0.0;
    proto.solver = SolverType::SundialsIDA;
    proto.model = CellModelType::SPMe;
    proto.spme_params = io::load_spme_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params"
        / "spme_chen2020.yml");
    proto.spme_params.solid.initial_soc = 0.0;

    CCStep cc;
    cc.mode = CCStep::Mode::Charge;
    cc.current = 5.0;
    cc.cutoff_voltage = 4.2;
    cc.cutoff_time = 7200.0;
    proto.steps = {cc};

    SimulationEngine engine;
    engine.start(SPMModel::Parameters{}, proto);
    const SimulationStatus st = run_until_done(engine, 600000);
    ASSERT_EQ(st, SimulationStatus::Completed)
        << "engine failed: " << engine.error_message();

    const auto snaps = engine.collect_snapshots();
    std::printf("\n---- SPMe 1C 充电（每 120 s）----\n");
    for (std::size_t i = 0; i < snaps.size(); i += 120) {
        const auto& s = snaps[i];
        std::printf("t=%7.0f s  V=%.4f  I=%+6.2f A\n",
                    s.time, s.voltage, s.current);
    }
    std::printf("末端: t=%.0f s  V=%.4f\n", snaps.back().time,
                snaps.back().voltage);
    engine.stop();
}

/// 回归：500 圈 1C 充放电在圈 8 放电切换时“反应分布迭代未收敛”
/// （固定阻尼的 Gummel 定点迭代在 SEI 膜电阻反馈环增益随 L/温升
/// 增大后振荡失效，详见 dfn_model.cpp solve_reaction_distribution
/// 的自适应阻尼注释）。修复后 10 圈应全程收敛；圈末打印同时
/// 充当多圈寿命趋势的人工检查输出
TEST(IdaEngine, DFNLongCycleDivergenceDiagnostic) {
    MultiStepProtocol proto;
    proto.time_step = 10.0;
    proto.real_time_factor = 0.0;
    proto.solver = SolverType::SundialsIDA;
    proto.model = CellModelType::DFN;
    proto.dfn_params = io::load_dfn_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params"
        / "dfn_chen2020.yml"); // 默认 lumped 热（GUI 加载口径）
    proto.dfn_params.degradation.enabled = true;
    proto.cycles = 10;

    CCStep d;
    d.mode = CCStep::Mode::Discharge;
    d.current = 5.0;         // 1C
    d.cutoff_voltage = 2.5;
    d.cutoff_time = 7200.0;
    CCStep c;
    c.mode = CCStep::Mode::Charge;
    c.current = 5.0;         // 1C
    c.cutoff_voltage = 4.2;
    c.cutoff_time = 7200.0;
    proto.steps = {d, c};

    SimulationEngine engine;
    // 单步逐步推进（异常在调用线程同步处理，error_ 无竞态窗口），
    // 同时避开后台线程模式下 Failed 时 error_ 读空的疑似问题
    engine.single_step(SPMModel::Parameters{}, proto);
    int guard = 0;
    while (engine.status() == SimulationStatus::Paused && guard < 200000) {
        engine.single_step(SPMModel::Parameters{}, proto);
        ++guard;
    }
    const SimulationStatus st = engine.status();
    std::printf("\nstatus=%d  steps=%d  err=%s\n", static_cast<int>(st),
                guard, engine.error_message().c_str());

    const auto snaps = engine.collect_snapshots();
    std::printf("---- 每圈末 ----\n");
    for (const auto& s : snaps) {
        if (s.cycle_end) {
            std::printf(
                "圈 %2d 末 t=%7.0f V=%.4f loss=%.5f L=%.3e soc=%.4f T=%.2f\n",
                s.cycle_number, s.time, s.voltage, s.li_loss_ah,
                s.sei_thickness, s.soc, s.temperature);
        }
    }
    std::printf("---- 末尾 12 个快照 ----\n");
    const std::size_t tail = snaps.size() >= 12 ? snaps.size() - 12 : 0;
    for (std::size_t i = tail; i < snaps.size(); ++i) {
        const auto& s = snaps[i];
        std::printf("t=%7.0f V=%.4f I=%+6.2f soc=%.4f L=%.3e T=%.2f\n",
                    s.time, s.voltage, s.current, s.soc, s.sei_thickness,
                    s.temperature);
    }
    engine.stop();
    EXPECT_EQ(st, SimulationStatus::Completed);
}

/// 回归：厚 SEI（l_sei_0 = 1.5 µm，膜压降 ~0.45 V @1C，对应 ~1000 圈
/// 1C 循环的生长量级）IDA 长循环稳定性（段切换 + 充电深段）。
///
/// 已知限制（细网格）：默认网格（n≈2105）下 IDA 用 SPGMR 无预条件
/// 迭代求解器 + 差分 Jv（vcpkg sundials 未编入 KLU），厚膜工况的
/// corrector Newton 在圈 2~3 充电段反复失败（flag -4，SUNDIALS 自报
/// "poor iterative algorithm performance"，V 呈 ±30~100 mV/s 锯齿
/// 重启动伪瞬态）。判别实验（dense LU vs GMRES 同工况对照）锁定
/// 根因为线性求解器而非 BDF/CalcIC/方程一致性：
/// - DAE 方程一致性已由残差注入验证（Gummel 解的 |res| ≤ 5e-9；
///   此前 dae_pack_state 液相源漏 j_sei 导致 ~1.7 mol/m³/s 初值
///   残差，已修复）；
/// - 粗网格（n≈260，dense LU）同工况三圈全过（本测试）。
/// maxl 提升到 300 无效——无预条件迭代法对此系统本质不足。
/// 根治已落地：sundials[klu]（项目 overlay port）+ 稀疏 Jacobian
/// （逐列扰动探测模式 + 距离-1 着色差分填充 CSC + KLU 直接分解，见
/// ida_solver.cpp），细网格恢复验证见 DFNThickSeiKluSparseDirect。
/// 本测试保留为 dense LU 路径的厚膜回归护栏（n≤300 路径不变）。
/// 在此之前细网格长循环用 SolverType::EigenCN（Gummel 路径全程稳健，
/// 见 test_dfn 回归）。
TEST(IdaEngine, DFNThickSeiCoarseGridDenseLU) {
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    proto.real_time_factor = 0.0;
    proto.solver = SolverType::SundialsIDA;
    proto.model = CellModelType::DFN;
    proto.dfn_params = io::load_dfn_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params"
        / "dfn_chen2020.yml");
    // 粗网格：径向 6 × 厚度 36 → dae_size ≈ 260 < 300（dense LU）
    proto.dfn_params.solid.particle_cells = 6;
    proto.dfn_params.electrolyte_cells = 36;
    proto.dfn_params.thermal_type = model::ThermalType::Isothermal;
    proto.dfn_params.degradation.enabled = true;
    proto.dfn_params.degradation.l_sei_0 = 1.5e-6;
    proto.cycles = 3;

    CCStep d;
    d.mode = CCStep::Mode::Discharge;
    d.current = 5.0;         // 1C
    d.cutoff_voltage = 2.5;
    d.cutoff_time = 7200.0;
    CCStep c;
    c.mode = CCStep::Mode::Charge;
    c.current = 5.0;         // 1C
    c.cutoff_voltage = 4.2;
    c.cutoff_time = 7200.0;
    proto.steps = {d, c};

    SimulationEngine engine;
    engine.start(SPMModel::Parameters{}, proto);
    const SimulationStatus st = run_until_done(engine, 600000);
    if (st != SimulationStatus::Completed) {
        std::printf("coarse: err=%s\n", engine.error_message().c_str());
    }
    ASSERT_EQ(st, SimulationStatus::Completed);
    const auto snaps = engine.collect_snapshots();
    int ends = 0;
    for (const auto& s : snaps) {
        if (s.cycle_end) {
            ++ends;
            std::printf("coarse 圈 %d 末 t=%.0f V=%.4f L=%.4e\n",
                        s.cycle_number, s.time, s.voltage, s.sei_thickness);
        }
    }
    EXPECT_EQ(ends, 3);
    engine.stop();
}

/// KLU 稀疏直解根治验证：细网格（默认网格 n≈2105 > 300）+ 厚 SEI
///（l_sei_0 = 1.5 µm）同工况三圈——此前 SPGMR 无预条件 GMRES 在圈
/// 2~3 充电段反复崩溃（flag -4），sundials[klu]（overlay port）+
/// 稀疏 Jacobian 差分填充接入后应完整跑通（PyBaMM IDAKLUSolver
/// 同构），验证细网格 IDA 厚膜长循环能力恢复
TEST(IdaEngine, DFNThickSeiKluSparseDirect) {
#if !defined(FLUXCELL_WITH_SUNDIALS_KLU)
    GTEST_SKIP() << "未编译 SUNDIALS KLU（需项目 overlay port 的 "
                    "sundials 端口），细网格回退 GMRES 已知厚膜崩溃";
#else
    MultiStepProtocol proto;
    proto.time_step = 1.0;
    proto.real_time_factor = 0.0;
    proto.solver = SolverType::SundialsIDA;
    proto.model = CellModelType::DFN;
    proto.dfn_params = io::load_dfn_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params"
        / "dfn_chen2020.yml");
    proto.dfn_params.thermal_type = model::ThermalType::Isothermal;
    proto.dfn_params.degradation.enabled = true;
    proto.dfn_params.degradation.l_sei_0 = 1.5e-6;
    proto.cycles = 3;

    CCStep d;
    d.mode = CCStep::Mode::Discharge;
    d.current = 5.0;         // 1C
    d.cutoff_voltage = 2.5;
    d.cutoff_time = 7200.0;
    CCStep c;
    c.mode = CCStep::Mode::Charge;
    c.current = 5.0;         // 1C
    c.cutoff_voltage = 4.2;
    c.cutoff_time = 7200.0;
    proto.steps = {d, c};

    SimulationEngine engine;
    engine.start(SPMModel::Parameters{}, proto);
    const SimulationStatus st = run_until_done(engine, 600000);
    if (st != SimulationStatus::Completed) {
        std::printf("klu: err=%s\n", engine.error_message().c_str());
    }
    ASSERT_EQ(st, SimulationStatus::Completed);
    const auto snaps = engine.collect_snapshots();
    int ends = 0;
    real_t l_last = 0.0;
    for (const auto& s : snaps) {
        if (s.cycle_end) {
            ++ends;
            l_last = s.sei_thickness;
            std::printf("klu 圈 %d 末 t=%.0f V=%.4f L=%.4e\n",
                        s.cycle_number, s.time, s.voltage, s.sei_thickness);
        }
    }
    EXPECT_EQ(ends, 3);
    // 厚膜继续生长：圈末厚度应高于初值
    EXPECT_GT(l_last, 1.5e-6);
    engine.stop();
#endif
}

/// CN 路径 50 圈长循环对标：与 PyBaMM okane2022_dfn_sei_cc1c_50cycles.csv
///（generate_sei_cycles_reference.py --cycles 50，1C CC 循环、
/// 2.5~4.2 V、60 s 采样）逐圈对比圈末锂损失与 SEI 厚度。细网格
///（n≈2105）+ SolverType::EigenCN（固定步长 + Gummel 逐格精确解膜环）
/// 全程稳健，为长循环（2000 圈目标）的稳定路径；细网格 IDA 厚膜
/// 能力已由 KLU 稀疏直解恢复（DFNThickSeiKluSparseDirect），见
/// DFNThickSeiCoarseGridDenseLU 注释。输出步
/// 2 s：电压截止错位 ≤2 s（陡坡区 ~12 mV/s → 局部 ~24 mV）不影响
/// 圈末寿命量（j_sei ~1e-3 A/m² 慢变量的圈累积），故仅对比圈末量
TEST(DfnCycles, SeiDegradationCn50MatchPyBaMM) {
    // 长循环对标耗时分钟级（Debug 构建下更长），默认跳过；
    // 设 FLUXCELL_LONG_TESTS=1 启用（夜间/发版前验证）
    if (std::getenv("FLUXCELL_LONG_TESTS") == nullptr) {
        GTEST_SKIP() << "长循环对标：设 FLUXCELL_LONG_TESTS=1 启用";
    }
    MultiStepProtocol proto;
    proto.time_step = 2.0;
    proto.real_time_factor = 0.0;
    proto.solver = SolverType::EigenCN;
    proto.model = CellModelType::DFN;
    proto.dfn_params = io::load_dfn_parameters(
        std::filesystem::path(FLUXCELL_SOURCE_DIR) / "params"
        / "dfn_chen2020.yml");
    proto.dfn_params.thermal_type = model::ThermalType::Isothermal;
    proto.dfn_params.degradation.enabled = true;
    proto.cycles = 50;

    CCStep d;
    d.mode = CCStep::Mode::Discharge;
    d.current = 5.0;         // 1C
    d.cutoff_voltage = 2.5;
    d.cutoff_time = 7200.0;
    CCStep c;
    c.mode = CCStep::Mode::Charge;
    c.current = 5.0;         // 1C
    c.cutoff_voltage = 4.2;
    c.cutoff_time = 7200.0;
    proto.steps = {d, c};

    SimulationEngine engine;
    engine.start(SPMModel::Parameters{}, proto);
    const SimulationStatus st = run_until_done(engine, 1200000);
    if (st != SimulationStatus::Completed) {
        std::printf("cn50: err=%s\n", engine.error_message().c_str());
    }
    ASSERT_EQ(st, SimulationStatus::Completed);
    const auto snaps = engine.collect_snapshots();

    // 参考圈末行（loss, L_sei）
    struct RefEnd {
        real_t loss = 0.0;
        real_t l_sei = 0.0;
    };
    std::vector<RefEnd> refs;
    {
        std::ifstream csv(std::filesystem::path(FLUXCELL_SOURCE_DIR)
                          / "python" / "validation" / "reference_data"
                          / "okane2022_dfn_sei_cc1c_50cycles.csv");
        ASSERT_TRUE(csv.is_open()) << "参考 CSV 缺失";
        std::string line;
        std::getline(csv, line); // header
        while (std::getline(csv, line)) {
            if (line.empty()) {
                continue;
            }
            std::istringstream ss(line);
            std::string cell;
            std::getline(ss, cell, ',');           // time_s
            std::getline(ss, cell, ',');           // voltage_V
            RefEnd row;
            std::getline(ss, cell, ',');
            row.loss = static_cast<real_t>(std::stod(cell));
            std::getline(ss, cell, ',');
            row.l_sei = static_cast<real_t>(std::stod(cell));
            std::getline(ss, cell, ',');           // cycle_end
            if (cell == "1") {
                refs.push_back(row);
            }
        }
    }
    ASSERT_EQ(refs.size(), 50u) << "参考圈末行数不符";

    std::vector<const SolutionSnapshot*> ends;
    for (const auto& s : snaps) {
        if (s.cycle_end) {
            ends.push_back(&s);
        }
    }
    ASSERT_EQ(ends.size(), refs.size());

    // 逐圈对比：锂损失容差 10%（j_sei 对 eta_SEI 的 Tafel 指数敏感，
    // 段截止时刻秒级错位 + Casadi/CN 离散差异的圈级累积）；SEI 厚度
    // 相对增长量 10%
    for (std::size_t i = 0; i < refs.size(); ++i) {
        EXPECT_NEAR(ends[i]->li_loss_ah, refs[i].loss,
                    0.10 * refs[i].loss)
            << "圈 " << i + 1 << " 锂损失";
        EXPECT_NEAR(ends[i]->sei_thickness, refs[i].l_sei,
                    0.10 * (refs[i].l_sei - 5e-9))
            << "圈 " << i + 1 << " SEI 厚度";
    }
    std::printf("cn50: 圈1 loss=%.4f/%.4f L=%.3e/%.3e | 圈50 "
                "loss=%.4f/%.4f L=%.3e/%.3e\n",
                ends[0]->li_loss_ah, refs[0].loss,
                ends[0]->sei_thickness, refs[0].l_sei,
                ends[49]->li_loss_ah, refs[49].loss,
                ends[49]->sei_thickness, refs[49].l_sei);
    engine.stop();
}
