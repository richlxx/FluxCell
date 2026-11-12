/// @file test_io.cpp
/// @brief 工况文件（save/load_protocol YAML 往返）单元测试
#include <fluxcell/core/io/io.hpp>
#include <fluxcell/core/simulation/simulation.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <stdexcept>

using namespace fluxcell;

namespace {

/// 临时目录内的唯一文件路径（测试结束自动清理由调用方 RAII 保证）
std::filesystem::path temp_yaml(const char* name) {
    return std::filesystem::temp_directory_path() / name;
}

} // namespace

// ══════════════════════════════════════════════════════════════
//  save_protocol → load_protocol 往返
// ══════════════════════════════════════════════════════════════

TEST(ProtocolIO, RoundTripFullConfiguration) {
    io::ProtocolFile pf;

    // 完整多段工况：CC-CV 充电（展开后形态）→ 静置 → CP 放电 → CV
    pf.protocol.steps = {
        {CCStep::Mode::Charge,    CCStep::Control::Current, 5.0, 4.2, 20.0,
         3600.0, 4.3, 0.05},
        {CCStep::Mode::Rest,      CCStep::Control::Current, 0.0, 4.2, 20.0,
         600.0, 0.0, 0.0},
        {CCStep::Mode::Discharge, CCStep::Control::Power,   5.0, 4.2, 20.0,
         7200.0, 2.5, 0.0},
        {CCStep::Mode::Discharge, CCStep::Control::Voltage, 5.0, 3.8, 20.0,
         0.0, 0.0, 0.05},
    };
    pf.protocol.cycles = 100;
    pf.protocol.time_step = 1.0;
    pf.protocol.real_time_factor = 0.0;
    pf.protocol.model = CellModelType::DFN;
    pf.protocol.solver = SolverType::SundialsIDA;
    pf.protocol.solver_rel_tol = 1e-6;
    pf.initial_soc = 0.0;
    pf.ambient_temp_c = 25.0;
    pf.thermal_type = model::ThermalType::Lumped;

    const auto path = temp_yaml("fluxcell_test_protocol_roundtrip.yml");
    io::save_protocol(path, pf);
    const io::ProtocolFile back = io::load_protocol(path);
    std::filesystem::remove(path);

    // 顶层运行条件
    EXPECT_EQ(back.protocol.model, CellModelType::DFN);
    EXPECT_EQ(back.protocol.solver, SolverType::SundialsIDA);
    EXPECT_DOUBLE_EQ(back.protocol.solver_rel_tol, 1e-6);
    EXPECT_EQ(back.thermal_type, model::ThermalType::Lumped);
    EXPECT_DOUBLE_EQ(back.initial_soc, 0.0);
    EXPECT_DOUBLE_EQ(back.ambient_temp_c, 25.0);
    EXPECT_EQ(back.protocol.cycles, 100);
    EXPECT_DOUBLE_EQ(back.protocol.time_step, 1.0);
    EXPECT_DOUBLE_EQ(back.protocol.real_time_factor, 0.0);

    // 逐段比对
    ASSERT_EQ(back.protocol.steps.size(), pf.protocol.steps.size());
    for (std::size_t i = 0; i < pf.protocol.steps.size(); ++i) {
        const CCStep& a = pf.protocol.steps[i];
        const CCStep& b = back.protocol.steps[i];
        EXPECT_EQ(b.mode, a.mode) << "step " << i;
        EXPECT_EQ(b.control, a.control) << "step " << i;
        EXPECT_DOUBLE_EQ(b.current, a.current) << "step " << i;
        EXPECT_DOUBLE_EQ(b.voltage, a.voltage) << "step " << i;
        EXPECT_DOUBLE_EQ(b.power, a.power) << "step " << i;
        EXPECT_DOUBLE_EQ(b.cutoff_time, a.cutoff_time) << "step " << i;
        EXPECT_DOUBLE_EQ(b.cutoff_voltage, a.cutoff_voltage) << "step " << i;
        EXPECT_DOUBLE_EQ(b.cutoff_current, a.cutoff_current) << "step " << i;
    }
}

TEST(ProtocolIO, RoundTripSpmDefaults) {
    // SPM + EigenCN + isothermal（默认口径）也能无损往返
    io::ProtocolFile pf;
    pf.protocol.steps = {
        {CCStep::Mode::Discharge, CCStep::Control::Current, 2.5, 4.2, 20.0,
         1800.0, 2.8, 0.0},
    };
    pf.protocol.model = CellModelType::SPM;
    pf.protocol.solver = SolverType::EigenCN;
    pf.thermal_type = model::ThermalType::Isothermal;
    pf.initial_soc = 0.5;
    pf.ambient_temp_c = 40.0;
    pf.protocol.cycles = 3;
    pf.protocol.time_step = 0.5;
    pf.protocol.real_time_factor = 10.0;

    const auto path = temp_yaml("fluxcell_test_protocol_spm.yml");
    io::save_protocol(path, pf);
    const io::ProtocolFile back = io::load_protocol(path);
    std::filesystem::remove(path);

    EXPECT_EQ(back.protocol.model, CellModelType::SPM);
    EXPECT_EQ(back.protocol.solver, SolverType::EigenCN);
    EXPECT_EQ(back.thermal_type, model::ThermalType::Isothermal);
    EXPECT_DOUBLE_EQ(back.initial_soc, 0.5);
    EXPECT_DOUBLE_EQ(back.ambient_temp_c, 40.0);
    EXPECT_EQ(back.protocol.cycles, 3);
    EXPECT_DOUBLE_EQ(back.protocol.time_step, 0.5);
    EXPECT_DOUBLE_EQ(back.protocol.real_time_factor, 10.0);
    ASSERT_EQ(back.protocol.steps.size(), 1u);
    EXPECT_EQ(back.protocol.steps[0].mode, CCStep::Mode::Discharge);
    EXPECT_EQ(back.protocol.steps[0].control, CCStep::Control::Current);
    EXPECT_DOUBLE_EQ(back.protocol.steps[0].current, 2.5);
    EXPECT_DOUBLE_EQ(back.protocol.steps[0].cutoff_time, 1800.0);
    EXPECT_DOUBLE_EQ(back.protocol.steps[0].cutoff_voltage, 2.8);
}

// ══════════════════════════════════════════════════════════════
//  容错：缺失键取默认、非法枚举报错、协议校验联动
// ══════════════════════════════════════════════════════════════

TEST(ProtocolIO, MissingKeysFallBackToDefaults) {
    const auto path = temp_yaml("fluxcell_test_protocol_empty.yml");
    {
        std::ofstream ofs(path);
        ofs << "# 只有段序列，其余键缺失\n";
        ofs << "steps:\n";
        ofs << "  - mode: charge\n";
        ofs << "    control: current\n";
        ofs << "    current: 1.5\n";
        ofs << "    cutoff_time: 900\n";
    }
    const io::ProtocolFile pf = io::load_protocol(path);
    std::filesystem::remove(path);

    // 缺失键取 ProtocolFile/MultiStepProtocol 默认值
    EXPECT_EQ(pf.protocol.model, CellModelType::SPM);
    EXPECT_EQ(pf.protocol.solver, SolverType::EigenCN);
    EXPECT_DOUBLE_EQ(pf.protocol.solver_rel_tol, 1e-6);
    EXPECT_EQ(pf.thermal_type, model::ThermalType::Isothermal);
    EXPECT_DOUBLE_EQ(pf.initial_soc, 0.0);
    EXPECT_DOUBLE_EQ(pf.ambient_temp_c, 25.0);
    EXPECT_EQ(pf.protocol.cycles, 1);
    ASSERT_EQ(pf.protocol.steps.size(), 1u);
    EXPECT_EQ(pf.protocol.steps[0].mode, CCStep::Mode::Charge);
    EXPECT_EQ(pf.protocol.steps[0].control, CCStep::Control::Current);
    EXPECT_DOUBLE_EQ(pf.protocol.steps[0].current, 1.5);
    // 段内其余字段取 CCStep 默认值
    EXPECT_DOUBLE_EQ(pf.protocol.steps[0].voltage, 4.2);
    EXPECT_DOUBLE_EQ(pf.protocol.steps[0].cutoff_voltage, 0.0);
}

TEST(ProtocolIO, InvalidEnumThrows) {
    const auto path = temp_yaml("fluxcell_test_protocol_bad.yml");

    const char* bad_mode = "steps:\n  - mode: floating\n";
    {
        std::ofstream ofs(path);
        ofs << bad_mode;
    }
    EXPECT_THROW(io::load_protocol(path), std::exception);

    const char* bad_control = "steps:\n  - mode: charge\n    control: pressure\n";
    {
        std::ofstream ofs(path);
        ofs << bad_control;
    }
    EXPECT_THROW(io::load_protocol(path), std::exception);

    const char* bad_model = "model: DFMN\n";
    {
        std::ofstream ofs(path);
        ofs << bad_model;
    }
    EXPECT_THROW(io::load_protocol(path), std::exception);

    const char* bad_solver = "solver: Euler\n";
    {
        std::ofstream ofs(path);
        ofs << bad_solver;
    }
    EXPECT_THROW(io::load_protocol(path), std::exception);

    const char* bad_thermal = "thermal_type: adiabatic\n";
    {
        std::ofstream ofs(path);
        ofs << bad_thermal;
    }
    EXPECT_THROW(io::load_protocol(path), std::exception);

    std::filesystem::remove(path);
}

TEST(ProtocolIO, LoadedProtocolPassesValidation) {
    // 往返后的协议应通过 validate（校验口径与引擎一致）
    io::ProtocolFile pf;
    pf.protocol.steps = {
        {CCStep::Mode::Charge,    CCStep::Control::Current, 5.0, 4.2, 20.0,
         3600.0, 4.3, 0.05},
        {CCStep::Mode::Rest,      CCStep::Control::Current, 0.0, 4.2, 20.0,
         600.0, 0.0, 0.0},
        {CCStep::Mode::Discharge, CCStep::Control::Current, 5.0, 4.2, 20.0,
         3600.0, 2.5, 0.0},
    };
    pf.protocol.cycles = 10;
    pf.protocol.time_step = 1.0;
    pf.protocol.model = CellModelType::DFN;
    pf.protocol.solver = SolverType::SundialsIDA;

    const auto path = temp_yaml("fluxcell_test_protocol_valid.yml");
    io::save_protocol(path, pf);
    const io::ProtocolFile back = io::load_protocol(path);
    std::filesystem::remove(path);

    EXPECT_NO_THROW(back.protocol.validate());
}
