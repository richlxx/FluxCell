/// @file cell_model.hpp
/// @brief 电化学模型抽象接口：SPM / SPMe / DFN 统一接入仿真引擎
#pragma once

#include <fluxcell/core/solver/dae_problem.hpp>
#include <fluxcell/core/types.hpp>

namespace fluxcell::model {

/// 电化学模型抽象接口
///
/// 仿真引擎仅依赖此接口驱动模型；具体模型（SPM、SPMe 等）负责各自的
/// 参数结构、初始化与时间步进。快照统一通过 SolutionSnapshot 输出，
/// 模型不支持的字段（如 SPM 的电解液分布）留空。
class CellModel {
public:
    virtual ~CellModel() = default;

    CellModel() = default;
    CellModel(const CellModel&) = delete;
    CellModel& operator=(const CellModel&) = delete;

    /// 以恒定电流推进一个时间步
    /// @param current  外加电流 [A]（正 = 放电）
    /// @param dt       时间步长 [s]
    virtual void step(real_t current, real_t dt) = 0;

    /// 当前端电压 [V]
    virtual real_t voltage() const = 0;

    /// 已仿真时间 [s]
    virtual real_t time() const = 0;

    /// 名义容量 [Ah]（SOC 库仑计数估算用）
    virtual real_t nominal_capacity() const = 0;

    /// 初始荷电状态 0~1（SOC 估算基线）
    virtual real_t initial_soc() const = 0;

    /// 生成当前状态的仿真快照
    /// @param current 本步工作电流 [A]（快照记录用）
    virtual SolutionSnapshot snapshot(real_t current) const = 0;

    /// DAE 问题视图：支持 SUNDIALS IDA 自适应求解的模型返回自身，
    /// 否则返回 nullptr
    virtual solver::DaeProblem* dae_problem() { return nullptr; }
};

} // namespace fluxcell::model
