/// @file simulation.hpp
/// @brief 仿真引擎：后台线程驱动模型，线程安全快照队列供 GUI 轮询
#pragma once

#include <fluxcell/core/model/cell_model.hpp>
#include <fluxcell/core/model/dfn_model.hpp>
#include <fluxcell/core/model/spm_model.hpp>
#include <fluxcell/core/model/spme_model.hpp>
#include <fluxcell/core/types.hpp>

#if defined(FLUXCELL_WITH_SUNDIALS)
#include <fluxcell/core/solver/ida_solver.hpp>
#endif

#include <atomic>
#include <cstddef>
#include <deque>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace fluxcell {

/// 电化学模型类型
enum class CellModelType : uint8_t {
    SPM,  ///< 单颗粒模型（电解液均匀）
    SPMe, ///< 扩展单颗粒模型（电解液 1D 传输）
    DFN,  ///< 全阶伪二维模型（逐颗粒 + 电解液 1D + 固相欧姆）
};

/// 时间积分求解器
enum class SolverType : uint8_t {
    EigenCN,     ///< 固定步长 Crank-Nicolson（各模型内置）
    SundialsIDA, ///< SUNDIALS IDA 自适应步长（需模型提供 DAE 视图）
};

/// 工况段：恒流（CC）/ 恒压（CV）/ 恒功率（CP）充放电或静置，
/// 按时长或对应电气量截止
struct CCStep {
    enum class Mode {
        Discharge, ///< 放电（锂流出负极）
        Charge,    ///< 充电（锂流入负极）
        Rest,      ///< 静置（电流为 0，仅按时间截止）
    };

    /// 控制量类型（静置段固定为电流 0，忽略本字段）
    enum class Control {
        Current, ///< 恒流 CC：电流恒定，截止 = 时间或电压
        Voltage, ///< 恒压 CV：端电压跟踪目标，截止 = 时间或尾电流
        Power,   ///< 恒功率 CP：P = U·I 恒定，截止 = 时间或电压
    };

    Mode mode = Mode::Discharge;
    Control control = Control::Current;
    real_t current = 5.0;        ///< CC 幅值 [A]（非负；方向由 mode 决定，Rest 忽略）
    real_t voltage = 4.2;        ///< CV 目标端电压 [V]（>0；充/放方向由电压差自然决定）
    real_t power = 20.0;         ///< CP 功率 [W]（非负；方向由 mode 决定）
    real_t cutoff_time = 3600.0; ///< 段截止时间 [s]（从段开始计时；<= 0 = 不按时间截止）
    real_t cutoff_voltage = 0.0; ///< CC/CP 截止电压 [V]：放电 = 下限，充电 = 上限；
                                 ///  <= 0 = 不启用（CV 段忽略，目标电压即电压约束）
    real_t cutoff_current = 0.0; ///< CV 截止尾电流 [A]：|I| <= 该值时段结束；
                                 ///  <= 0 = 不启用
};

/// 多段工况协议
///
/// 按顺序执行各段：段的截止条件满足后进入下一段，
/// 全部段完成即仿真结束。浓度场等状态在段间连续传递（静置段电压自然弛豫）。
/// CV/CP 段的电流依赖模型状态，由引擎每步反馈计算（见 SimulationEngine::
/// segment_current），协议层仅提供控制目标与截止条件。
struct MultiStepProtocol {
    std::vector<CCStep> steps;
    /// 循环圈数：steps 序列作为一个循环重复 N 次（机理寿命仿真工况，
    /// 如 [CC 放电, 静置, CC-CV 充电] × 100 圈）；状态跨圈连续传递，
    /// 降解累积反映在每圈末的锂损失/容量诊断（快照 cycle_end 标记采样点）
    int cycles = 1;
    real_t time_step = 1.0;        ///< 积分/输出步长 [s]（IDA 模式 = 输出采样间隔）
    real_t real_time_factor = 0.0; ///< 回放速率因子：0 = 不限速，1.0 = 实时，N = N 倍速
    CellModelType model = CellModelType::SPM; ///< 电化学模型类型
    SolverType solver = SolverType::EigenCN;  ///< 时间积分求解器
    real_t solver_rel_tol = 1e-6;  ///< IDA 相对容差（EigenCN 忽略）
    model::SPMeModel::Parameters spme_params; ///< SPMe 参数（SPM 忽略；
                                              ///  solid 部分由 GUI 用运行条件覆盖）
    model::DFNModel::Parameters dfn_params;   ///< DFN 参数（SPM/SPMe 忽略；
                                              ///  solid 部分由 GUI 用运行条件覆盖）

    /// 段的标称电流：CC = 放电 +|I| / 充电 -|I|，静置 = 0，CV/CP = 0（平衡态）
    real_t step_current(const CCStep& step) const;

    /// 全部段的时间预算之和 × 循环圈数（进度估算用；
    /// 纯电压截止段按时长 0 计）
    real_t total_time_budget() const;

    /// 校验协议合法性（段非空、步长为正、控制量有效、每段至少一种截止条件），
    /// 非法时抛出
    void validate() const;
};

/// 后台仿真引擎
///
/// 线程模型：start() 派生工作线程按协议逐段推进模型（SPM/SPMe/DFN）；
/// GUI 线程每帧调用 collect_snapshots() 取走新快照。
/// 单步模式：single_step() 由主线程直接推进一步（不依赖工作线程），
/// 可与连续模式混合使用（运行中单步 = 自动暂停后续走当前进度）。
class SimulationEngine {
public:
    SimulationEngine();
    ~SimulationEngine();

    SimulationEngine(const SimulationEngine&) = delete;
    SimulationEngine& operator=(const SimulationEngine&) = delete;

    /// 启动仿真（会先终止未完成的旧仿真）
    void start(model::SPMModel::Parameters params, const MultiStepProtocol& protocol);

    /// 请求停止并等待线程退出，终止本次仿真
    void stop();

    /// 暂停 / 恢复（恢复时若工作线程不存在则自动补起，从当前进度继续连续运行）
    void pause();
    void resume();

    /// 单步执行：主线程推进一个时间步
    ///
    /// - 首次调用（无模型）时用传入参数初始化并产生 t = 0 快照；
    /// - Running 状态下调用会先停线程再步进（保留当前进度）；
    /// - 步进后状态为 Paused（段截止满足则切换下一段，全部段完成为 Completed）。
    void single_step(model::SPMModel::Parameters params, const MultiStepProtocol& protocol);

    /// 取走自上次调用以来积累的快照
    std::vector<SolutionSnapshot> collect_snapshots();

    // ── 状态查询（线程安全）──
    SimulationStatus status() const { return status_.load(std::memory_order_acquire); }
    real_t progress() const;           ///< 0 ~ 1（按全部段的时间预算）
    std::string error_message() const; ///< 仅 Failed 时有效

private:
    void run();

    /// 推进一个输出步：计算本步段电流（CV/CP 反馈）→ 积分 → 库仑计数 →
    /// 快照 → 段截止处理。全部段完成时置 Completed 并返回 true。
    ///（IDA 模式下段电流切换时自动 reinit + 代数一致化）
    bool advance_one_step();

    /// 当前段在本步应施加的电流 [A]（带符号；CP 按 I = P/U 换算），
    /// 依赖 model_ 已初始化
    real_t segment_current(const CCStep& step);

    /// CV 段推进（按模型能力分派）：
    /// - 增广 DAE（DFN）：端口电流为状态变量，外电路代数行锁端电压，
    ///   直推无外层反解（pybamm.external_circuit.ConstantVoltage 同构）
    /// - 其余（SPM/SPMe 的 IDA、CN 路径）：回退反馈修正反解
    /// 返回本步实际电流。
    real_t advance_cv(const CCStep& step);

    /// CV 段反馈修正反解（PyBaMM 式外电路求解，参考
    /// pybamm.external_circuit；增广端口不可用时的回退路径）：
    /// 电流视为代数量，每步按步首电压误差经缓存灵敏度修正后单次推进
    /// （快路径）；误差超限时回滚差分校准灵敏度重推（慢路径）。
    real_t advance_cv_step(const CCStep& step);

    /// 积分推进一个输出步：EigenCN = model->step；IDA = 自适应推进到 t + dt
    ///（IDA 模式下段电流大跳变时 reinit + 代数一致化，微调免 reinit）
    void advance_output(real_t current, real_t dt);

    /// 最近一个求解器内部步长（IDA 诊断；CN = dt）
    real_t last_internal_step() const;

    /// 生成快照入队（填充 SOC 库仑计数/循环圈数与内部步长，更新进度）
    void push_snapshot(real_t current, bool cycle_end = false);

    /// 释放模型与求解器（IDA 求解器持有模型 DAE 视图引用，必须先释放）
    void release_model();

    /// 失败时刻的模型状态摘要（拼进错误信息）：GUI 无 stderr，
    /// dae_dump_state 的详细输出不可见，状态栏错误是唯一诊断通道
    std::string failure_context() const;

    /// 停止并等待工作线程退出（保留模型与进度，供单步/恢复继续）
    void stop_thread();

    /// 初始化模型、段游标与输出队列（start 与首次 single_step 共用）
    void init_protocol_state();

    /// 当前段是否已满足截止条件（时间；CC/CP 电压；CV 尾电流）
    bool step_finished(const CCStep& step) const;

    model::SPMModel::Parameters params_;
    MultiStepProtocol protocol_;

    /// 当前模型实例：由 start()/single_step() 创建。
    /// 访问纪律——工作线程运行期间仅线程本身访问；
    /// 主线程仅在无线程时（single_step）访问，无并发。
    /// step_index_ / step_start_time_ 遵循同一纪律。
    std::unique_ptr<model::CellModel> model_;

#if defined(FLUXCELL_WITH_SUNDIALS)
    /// IDA 自适应求解器（SundialsIDA 模式；与 model_ 同一线程纪律）
    std::unique_ptr<solver::IdaSolver> ida_solver_;
#endif
    real_t ida_current_ = 0.0;       ///< IDA 已同步的段电流（大跳变/回退检测；
                                     ///  NaN = 下次推进强制 reinit）
    real_t charge_throughput_ = 0.0; ///< 累计电荷 [Ah]（SOC 库仑计数）

    /// CV 反解的试探初值：上一步段电流（段切换时继承，如 CC→CV 自然衔接）
    real_t seg_i_last_ = 0.0;
    /// CV 反馈修正的端电压-电流灵敏度 dV/dI [V/A]（< 0；一阶递推更新，
    /// 初始典型内阻 −0.02Ω 兜底）
    real_t cv_dvdi_ = -0.02;

    /// 增广端口 CV 直推的生效状态（段切换/目标变化时重新 set_voltage
    /// + reinit；CC/CP 段进入时复位）
    bool cv_augmented_active_ = false;
    real_t cv_target_applied_ = std::numeric_limits<real_t>::quiet_NaN();
                               ///< 已应用的恒压目标（NaN = 未应用）
    real_t last_current_ = 0.0; ///< 本步实际施加电流 [A]（CV 尾电流截止判据）

    std::size_t step_index_ = 0;   ///< 当前段索引（圈内）
    int cycle_index_ = 0;          ///< 当前循环圈索引（0-based）
    real_t step_start_time_ = 0.0; ///< 当前段开始时刻的模型时间 [s]

    std::thread thread_;
    std::atomic<SimulationStatus> status_{SimulationStatus::Idle};
    std::atomic<bool> pause_flag_{false};
    std::atomic<bool> stop_flag_{false};

    // 由 queue_mutex_ 保护
    mutable std::mutex queue_mutex_;
    std::deque<SolutionSnapshot> queue_;
    real_t progress_ = 0.0;
    std::string error_;
};

} // namespace fluxcell
