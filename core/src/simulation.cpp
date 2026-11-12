/// @file simulation.cpp
/// @brief 仿真引擎实现
#include <fluxcell/core/simulation/simulation.hpp>

#include <fluxcell/core/solver/dae_problem.hpp>

#include <chrono>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fluxcell {

using model::SPMModel;

/// IDA 免重初始化的电流变化容差 [A]：段切换（数 A 级跳变）必须 reinit +
/// CalcIC 一致化代数变量；CV 反馈修正的步间微调（误差/灵敏度，典型
/// mA 级）保留 BDF 历史直接推进，Newton 修正吸收参数微扰，步内误差
/// 估计兜底（超限时 IDA 自动减步重试）。NaN 标记（时间回退）不在此列
constexpr real_t kReinitCurrentTol = 0.2;

/// CV 快路径电压锁定容差 [V]：反馈修正的稳态误差 = 步内电压动态漂移
/// （μV~百 μV 级），超限说明灵敏度失准/工况突变，转入差分校准慢路径
constexpr real_t kCvVoltTolS = 1e-3;

real_t MultiStepProtocol::step_current(const CCStep& step) const {
    switch (step.mode) {
        case CCStep::Mode::Discharge:
            return std::abs(step.current); // 正 = 放电（锂流出负极）
        case CCStep::Mode::Charge:
            return -std::abs(step.current);
        case CCStep::Mode::Rest:
            return 0.0;
    }
    return 0.0;
}

real_t MultiStepProtocol::total_time_budget() const {
    real_t total = 0.0;
    for (const auto& step : steps) {
        total += std::max(step.cutoff_time, 0.0);
    }
    // 全部为纯电压/尾电流截止段时预算为 0，用 1.0 兜底避免进度除零
    return (total > 0.0 ? total : 1.0) * static_cast<real_t>(std::max(cycles, 1));
}

void MultiStepProtocol::validate() const {
    if (steps.empty()) {
        throw std::invalid_argument("工况段为空");
    }
    if (cycles < 1) {
        throw std::invalid_argument("循环圈数必须 ≥ 1");
    }
    if (time_step <= 0.0) {
        throw std::invalid_argument("时间步长必须为正");
    }
    for (std::size_t i = 0; i < steps.size(); ++i) {
        const CCStep& step = steps[i];
        const std::string tag = "段 " + std::to_string(i + 1) + ": ";
        switch (step.control) {
            case CCStep::Control::Current:
                if (step.mode != CCStep::Mode::Rest && step.current <= 0.0) {
                    throw std::invalid_argument(tag + "恒流段电流必须为正");
                }
                break;
            case CCStep::Control::Voltage:
                if (step.mode == CCStep::Mode::Rest) {
                    throw std::invalid_argument(tag + "静置段不支持恒压控制");
                }
                if (step.voltage <= 0.0) {
                    throw std::invalid_argument(tag + "恒压段目标电压必须为正");
                }
                break;
            case CCStep::Control::Power:
                if (step.mode == CCStep::Mode::Rest) {
                    throw std::invalid_argument(tag + "静置段不支持恒功率控制");
                }
                if (step.power <= 0.0) {
                    throw std::invalid_argument(tag + "恒功率段功率必须为正");
                }
                break;
        }
        // 每段至少一种截止条件：时间 / 电压（CC/CP）/ 尾电流（CV）
        const bool has_current_cutoff = step.control == CCStep::Control::Voltage
                                     && step.cutoff_current > 0.0;
        if (step.cutoff_time <= 0.0 && step.cutoff_voltage <= 0.0
            && !has_current_cutoff) {
            throw std::invalid_argument(tag + "缺少截止条件（时间/电压/电流至少一项）");
        }
    }
}

void SimulationEngine::start(model::SPMModel::Parameters params,
                             const MultiStepProtocol& protocol) {
    stop(); // 终止旧线程并清空旧模型

    params_ = std::move(params);
    protocol_ = protocol;
    protocol_.validate();

    stop_flag_ = false;
    pause_flag_ = false;
    init_protocol_state();

    status_.store(SimulationStatus::Running, std::memory_order_release);
    thread_ = std::thread(&SimulationEngine::run, this);
}

SimulationEngine::SimulationEngine() = default;

SimulationEngine::~SimulationEngine() {
    stop();
}

/// 创建 IDA 求解器（仅 SundialsIDA 模式；模型须提供 DAE 视图）
static void setup_ida_solver(model::CellModel& model, const MultiStepProtocol& protocol,
                             real_t initial_current,
                             std::unique_ptr<solver::IdaSolver>& out,
                             real_t& ida_current) {
#if defined(FLUXCELL_WITH_SUNDIALS)
    solver::DaeProblem* problem = model.dae_problem();
    if (problem == nullptr) {
        throw std::invalid_argument("所选模型不支持 SUNDIALS IDA 求解器");
    }
    problem->dae_set_current(initial_current);
    out = std::make_unique<solver::IdaSolver>(*problem, model.time(),
                                              protocol.solver_rel_tol);
    ida_current = initial_current;
#else
    (void)model; (void)protocol; (void)initial_current;
    (void)out; (void)ida_current;
    throw std::runtime_error("本构建未启用 SUNDIALS（FLUXCELL_WITH_SUNDIALS），"
                             "无法使用 IDA 自适应求解器");
#endif
}

void SimulationEngine::init_protocol_state() {
    {
        std::lock_guard<std::mutex> lock(queue_mutex_);
        queue_.clear();
        progress_ = 0.0;
        error_.clear();
    }

    step_index_ = 0;
    cycle_index_ = 0;
    step_start_time_ = 0.0;
    charge_throughput_ = 0.0;
    seg_i_last_ = 0.0;
    last_current_ = 0.0;
    cv_augmented_active_ = false;
    cv_target_applied_ = std::numeric_limits<real_t>::quiet_NaN();
    const real_t i0 = protocol_.step_current(protocol_.steps.front());

    // 按协议选择模型：SPM 用引擎参数，SPMe/DFN 用协议携带的参数；
    // 首段 CV/CP 时 i0 = 0（平衡态起步，段内由反馈逐步建立电流）
    switch (protocol_.model) {
        case CellModelType::SPM:
            model_ = std::make_unique<SPMModel>(params_, i0);
            break;
        case CellModelType::SPMe:
            model_ = std::make_unique<model::SPMeModel>(protocol_.spme_params, i0);
            break;
        case CellModelType::DFN:
            model_ = std::make_unique<model::DFNModel>(protocol_.dfn_params, i0);
            break;
        default:
            throw std::invalid_argument("未知模型类型");
    }

    // IDA 模式：构造自适应求解器（段切换时 reinit）
    if (protocol_.solver == SolverType::SundialsIDA) {
        setup_ida_solver(*model_, protocol_, i0, ida_solver_, ida_current_);
    } else {
        ida_current_ = 0.0;
    }

    // 初始快照（t = 0）：统一经 push_snapshot 填充 SOC/诊断字段
    push_snapshot(i0);
}

void SimulationEngine::advance_output(real_t current, real_t dt) {
    if (protocol_.solver == SolverType::SundialsIDA) {
#if defined(FLUXCELL_WITH_SUNDIALS)
        cv_augmented_active_ = false; // 离开 CV 增广路径（set_current 切回
                                      // CurrentControl；不清除会使下一
                                      // 圈 CV 段误判已生效）
        solver::DaeProblem* problem = model_->dae_problem();
        problem->dae_set_current(current);
        if (!std::isfinite(ida_current_)
            || std::abs(current - ida_current_) > kReinitCurrentTol) {
            // 段切换/大跳变/时间回退（NaN 标记）：重置历史并一致化；
            // 小变化保留历史直接推进（见 kReinitCurrentTol 注释）
            ida_solver_->reinit(model_->time(), dt);
            ida_current_ = current;
        }
        ida_solver_->advance_to(model_->time() + dt);
#else
        (void)current;
        throw std::runtime_error("本构建未启用 SUNDIALS，无法使用 IDA 求解器");
#endif
    } else {
        model_->step(current, dt);
    }
}

real_t SimulationEngine::segment_current(const CCStep& step) {
    const real_t v_now = model_->voltage();
    switch (step.control) {
        case CCStep::Control::Current:
            return protocol_.step_current(step); // CC / 静置
        case CCStep::Control::Voltage:
            return seg_i_last_; // 由 advance_cv_step 试探反解，此值仅为初值
        case CCStep::Control::Power: {
            // I = P/U：用步初端电压换算（步内 U 变化带来的功率偏差为步级误差）；
            // 端电压过低时钳位防电流爆炸
            const real_t v_safe = std::max(std::abs(v_now), 0.5);
            const real_t p = std::abs(step.power);
            return step.mode == CCStep::Mode::Charge ? -p / v_safe : p / v_safe;
        }
    }
    return 0.0;
}

real_t SimulationEngine::advance_cv(const CCStep& step) {
#if defined(FLUXCELL_WITH_SUNDIALS)
    // 增广 DAE（DFN）：pybamm.external_circuit.ConstantVoltage 同构——
    // 端口电流是状态中的代数变量，外电路行锁端电压目标，Newton 与
    // 内部方程（BV/守恒/热）一致求解，直推无外层反解
    if (protocol_.solver == SolverType::SundialsIDA && ida_solver_) {
        solver::DaeProblem* prob = model_->dae_problem();
        if (prob != nullptr && prob->dae_port_augmented()) {
            if (!cv_augmented_active_
                || cv_target_applied_ != step.voltage) {
                // 段切换/目标变化：外电路行表达式突变（残差不连续），
                // 重置历史并一致化（端口电流继承 CC 段末值，端电压
                // 已 ≈ 目标，CalcIC 最小修正收敛）
                prob->dae_set_voltage(step.voltage);
                cv_augmented_active_ = true;
                cv_target_applied_ = step.voltage;
                ida_current_ = std::numeric_limits<real_t>::quiet_NaN();
                ida_solver_->reinit(model_->time(), protocol_.time_step);
            }
            ida_solver_->advance_to(model_->time() + protocol_.time_step);
            return prob->dae_port_current();
        }
    }
#endif
    // 回退：反馈修正反解（SPM/SPMe 的 IDA、CN 路径）
    return advance_cv_step(step);
}

real_t SimulationEngine::advance_cv_step(const CCStep& step) {
    // PyBaMM 式外电路求解（参考 pybamm.external_circuit.ConstantVoltage）：
    // 电流视为代数量，端电压跟踪 V_target。
    // 反馈修正范式：每步按步首电压误差经缓存灵敏度 dV/dI 修正电流后
    // 单次推进（快路径，稳态误差 = 步内电压动态漂移，mV 以下）；误差
    // 超限才回滚差分校准灵敏度重推（慢路径）。原「每步两次试探 + 回滚
    // 重推」范式的时间倒退必然 IDAReInit+CalcIC（每步 2 次）是 CV 段
    // 慢的主因；快路径配合小电流微调免 reinit（kReinitCurrentTol），
    // CV 每步成本与 CC 段同级。逐步割线法在 SOC 平台区灵敏度失配时
    // 会振荡的问题由一阶平滑 + 慢路径实时校准抑制。
    solver::DaeProblem* prob = model_->dae_problem();
    Vector y0;
    Vector yp0;
    prob->dae_pack_state(y0, yp0); // 步首状态（慢路径回滚用）
    const real_t t0 = model_->time();
    const real_t dt = protocol_.time_step;
    const real_t v_prev = model_->voltage();
    const real_t i_prev = seg_i_last_;

    // 快路径：反馈修正电流，单次推进
    const real_t lim = std::max(std::abs(i_prev) * 2.0, 10.0);
    real_t i_try = i_prev;
    if (std::abs(v_prev - step.voltage) > 1e-6) {
        i_try = std::clamp(i_prev - (v_prev - step.voltage) / cv_dvdi_,
                           -lim, lim);
    }
    advance_output(i_try, dt);
    const real_t v1 = model_->voltage();
    const real_t err = v1 - step.voltage;
    seg_i_last_ = i_try;

    // 灵敏度递推更新（ΔV/ΔI，一阶平滑；分母过小或比值异常时保持）
    const real_t di_step = i_try - i_prev;
    if (std::abs(di_step) > 1e-3) {
        const real_t dvdi = (v1 - v_prev) / di_step;
        if (dvdi < -1e-4 && dvdi > -10.0) {
            cv_dvdi_ = 0.5 * (cv_dvdi_ + dvdi);
        }
    }
    if (std::abs(err) <= kCvVoltTolS) {
        return i_try; // 电压锁定达标，本步完成
    }

    // 慢路径：差分校准灵敏度（回退试探）→ Newton 反解 → 回滚重推
    const real_t di = (err > 0.0 ? 1.0 : -1.0)
                    * std::max(0.1, 0.05 * std::abs(i_try));
    prob->dae_load_state(t0, y0);
#if defined(FLUXCELL_WITH_SUNDIALS)
    // IDA 历史/时间已被快路径推进污染，强制重建（NaN != 任何电流）
    ida_current_ = std::numeric_limits<real_t>::quiet_NaN();
#endif
    advance_output(i_try + di, dt);
    const real_t v2 = model_->voltage();
    const real_t dvdi2 = (v2 - v1) / di;
    if (dvdi2 < -1e-6) {
        cv_dvdi_ = 0.5 * (cv_dvdi_ + dvdi2);
    }

    real_t i2 = i_try - err / cv_dvdi_;
    i2 = std::clamp(i2, -lim, lim);

    // 回滚到步首按校准电流正式推进（时间回退强制重建）
    prob->dae_load_state(t0, y0);
#if defined(FLUXCELL_WITH_SUNDIALS)
    ida_current_ = std::numeric_limits<real_t>::quiet_NaN();
#endif
    advance_output(i2, dt);
    seg_i_last_ = i2;
    return i2;
}

bool SimulationEngine::advance_one_step() {
    const CCStep& step = protocol_.steps[step_index_];
    real_t current;
    if (step.control == CCStep::Control::Voltage) {
        // CV：增广端口直推或反馈反解（内部完成最终推进，模型已处于本步末状态）
        current = advance_cv(step);
    } else {
        current = segment_current(step);
        advance_output(current, protocol_.time_step);
    }
    // 带符号库仑计数：放电（+I）累积、充电（−I）抵扣、静置不变
    charge_throughput_ += current * protocol_.time_step / 3600.0;
    last_current_ = current;
    seg_i_last_ = current;

    // 段截止：切换下一段/下一圈（快照标记圈末，供寿命曲线采样）
    const bool finished = step_finished(step);
    const bool last_step_of_cycle = step_index_ + 1 >= protocol_.steps.size();
    push_snapshot(current, finished && last_step_of_cycle);

    if (finished) {
        ++step_index_;
        step_start_time_ = model_->time();
        if (step_index_ >= protocol_.steps.size()) {
            step_index_ = 0;
            ++cycle_index_;
            if (cycle_index_ >= protocol_.cycles) {
                {
                    std::lock_guard<std::mutex> lock(queue_mutex_);
                    progress_ = 1.0;
                }
                release_model();
                status_.store(SimulationStatus::Completed, std::memory_order_release);
                return true;
            }
        }
    }
    return false;
}

real_t SimulationEngine::last_internal_step() const {
#if defined(FLUXCELL_WITH_SUNDIALS)
    if (protocol_.solver == SolverType::SundialsIDA && ida_solver_) {
        return ida_solver_->last_step();
    }
#else
    (void)this;
#endif
    return protocol_.time_step;
}

void SimulationEngine::push_snapshot(real_t current, bool cycle_end) {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    SolutionSnapshot snap = model_->snapshot(current);
    snap.internal_step = last_internal_step();
    snap.cycle_number = cycle_index_ + 1; // 1-based 圈数
    snap.cycle_end = cycle_end;
    // SOC 库仑计数：soc = initial_soc − throughput/cap（放电 +I 下降、充电 −I 回升）；
    // 容量无效时保持默认 0
    const real_t cap = model_->nominal_capacity();
    if (cap > 0.0) {
        snap.soc = model_->initial_soc() - charge_throughput_ / cap;
    }
    queue_.push_back(std::move(snap));
    progress_ = std::clamp(model_->time() / protocol_.total_time_budget(),
                           0.0, 1.0);
}

void SimulationEngine::release_model() {
#if defined(FLUXCELL_WITH_SUNDIALS)
    ida_solver_.reset(); // 先于 model_：求解器持有模型 DAE 视图引用
#endif
    model_.reset();
}

std::string SimulationEngine::failure_context() const {
    if (model_ == nullptr) {
        return {};
    }
    // 求解器 throw 已带 t/flag/步长，这里补电芯级状态（CV 段崩溃定位用）
    return " [失败时刻 t=" + std::to_string(model_->time())
         + " s, V=" + std::to_string(model_->voltage())
         + " V, I=" + std::to_string(last_current_) + " A]";
}

void SimulationEngine::stop_thread() {
    if (thread_.joinable()) {
        stop_flag_ = true;
        pause_flag_ = false;
        thread_.join();
    }
}

void SimulationEngine::stop() {
    stop_thread();
    release_model(); // 终止本次仿真，下次 start/single_step 重新初始化
    status_.store(SimulationStatus::Idle, std::memory_order_release);
}

void SimulationEngine::pause() {
    // 仅连续运行中可暂停：置暂停标志并同步状态，
    // 使 GUI 状态栏/按钮与线程实际行为一致
    if (status_.load(std::memory_order_acquire) == SimulationStatus::Running) {
        pause_flag_ = true;
        status_.store(SimulationStatus::Paused, std::memory_order_release);
    }
}

void SimulationEngine::resume() {
    if (status_.load(std::memory_order_acquire) != SimulationStatus::Paused) {
        pause_flag_ = false;
        return;
    }

    if (!thread_.joinable()) {
        // 纯单步产生的暂停：补起工作线程，从当前进度继续连续运行
        stop_flag_ = false;
        pause_flag_ = false;
        status_.store(SimulationStatus::Running, std::memory_order_release);
        thread_ = std::thread(&SimulationEngine::run, this);
    } else {
        // 唤醒 start() 后被暂停的线程，恢复运行状态
        pause_flag_ = false;
        status_.store(SimulationStatus::Running, std::memory_order_release);
    }
}

void SimulationEngine::single_step(model::SPMModel::Parameters params,
                                   const MultiStepProtocol& protocol) {
    // 连续运行中：先停线程（保留模型与进度），转为单步接管
    if (status_.load(std::memory_order_acquire) == SimulationStatus::Running) {
        stop_thread();
        status_.store(SimulationStatus::Paused, std::memory_order_release);
    }
    if (status_.load(std::memory_order_acquire) != SimulationStatus::Idle
        && status_.load(std::memory_order_acquire) != SimulationStatus::Paused) {
        return; // Completed / Failed 无步可走
    }

    try {
        // 首次单步：初始化模型、段游标与输出队列
        if (!model_) {
            params_ = std::move(params);
            protocol_ = protocol;
            protocol_.validate();
            stop_flag_ = false;
            pause_flag_ = false;
            init_protocol_state();
        }

        // 推进一步；段截止/全部完成在 advance_one_step 内处理
        if (!advance_one_step()) {
            status_.store(SimulationStatus::Paused, std::memory_order_release);
        }
    } catch (const std::exception& e) {
        const std::string ctx = failure_context();
        release_model();
        std::lock_guard<std::mutex> lock(queue_mutex_);
        error_ = e.what() + ctx;
        status_.store(SimulationStatus::Failed, std::memory_order_release);
    }
}

std::vector<SolutionSnapshot> SimulationEngine::collect_snapshots() {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    // 快照含 14 个分布向量（DFN ~22 KB），移动而非拷贝转出队列
    std::vector<SolutionSnapshot> out;
    out.reserve(queue_.size());
    std::move(queue_.begin(), queue_.end(), std::back_inserter(out));
    queue_.clear();
    return out;
}

real_t SimulationEngine::progress() const {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    return progress_;
}

std::string SimulationEngine::error_message() const {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    return error_;
}

bool SimulationEngine::step_finished(const CCStep& step) const {
    // 段内计时（1e-9 容差吸收浮点累计误差）
    if (step.cutoff_time > 0.0
        && model_->time() - step_start_time_ + 1e-9 >= step.cutoff_time) {
        return true;
    }
    // 电压截止（CC/CP）：放电看下限，充电看上限；静置无方向语义，不检查
    const real_t v = model_->voltage();
    if (step.mode == CCStep::Mode::Discharge && step.cutoff_voltage > 0.0
        && v <= step.cutoff_voltage) {
        return true;
    }
    if (step.mode == CCStep::Mode::Charge && step.cutoff_voltage > 0.0
        && v >= step.cutoff_voltage) {
        return true;
    }
    // 尾电流截止（CV）：|I| 衰减到阈值（CC-CV 充电末段惯例）
    if (step.control == CCStep::Control::Voltage && step.cutoff_current > 0.0
        && std::abs(last_current_) <= step.cutoff_current) {
        return true;
    }
    return false;
}

void SimulationEngine::run() {
    try {
        // 限速基准：把时钟起点前移已仿真时间，
        // 保证 resume 补起线程时（model_->time() > 0）回放节奏正确
        auto t_start = std::chrono::steady_clock::now()
                     - std::chrono::duration<double>(
                             protocol_.real_time_factor > 0.0
                                 ? model_->time() / protocol_.real_time_factor
                                 : 0.0);

        while (!stop_flag_) {
            // 暂停等待
            const auto pause_enter = std::chrono::steady_clock::now();
            while (pause_flag_ && !stop_flag_) {
                std::this_thread::sleep_for(std::chrono::milliseconds(20));
            }
            if (stop_flag_) {
                break;
            }
            // 暂停时长不计入回放时钟：基准整体后移，
            // 避免恢复后 sleep_until 因基准过期而全速追帧
            t_start += std::chrono::steady_clock::now() - pause_enter;

            // 回放限速：仿真时间 / rtf = 应消耗的墙钟时间
            if (protocol_.real_time_factor > 0.0) {
                const auto scheduled =
                    t_start + std::chrono::duration<double>(model_->time() / protocol_.real_time_factor);
                std::this_thread::sleep_until(scheduled);
            }

            if (advance_one_step()) {
                return; // 全部段完成（Completed 已在内部置位）
            }
        }
    } catch (const std::exception& e) {
        const std::string ctx = failure_context();
        release_model();
        std::lock_guard<std::mutex> lock(queue_mutex_);
        error_ = e.what() + ctx;
        status_.store(SimulationStatus::Failed, std::memory_order_release);
        return;
    }

    // stop_thread() 被调用：线程自然退出，模型与进度保留
}

} // namespace fluxcell
