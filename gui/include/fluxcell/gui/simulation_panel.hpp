/// @file simulation_panel.hpp
/// @brief 仿真控制面板：多段工况编辑 + Run/Pause/Stop 控制 + 进度显示
#pragma once

#include <fluxcell/core/model/spm_model.hpp>
#include <fluxcell/core/simulation/simulation.hpp>
#include <fluxcell/core/types.hpp>

#include <functional>
#include <string>
#include <vector>

namespace fluxcell::gui {

/// 单段工况的 UI 缓存
struct StepUI {
    int mode = 0;                 ///< 0 = 放电，1 = 充电，2 = 静置
    int control = 0;              ///< 0 = 恒流 CC，1 = 恒压 CV，2 = 恒功率 CP
                                  ///  （静置段忽略，恒为电流 0）
    float current = 5.0f;         ///< 恒流幅值 [A]（静置时忽略）
    float voltage = 4.2f;         ///< 恒压目标端电压 [V]
    float power = 20.0f;          ///< 恒功率 [W]
    float cutoff_time = 3600.0f;  ///< 段截止时间 [s]；<= 0 = 不按时间截止
    float cutoff_voltage = 2.5f;  ///< CC/CP 段截止电压 [V]；<= 0 = 不启用
    float cutoff_current = 0.05f; ///< CV 段截止尾电流 [A]；<= 0 = 不启用
    bool cccv = false;            ///< CC 段：电压截止后转恒压（CC-CV 一键工况）
};

/// 仿真控制面板
///
/// 持有 SimulationEngine 与 SPM 参数的引用，编辑多段工况后启动后台仿真。
/// 参数区随模型选择联动（SPM/SPMe/DFN 各自的可调参数与求解器选项）。
class SimulationPanel {
public:
    SimulationPanel(SimulationEngine& engine, model::SPMModel::Parameters& params);

    /// 渲染面板（每帧调用）
    void render();

    /// 供菜单栏快捷操作
    void start();
    void pause_resume();
    void stop();
    void reset();
    void step_once();

    /// 打开参数文件选择对话框（面板“加载…”按钮 / 菜单“加载参数…”/ Ctrl+O 共用）
    void open_file_dialog();

    /// 打开参数保存对话框（菜单“保存参数…”/ Ctrl+S 共用）；
    /// 保存按 UI 当前模型类型写出对应层级的参数（文件口径，不含面板运行条件）
    void open_save_file_dialog();

    /// 打开工况文件选择对话框（面板“载入工况…”按钮 / 菜单“载入工况…”共用）
    void open_protocol_dialog();

    /// 打开工况保存对话框（面板“保存工况…”按钮 / 菜单“保存工况…”共用）：
    /// 写出当前段序列 + 循环圈数 + 模型/求解器选择 + 初始 SOC/环境温度
    void open_save_protocol_dialog();

    /// 注册重启回调：start()/reset() 时触发，用于清空旧仿真曲线
    void set_restart_callback(std::function<void()> cb);

    /// 当前选中的模型类型（图表/输出联动）
    CellModelType current_model() const { return static_cast<CellModelType>(model_type_); }

private:
    void render_protocol_section();
    void render_model_params_section();
    void render_cell_params();
    void render_progress_section();

    /// 从 params_path_ 加载 YAML 参数（UI 手输回车 / 文件对话框选中后共用）
    void load_params();

    /// 将 UI 当前模型口径的参数保存到 params_path_（保存对话框选中后共用）
    void save_params();

    /// 从 protocol_path_ 加载工况并覆盖 UI 缓存（手输回车 / 文件对话框选中后共用）
    void load_protocol();

    /// 将 UI 当前工况保存到 protocol_path_（保存对话框选中后共用）
    void save_protocol();

    /// 从 UI 缓存构造多段工况协议（start 与单步共用）
    MultiStepProtocol build_protocol() const;

    /// 参数文件加载值 + 仿真条件（初始 SOC / 环境温度）合并后的模型参数
    model::SPMModel::Parameters effective_params() const;

    SimulationEngine& engine_;
    model::SPMModel::Parameters& params_;
    std::function<void()> on_restart_; ///< 重新开始仿真前的清理回调

    // ── 仿真条件（UI 缓存）──
    float initial_soc_pct_ = 100.0f; ///< 初始 SOC [%]
    float ambient_temp_c_ = 25.0f;   ///< 环境温度 [℃]（等温假设）

    // ── 模型选择（UI 缓存）──
    int model_type_ = 2;             ///< 0 = SPM，1 = SPMe，2 = DFN（默认 DFN）
    int thermal_type_ = 0;           ///< 热模型：0 = 等温，1 = 集总热（仅 DFN）
    model::SPMeModel::Parameters spme_params_; ///< SPMe 电解液参数缓存（YAML 加载）
    model::DFNModel::Parameters dfn_params_;   ///< DFN 固相电导等参数缓存（YAML 加载）

    // ── 求解器选择（UI 缓存）──
    int solver_type_ = 1;            ///< 0 = CN 固定步长，1 = SUNDIALS IDA 自适应（默认 IDA）
    float solver_rel_tol_ = 1e-6f;   ///< IDA 相对容差

    // ── 工况设置（UI 缓存）──
    std::vector<StepUI> steps_;
    int cycles_ = 1;                 ///< 循环圈数（段序列整体重复 N 次）
    float time_step_ = 1.0f;         ///< 积分步长 [s]
    float real_time_factor_ = 0.0f;  ///< 回放速率（0 = 不限速）

    char params_path_[512] = "params/dfn_chen2020.yml";
    std::string load_error_;
    std::string save_error_; ///< 参数保存失败信息（成功时清空）
    std::string protocol_error_; ///< 协议校验错误提示（如缺少截止条件）
    std::string error_select_buf_; ///< 失败信息只读多行框缓冲（std::string 扩容配方）

    // ── 工况文件（多段协议 + 运行条件的存取）──
    char protocol_path_[512] = ""; ///< 工况文件路径（空 = 保存/载入对话框落到 params/）
    std::string protocol_load_error_;  ///< 工况载入失败信息
    std::string protocol_save_error_;  ///< 工况保存失败信息（成功时清空）
};

} // namespace fluxcell::gui
