/// @file main_window.hpp
/// @brief 主窗口布局（停靠空间、菜单栏、仿真集成）
#pragma once

#include <fluxcell/core/model/spm_model.hpp>
#include <fluxcell/core/simulation/simulation.hpp>
#include <fluxcell/core/types.hpp>

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace fluxcell::gui {

class PlotPanel;
class SimulationPanel;

/// 主窗口：管理所有面板和仿真状态
class MainWindow {
public:
    MainWindow();
    ~MainWindow();

    /// 渲染一帧（在主循环中调用）
    void render();

    /// 注册退出回调（菜单“退出”触发；Application 注入 request_exit）
    void set_exit_callback(std::function<void()> cb) { on_exit_ = std::move(cb); }

private:
    void render_menu_bar();
    void render_toolbar();
    void render_view_menu();
    void render_dockspace();
    void render_status_bar();

    /// 打开 CSV 导出保存对话框（菜单“导出结果 (CSV)…”共用）
    void open_export_csv_dialog();
    /// 每帧轮询导出对话框，选中路径后导出 snapshots_ 并弹窗反馈
    void render_export_csv_dialog();
    /// 结果反馈弹窗（导出成功/失败提示，模态）
    void render_notify_popup();
    /// “关于 FluxCell” 弹窗（模态）
    void render_about_popup();

    /// 从引擎取走新快照并刷新绘图/状态
    void pump_simulation();

    // 面板
    std::unique_ptr<PlotPanel> plot_panel_;
    std::unique_ptr<SimulationPanel> simulation_panel_;

    // 仿真
    SimulationEngine engine_;
    model::SPMModel::Parameters spm_params_; ///< 默认 Chen2020，可经面板加载

    // 仿真数据
    std::vector<SolutionSnapshot> snapshots_;
    SimulationStatus sim_status_ = SimulationStatus::Idle;

    // 面板可见性（查看菜单控制；图表面板数量与 PlotPanel::kPanelCount 对应）
    bool show_simulation_panel_ = true;
    bool show_plot_panels_[4] = {true, true, true, true};
    bool show_status_bar_ = true;
    bool show_toolbar_ = true;
    bool rebuild_layout_ = false; ///< 查看菜单请求的默认停靠布局重建

    // 菜单功能：退出 / 关于 / 导出反馈
    std::function<void()> on_exit_; ///< 菜单“退出”回调
    bool show_about_ = false;       ///< 请求打开“关于”弹窗（下一帧消费）
    std::string notify_message_;    ///< 结果反馈弹窗内容（空 = 无）
    bool notify_pending_ = false;   ///< 下一帧打开结果反馈弹窗
};

} // namespace fluxcell::gui
