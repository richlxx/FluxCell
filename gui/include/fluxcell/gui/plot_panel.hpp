/// @file plot_panel.hpp
/// @brief ImPlot 实时绘图面板
#pragma once

#include <fluxcell/core/simulation/simulation.hpp>
#include <fluxcell/core/types.hpp>

#include <string>
#include <vector>

namespace fluxcell::gui {

/// 实时绘图面板：图表清单随模型类型联动（SPM/SPMe/DFN 各自的可用输出）。
/// 提供 4 路独立图表窗口（可停靠为 2×2 网格），每路可自选图表类型。
class PlotPanel {
public:
    /// 独立图表面板数量
    static constexpr int kPanelCount = 4;

    PlotPanel();
    ~PlotPanel() = default;

    /// 渲染第 index 路（0 起）图表面板窗口
    void render_panel(int index);

    /// 绑定外部快照历史（只读共享，所有权留在外部）：长循环历史若在
    /// 每帧 update 中全量深拷贝，GB 级数据的拷贝本身就是性能灾难
    void bind_history(const std::vector<SolutionSnapshot>* history);

    /// 设置当前模型类型（切换图表清单）
    void set_model_type(CellModelType type) { model_type_ = type; }

    // 图表渲染方法：注册表（cpp 匿名 namespace）需取其成员函数指针，
    // 语义上是可选择的输出视图，故公开
    // ── 全模型 ──
    void plot_voltage_curve();
    void plot_current_curve();
    void plot_voltage_breakdown();
    void plot_particle_concentration_profile();
    void plot_particle_surface_concentration();
    void plot_stoichiometry_history();
    void plot_soc_history();
    void plot_overpotential_history();

    // ── SPMe / DFN（液相）──
    void plot_concentration_profile();
    void plot_potential_profile();

    // ── DFN 分布输出 ──
    void plot_reaction_current_profile();
    void plot_overpotential_profile();
    void plot_ocp_profile();
    void plot_surface_conc_profile();
    void plot_liquid_current_profile();

    // ── DFN 热模型 ──
    void plot_temperature_history();
    void plot_heat_generation_history();

    // ── DFN 机理寿命 ──
    void plot_degradation_history();

private:
    const std::vector<SolutionSnapshot>* snapshots_ = nullptr; ///< 外部历史（只读共享）
    CellModelType model_type_ = CellModelType::SPM;
    int selected_plot_[kPanelCount] = {0, 1, 2, 3}; ///< 各面板选中的图表（全清单索引）
};

} // namespace fluxcell::gui
