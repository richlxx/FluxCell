/// @file plot_panel.cpp
/// @brief ImPlot 实时绘图面板实现
#include <fluxcell/gui/plot_panel.hpp>

#include <fluxcell/gui/plot_util.hpp>

#include <imgui.h>
#include <implot.h>

#include <cstddef>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace fluxcell::gui {

namespace {

/// 图表可见性掩码位（与 CellModelType 枚举值对应）
constexpr unsigned kSPM = 1u << 0;
constexpr unsigned kSPMe = 1u << 1;
constexpr unsigned kDFN = 1u << 2;
constexpr unsigned kAll = kSPM | kSPMe | kDFN;
constexpr unsigned kLiquid = kSPMe | kDFN; ///< 有液相 1D 分布的模型

/// 图表定义：标签 + 适用模型掩码 + 渲染函数
struct PlotSpec {
    const char* label;
    unsigned model_mask;
    void (PlotPanel::*render_fn)();
};

const PlotSpec kPlots[] = {
    {"电压曲线", kAll, &PlotPanel::plot_voltage_curve},
    {"电流曲线", kAll, &PlotPanel::plot_current_curve},
    {"电压分解", kAll, &PlotPanel::plot_voltage_breakdown},
    {"颗粒径向浓度分布", kAll, &PlotPanel::plot_particle_concentration_profile},
    {"颗粒表面浓度（时序）", kAll, &PlotPanel::plot_particle_surface_concentration},
    {"嵌入率 x/y（时序）", kAll, &PlotPanel::plot_stoichiometry_history},
    {"SOC（时序）", kAll, &PlotPanel::plot_soc_history},
    {"过电位 η（时序）", kAll, &PlotPanel::plot_overpotential_history},
    {"电解液浓度分布", kLiquid, &PlotPanel::plot_concentration_profile},
    {"电解液/固相电位分布", kLiquid, &PlotPanel::plot_potential_profile},
    {"反应电流密度 i_surf(x)", kDFN, &PlotPanel::plot_reaction_current_profile},
    {"过电位分布 η(x)", kDFN, &PlotPanel::plot_overpotential_profile},
    {"开路电位 U(x)", kDFN, &PlotPanel::plot_ocp_profile},
    {"表面浓度沿极片 c_surf(x)", kDFN, &PlotPanel::plot_surface_conc_profile},
    {"液相面电流 i_e(x)", kDFN, &PlotPanel::plot_liquid_current_profile},
    {"电芯温度（时序）", kDFN, &PlotPanel::plot_temperature_history},
    {"产热功率分解（时序）", kDFN, &PlotPanel::plot_heat_generation_history},
    {"寿命曲线（降解）", kDFN, &PlotPanel::plot_degradation_history},
};

/// 时序图表滑动窗口长度：只显示最后 5 小时的数据点。长循环仿真的
/// 历史无限增长，全量绘制会把近期细节压缩成密集竖线；空间分布图
/// 与寿命曲线（x = 圈数）不受此窗口影响
constexpr real_t kHistoryWindowS = 18000.0;

/// 时序窗口起始索引：快照 time 单调不减，lower_bound 定位首个
/// time >= 末时刻 − 窗口 的快照；回退一个采样点使曲线从窗口左缘
/// 自然延伸入内（不足窗口长度时覆盖全部数据）
int history_window_begin(const std::vector<SolutionSnapshot>& snaps) {
    const real_t t_min = snaps.back().time - kHistoryWindowS;
    auto it = std::lower_bound(
        snaps.begin(), snaps.end(), t_min,
        [](const SolutionSnapshot& s, real_t v) { return s.time < v; });
    int begin = static_cast<int>(it - snaps.begin());
    if (begin > 0) {
        --begin; // 左缘锚点
    }
    return begin;
}

/// 电极内控制体中心无量纲位置（0 = 集流体侧，1 = 隔膜侧）
std::vector<real_t> electrode_x(int n) {
    std::vector<real_t> x(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        x[static_cast<std::size_t>(i)] =
            (static_cast<real_t>(i) + 0.5) / static_cast<real_t>(n);
    }
    return x;
}

} // namespace

PlotPanel::PlotPanel() {
    ImPlot::CreateContext();
    // 显式跟随当前 ImGui 主题（亮色）——避免依赖 CreateContext 时的隐式感知顺序
    ImPlot::StyleColorsAuto();
}

void PlotPanel::render_panel(int index) {
    char title[32];
    std::snprintf(title, sizeof(title), "图表 %d", index + 1);
    ImGui::Begin(title);

    // 图表类型选择：仅列出当前模型支持的图表；
    // 模型切换后若当前选中不可用则回落到第一个可用项。
    // 各面板独立记忆选择，互不影响。
    const unsigned mask = 1u << static_cast<unsigned>(model_type_);
    std::vector<const char*> labels;
    std::vector<int> indices;
    int selected_pos = -1;
    for (int i = 0; i < static_cast<int>(IM_ARRAYSIZE(kPlots)); ++i) {
        if ((kPlots[i].model_mask & mask) == 0) {
            continue;
        }
        if (i == selected_plot_[index]) {
            selected_pos = static_cast<int>(labels.size());
        }
        labels.push_back(kPlots[i].label);
        indices.push_back(i);
    }
    if (selected_pos < 0) {
        selected_pos = 0;
    }
    ImGui::TextUnformatted("图表类型");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::Combo("##plot_type", &selected_pos,
                 selected_pos < static_cast<int>(labels.size()) ? labels.data()
                                                                : nullptr,
                 static_cast<int>(labels.size()));
    selected_plot_[index] = indices[static_cast<std::size_t>(selected_pos)];

    ImGui::Separator();

    (this->*kPlots[selected_plot_[index]].render_fn)();

    ImGui::End();
}

void PlotPanel::plot_voltage_curve() {
    if (snapshots_ == nullptr || snapshots_->empty()) {
        ImGui::TextDisabled("暂无仿真数据");
        return;
    }

    const int begin = history_window_begin(*snapshots_);
    const std::size_t count = snapshots_->size() - static_cast<std::size_t>(begin);
    std::vector<real_t> times, voltages;
    times.reserve(count);
    voltages.reserve(count);

    for (int i = begin; i < static_cast<int>(snapshots_->size()); ++i) {
        const auto& snap = (*snapshots_)[i];
        times.push_back(snap.time);
        voltages.push_back(snap.voltage);
    }

    if (ImPlot::BeginPlot("电压-时间曲线", ImVec2(-1, -1))) {
        // AutoFit：两轴随数据自适应缩放（含边距）
        ImPlot::SetupAxes("时间 [s]", "电压 [V]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(voltages);
        ImPlot::SetupFinish();
        ImPlot::PlotLine("电压", times.data(), voltages.data(), (int)times.size());
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_current_curve() {
    if (snapshots_ == nullptr || snapshots_->empty()) {
        ImGui::TextDisabled("暂无仿真数据");
        return;
    }

    const int begin = history_window_begin(*snapshots_);
    const std::size_t count = snapshots_->size() - static_cast<std::size_t>(begin);
    std::vector<real_t> times, currents;
    times.reserve(count);
    currents.reserve(count);

    for (int i = begin; i < static_cast<int>(snapshots_->size()); ++i) {
        const auto& snap = (*snapshots_)[i];
        times.push_back(snap.time);
        currents.push_back(snap.current); // 带符号：正 = 放电，负 = 充电
    }

    if (ImPlot::BeginPlot("电流-时间曲线", ImVec2(-1, -1))) {
        ImPlot::SetupAxes("时间 [s]", "电流 [A]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(currents);
        ImPlot::SetupFinish();
        ImPlot::PlotLine("电流", times.data(), currents.data(),
                         (int)times.size());
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_voltage_breakdown() {
    if (snapshots_ == nullptr || snapshots_->empty()) {
        ImGui::TextDisabled("暂无仿真数据");
        return;
    }

    // 堆叠分解：V = U_ocv + η_p − η_n（过电压 = η_n − η_p）。
    // 放电时 η_n > 0、η_p < 0 → 曲线逐层下移，端电压低于开路；
    // 充电时符号反转，层叠向上。数学关系对任意电流方向均成立。
    const int begin = history_window_begin(*snapshots_);
    const std::size_t count = snapshots_->size() - static_cast<std::size_t>(begin);
    std::vector<real_t> times, ocv, mid, volts;
    times.reserve(count);
    ocv.reserve(count);
    mid.reserve(count);
    volts.reserve(count);

    for (int i = begin; i < static_cast<int>(snapshots_->size()); ++i) {
        const auto& snap = (*snapshots_)[i];
        const real_t u = snap.u_pos - snap.u_neg;
        times.push_back(snap.time);
        ocv.push_back(u);
        mid.push_back(u - snap.eta_neg);   // 扣除负极活化损失
        volts.push_back(snap.voltage);     // = mid + eta_pos（端电压）
    }
    const int n = static_cast<int>(times.size());

    if (ImPlot::BeginPlot("电压堆叠分解（V = OCP + η_p − η_n）", ImVec2(-1, -1))) {
        ImPlot::SetupAxes("时间 [s]", "电压 [V]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(ocv, mid, volts);
        ImPlot::SetupFinish();
        ImPlot::PlotShaded("负极活化损失 η_n", times.data(), ocv.data(), mid.data(), n);
        ImPlot::PlotShaded("正极活化损失 η_p", times.data(), mid.data(), volts.data(), n);
        // 顶/底边界线：开路电压与端电压
        ImPlot::PlotLine("开路电压 U_ocv", times.data(), ocv.data(), n);
        ImPlot::SetNextLineStyle(IMPLOT_AUTO_COL, 2.5f);
        ImPlot::PlotLine("端电压 V", times.data(), volts.data(), n);
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_particle_concentration_profile() {
    if (snapshots_->empty() || snapshots_->back().c_s_neg.size() == 0) {
        ImGui::TextDisabled("暂无颗粒浓度数据");
        return;
    }

    const auto& snap = snapshots_->back();
    const int n = static_cast<int>(snap.c_s_neg.size());

    // 控制体中心无量纲半径 r/R = (i + 0.5)/n（0 = 颗粒中心，1 = 表面）
    std::vector<real_t> x(n);
    for (int i = 0; i < n; ++i) {
        x[i] = (static_cast<real_t>(i) + 0.5) / static_cast<real_t>(n);
    }

    if (ImPlot::BeginPlot("颗粒内径向浓度分布（最新时刻）", ImVec2(-1, -1))) {
        ImPlot::SetupAxes("无量纲半径 r/R", "c_s [mol/m^3]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(snap.c_s_neg, snap.c_s_pos);
        ImPlot::SetupFinish();
        ImPlot::PlotLine("负极（石墨）", x.data(), snap.c_s_neg.data(), n);
        ImPlot::PlotLine("正极（NMC）", x.data(), snap.c_s_pos.data(), n);
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_concentration_profile() {
    if (snapshots_->empty() || snapshots_->back().c_e.size() == 0) {
        ImGui::TextDisabled("暂无浓度数据（SPM 假设电解液均匀，请切换到 SPMe 模型）");
        return;
    }

    const auto& snap = snapshots_->back();
    const int n = static_cast<int>(snap.c_e.size());

    // 均匀网格控制体中心的无量纲位置 x/L
    std::vector<real_t> x(n);
    for (int i = 0; i < n; ++i) {
        x[i] = (static_cast<real_t>(i) + 0.5) / static_cast<real_t>(n);
    }

    if (ImPlot::BeginPlot("电解液浓度（最新时刻）", ImVec2(-1, 300))) {
        ImPlot::SetupAxes("位置 [x/L]", "c_e [mol/m^3]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(snap.c_e);
        ImPlot::SetupFinish();
        ImPlot::PlotLine("c_e", x.data(), snap.c_e.data(), n);
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_potential_profile() {
    if (snapshots_->empty() || snapshots_->back().phi_e.size() == 0) {
        ImGui::TextDisabled("暂无电位数据（SPM 不建模液相电位，请切换到 SPMe/DFN 模型）");
        return;
    }

    const auto& snap = snapshots_->back();
    const int n = static_cast<int>(snap.phi_e.size());

    // 均匀网格控制体中心的无量纲位置 x/L
    std::vector<real_t> x(n);
    for (int i = 0; i < n; ++i) {
        x[i] = (static_cast<real_t>(i) + 0.5) / static_cast<real_t>(n);
    }

    if (ImPlot::BeginPlot("电位分布（最新时刻）", ImVec2(-1, 300))) {
        ImPlot::SetupAxes("位置 [x/L]", "phi [V]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(snap.phi_e, snap.phi_s);
        ImPlot::SetupFinish();
        // 固相电位：负极/正极是两个独立导体，隔膜区无固相——分段绘制，
        // 曲线在隔膜区断开（快照缺区间信息时退化为整条绘制）。
        // 注意两段电位绝对值差异大（负极 ≈ U_n+η ~ 0.1 V，正极 ≈ 4 V，
        // 参考 φ_e(0)=0），端电压 = 正极集流体端 − 负极集流体端固相电位差
        if (static_cast<int>(snap.phi_s.size()) == n) {
            if (snap.sep_start >= 0 && snap.pos_start > snap.sep_start
                && snap.pos_start < n) {
                ImPlot::PlotLine("phi_s（负极）", x.data(),
                                 snap.phi_s.data(), snap.sep_start);
                ImPlot::PlotLine("phi_s（正极）", x.data() + snap.pos_start,
                                 snap.phi_s.data() + snap.pos_start,
                                 n - snap.pos_start);
            } else {
                ImPlot::PlotLine("phi_s", x.data(), snap.phi_s.data(), n);
            }
        }
        ImPlot::PlotLine("phi_e", x.data(), snap.phi_e.data(), n);
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_particle_surface_concentration() {
    if (snapshots_->empty() || snapshots_->back().c_s_surf_neg.size() == 0) {
        ImGui::TextDisabled("暂无颗粒表面浓度数据");
        return;
    }

    // SPM 每电极仅一个代表颗粒：表面浓度为单点，按时序展示才有意义
    const int begin = history_window_begin(*snapshots_);
    const std::size_t count = snapshots_->size() - static_cast<std::size_t>(begin);
    std::vector<real_t> times, c_neg, c_pos;
    times.reserve(count);
    c_neg.reserve(count);
    c_pos.reserve(count);

    for (int i = begin; i < static_cast<int>(snapshots_->size()); ++i) {
        const auto& snap = (*snapshots_)[i];
        times.push_back(snap.time);
        c_neg.push_back(snap.c_s_surf_neg(0));
        c_pos.push_back(snap.c_s_surf_pos(0));
    }

    if (ImPlot::BeginPlot("颗粒表面浓度时序", ImVec2(-1, -1))) {
        ImPlot::SetupAxes("时间 [s]", "c_s_surf [mol/m^3]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(c_neg, c_pos);
        ImPlot::SetupFinish();
        ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 2.5f);
        ImPlot::PlotLine("负极（石墨）", times.data(), c_neg.data(), (int)times.size());
        ImPlot::SetNextMarkerStyle(ImPlotMarker_Circle, 2.5f);
        ImPlot::PlotLine("正极（NMC）", times.data(), c_pos.data(), (int)times.size());
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_stoichiometry_history() {
    if (snapshots_ == nullptr || snapshots_->empty()) {
        ImGui::TextDisabled("暂无仿真数据");
        return;
    }

    const int begin = history_window_begin(*snapshots_);
    const std::size_t count = snapshots_->size() - static_cast<std::size_t>(begin);
    std::vector<real_t> times, x, y;
    times.reserve(count);
    x.reserve(count);
    y.reserve(count);

    for (int i = begin; i < static_cast<int>(snapshots_->size()); ++i) {
        const auto& snap = (*snapshots_)[i];
        times.push_back(snap.time);
        x.push_back(snap.x_surf);
        y.push_back(snap.y_surf);
    }

    if (ImPlot::BeginPlot("颗粒表面嵌入率时序", ImVec2(-1, -1))) {
        ImPlot::SetupAxes("时间 [s]", "嵌入率 x / y",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(x, y);
        ImPlot::SetupFinish();
        ImPlot::PlotLine("负极 x", times.data(), x.data(), (int)times.size());
        ImPlot::PlotLine("正极 y", times.data(), y.data(), (int)times.size());
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_soc_history() {
    if (snapshots_ == nullptr || snapshots_->empty()) {
        ImGui::TextDisabled("暂无仿真数据");
        return;
    }

    const int begin = history_window_begin(*snapshots_);
    const std::size_t count = snapshots_->size() - static_cast<std::size_t>(begin);
    std::vector<real_t> times, soc;
    times.reserve(count);
    soc.reserve(count);

    for (int i = begin; i < static_cast<int>(snapshots_->size()); ++i) {
        const auto& snap = (*snapshots_)[i];
        times.push_back(snap.time);
        soc.push_back(snap.soc);
    }

    if (ImPlot::BeginPlot("荷电状态 SOC 时序", ImVec2(-1, -1))) {
        ImPlot::SetupAxes("时间 [s]", "SOC",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(soc);
        ImPlot::SetupFinish();
        ImPlot::PlotLine("SOC", times.data(), soc.data(), (int)times.size());
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_overpotential_history() {
    if (snapshots_ == nullptr || snapshots_->empty()) {
        ImGui::TextDisabled("暂无仿真数据");
        return;
    }

    const int begin = history_window_begin(*snapshots_);
    const std::size_t count = snapshots_->size() - static_cast<std::size_t>(begin);
    std::vector<real_t> times, eta_n, eta_p, eta_e;
    times.reserve(count);
    eta_n.reserve(count);
    eta_p.reserve(count);
    eta_e.reserve(count);

    for (int i = begin; i < static_cast<int>(snapshots_->size()); ++i) {
        const auto& snap = (*snapshots_)[i];
        times.push_back(snap.time);
        eta_n.push_back(snap.eta_neg);
        eta_p.push_back(snap.eta_pos);
        eta_e.push_back(snap.eta_e);
    }

    if (ImPlot::BeginPlot("过电位分解时序", ImVec2(-1, -1))) {
        ImPlot::SetupAxes("时间 [s]", "η [V]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(eta_n, eta_p, eta_e);
        ImPlot::SetupFinish();
        ImPlot::PlotLine("η_n（负极活化）", times.data(), eta_n.data(), (int)times.size());
        ImPlot::PlotLine("η_p（正极活化）", times.data(), eta_p.data(), (int)times.size());
        ImPlot::PlotLine("η_e（液相欧姆+浓差）", times.data(), eta_e.data(), (int)times.size());
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_reaction_current_profile() {
    if (snapshots_->empty() || snapshots_->back().i_surf_neg.size() == 0) {
        ImGui::TextDisabled("暂无孔壁电流数据（需 DFN 模型）");
        return;
    }

    const auto& snap = snapshots_->back();
    const int n_neg = static_cast<int>(snap.i_surf_neg.size());
    const int n_pos = static_cast<int>(snap.i_surf_pos.size());
    const std::vector<real_t> x_neg = electrode_x(n_neg);
    const std::vector<real_t> x_pos = electrode_x(n_pos);

    if (ImPlot::BeginPlot("孔壁电流密度分布（最新时刻）", ImVec2(-1, -1))) {
        ImPlot::SetupAxes("电极内位置（0 = 集流体，1 = 隔膜）", "i_surf [A/m^2]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(snap.i_surf_neg, snap.i_surf_pos);
        ImPlot::SetupFinish();
        ImPlot::PlotLine("负极 i_surf", x_neg.data(), snap.i_surf_neg.data(), n_neg);
        ImPlot::PlotLine("正极 i_surf", x_pos.data(), snap.i_surf_pos.data(), n_pos);
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_overpotential_profile() {
    if (snapshots_->empty() || snapshots_->back().eta_neg_x.size() == 0) {
        ImGui::TextDisabled("暂无过电位分布数据（需 DFN 模型）");
        return;
    }

    const auto& snap = snapshots_->back();
    const int n_neg = static_cast<int>(snap.eta_neg_x.size());
    const int n_pos = static_cast<int>(snap.eta_pos_x.size());
    const std::vector<real_t> x_neg = electrode_x(n_neg);
    const std::vector<real_t> x_pos = electrode_x(n_pos);

    if (ImPlot::BeginPlot("活化过电位分布（最新时刻）", ImVec2(-1, -1))) {
        ImPlot::SetupAxes("电极内位置（0 = 集流体，1 = 隔膜）", "η [V]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(snap.eta_neg_x, snap.eta_pos_x);
        ImPlot::SetupFinish();
        ImPlot::PlotLine("负极 η", x_neg.data(), snap.eta_neg_x.data(), n_neg);
        ImPlot::PlotLine("正极 η", x_pos.data(), snap.eta_pos_x.data(), n_pos);
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_ocp_profile() {
    if (snapshots_->empty() || snapshots_->back().u_neg_x.size() == 0) {
        ImGui::TextDisabled("暂无开路电位分布数据（需 DFN 模型）");
        return;
    }

    const auto& snap = snapshots_->back();
    const int n_neg = static_cast<int>(snap.u_neg_x.size());
    const int n_pos = static_cast<int>(snap.u_pos_x.size());
    const std::vector<real_t> x_neg = electrode_x(n_neg);
    const std::vector<real_t> x_pos = electrode_x(n_pos);

    if (ImPlot::BeginPlot("开路电位分布（最新时刻）", ImVec2(-1, -1))) {
        ImPlot::SetupAxes("电极内位置（0 = 集流体，1 = 隔膜）", "U_ocp [V]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(snap.u_neg_x, snap.u_pos_x);
        ImPlot::SetupFinish();
        ImPlot::PlotLine("负极 U_n", x_neg.data(), snap.u_neg_x.data(), n_neg);
        ImPlot::PlotLine("正极 U_p", x_pos.data(), snap.u_pos_x.data(), n_pos);
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_surface_conc_profile() {
    if (snapshots_->empty() || snapshots_->back().c_surf_neg_x.size() == 0) {
        ImGui::TextDisabled("暂无表面浓度分布数据（需 DFN 模型）");
        return;
    }

    const auto& snap = snapshots_->back();
    const int n_neg = static_cast<int>(snap.c_surf_neg_x.size());
    const int n_pos = static_cast<int>(snap.c_surf_pos_x.size());
    const std::vector<real_t> x_neg = electrode_x(n_neg);
    const std::vector<real_t> x_pos = electrode_x(n_pos);

    if (ImPlot::BeginPlot("颗粒表面浓度沿极片分布（最新时刻）", ImVec2(-1, -1))) {
        ImPlot::SetupAxes("电极内位置（0 = 集流体，1 = 隔膜）", "c_surf [mol/m^3]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(snap.c_surf_neg_x, snap.c_surf_pos_x);
        ImPlot::SetupFinish();
        ImPlot::PlotLine("负极 c_surf", x_neg.data(), snap.c_surf_neg_x.data(), n_neg);
        ImPlot::PlotLine("正极 c_surf", x_pos.data(), snap.c_surf_pos_x.data(), n_pos);
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_liquid_current_profile() {
    if (snapshots_->empty() || snapshots_->back().i_e_face.size() == 0) {
        ImGui::TextDisabled("暂无液相电流数据（需 DFN 模型）");
        return;
    }

    const auto& snap = snapshots_->back();
    const int n = static_cast<int>(snap.i_e_face.size());

    // 面电流定义在格界面上：均匀网格界面位置 i/(m)（0 → 1，含两端边界）
    std::vector<real_t> x(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        x[static_cast<std::size_t>(i)] =
            static_cast<real_t>(i) / static_cast<real_t>(n - 1);
    }

    if (ImPlot::BeginPlot("液相面电流分布（最新时刻）", ImVec2(-1, -1))) {
        ImPlot::SetupAxes("电池位置 [x/L]", "i_e [A/m^2]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(snap.i_e_face);
        ImPlot::SetupFinish();
        ImPlot::PlotLine("i_e", x.data(), snap.i_e_face.data(), n);
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_temperature_history() {
    if (snapshots_ == nullptr || snapshots_->empty()) {
        ImGui::TextDisabled("暂无仿真数据");
        return;
    }

    const int begin = history_window_begin(*snapshots_);
    const std::size_t count = snapshots_->size() - static_cast<std::size_t>(begin);
    std::vector<real_t> times, t_cell;
    times.reserve(count);
    t_cell.reserve(count);

    for (int i = begin; i < static_cast<int>(snapshots_->size()); ++i) {
        const auto& snap = (*snapshots_)[i];
        times.push_back(snap.time);
        t_cell.push_back(snap.temperature - 273.15); // K → ℃ 展示
    }

    if (ImPlot::BeginPlot("电芯温度时序", ImVec2(-1, -1))) {
        ImPlot::SetupAxes("时间 [s]", "T [℃]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(t_cell);
        ImPlot::SetupFinish();
        ImPlot::PlotLine("电芯温度", times.data(), t_cell.data(),
                         (int)times.size());
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_heat_generation_history() {
    if (snapshots_ == nullptr || snapshots_->empty()) {
        ImGui::TextDisabled("暂无仿真数据");
        return;
    }

    const int begin = history_window_begin(*snapshots_);
    const std::size_t count = snapshots_->size() - static_cast<std::size_t>(begin);
    std::vector<real_t> times, q_ohm, q_rxn, q_rev, q_total;
    times.reserve(count);
    q_ohm.reserve(count);
    q_rxn.reserve(count);
    q_rev.reserve(count);
    q_total.reserve(count);

    for (int i = begin; i < static_cast<int>(snapshots_->size()); ++i) {
        const auto& snap = (*snapshots_)[i];
        times.push_back(snap.time);
        q_ohm.push_back(snap.q_ohm);
        q_rxn.push_back(snap.q_rxn);
        q_rev.push_back(snap.q_rev);
        q_total.push_back(snap.q_total);
    }

    if (ImPlot::BeginPlot("产热功率分解（PyBaMM 三项口径）", ImVec2(-1, -1))) {
        ImPlot::SetupAxes("时间 [s]", "Q [W]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(q_ohm, q_rxn, q_rev, q_total);
        ImPlot::SetupFinish();
        ImPlot::SetNextLineStyle(IMPLOT_AUTO_COL, 2.5f);
        ImPlot::PlotLine("总产热 Q_total", times.data(), q_total.data(),
                         (int)times.size());
        ImPlot::PlotLine("欧姆热 Q_ohm", times.data(), q_ohm.data(),
                         (int)times.size());
        ImPlot::PlotLine("界面热 Q_rxn", times.data(), q_rxn.data(),
                         (int)times.size());
        ImPlot::PlotLine("熵热 Q_rev", times.data(), q_rev.data(),
                         (int)times.size());
        ImPlot::EndPlot();
    }
}

void PlotPanel::plot_degradation_history() {
    // 寿命曲线：每圈末快照（cycle_end）的可用容量与 SEI 厚度 vs 圈数。
    // 降解未启用时锂损失恒 0（容量恒标称、SEI 恒 0），需在仿真面板开启
    if (snapshots_ == nullptr || snapshots_->empty()) {
        ImGui::TextDisabled("暂无仿真数据");
        return;
    }

    std::vector<real_t> cycles, capacity, sei_nm;
    cycles.reserve(snapshots_->size());
    capacity.reserve(snapshots_->size());
    sei_nm.reserve(snapshots_->size());
    for (const auto& snap : *snapshots_) {
        if (!snap.cycle_end) {
            continue;
        }
        cycles.push_back(static_cast<real_t>(snap.cycle_number));
        capacity.push_back(snap.degraded_capacity_ah);
        sei_nm.push_back(snap.sei_thickness * 1e9f); // m → nm 展示
    }
    if (cycles.empty()) {
        ImGui::TextDisabled("暂无圈末采样点（需循环圈数 > 1 或等待首圈完成）");
        return;
    }

    if (ImPlot::BeginPlot("寿命曲线（每圈末采样）", ImVec2(-1, -1))) {
        ImPlot::SetupAxes("循环圈数", "可用容量 [Ah]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
        ImPlot::SetupAxis(ImAxis_Y2, "SEI 厚度 [nm]",
                          ImPlotAxisFlags_AutoFit);
        ImPlot::SetupFinish();
        ImPlot::SetAxes(ImAxis_X1, ImAxis_Y1);
        ImPlot::PlotLine("可用容量", cycles.data(), capacity.data(),
                         (int)cycles.size());
        ImPlot::SetAxes(ImAxis_X1, ImAxis_Y2);
        ImPlot::PlotLine("SEI 厚度", cycles.data(), sei_nm.data(),
                         (int)cycles.size());
        ImPlot::EndPlot();
    }
}

void PlotPanel::bind_history(
    const std::vector<SolutionSnapshot>* history) {
    // 只保存外部向量的地址：数据所有权与生命周期归调用方（MainWindow），
    // 两者同生命周期且同线程顺序访问（pump 后 render），无悬垂风险
    snapshots_ = history;
}

} // namespace fluxcell::gui
