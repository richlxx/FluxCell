/// @file main_window.cpp
/// @brief 主窗口实现
#include <fluxcell/gui/main_window.hpp>
#include <fluxcell/gui/plot_panel.hpp>
#include <fluxcell/gui/simulation_panel.hpp>
#include <fluxcell/gui/text_encoding.hpp>

#include <fluxcell/core/io/io.hpp>

#include <imgui.h>
#include <imgui_internal.h> // DockBuilder（默认停靠布局）
#include <ImGuiFileDialog.h>

#include <cstdio>

#include <algorithm>
#include <exception>
#include <iterator>
#include <utility>

namespace fluxcell::gui {

namespace {

/// CSV 导出保存对话框的唯一 ID
constexpr const char* kExportCsvDialogId = "ExportCsvDialog";

/// 历史快照点数上限：超过后对旧点隔点抽稀。100k 点 × ~1 KB
///（标量 + 保留的表面浓度向量）≈ 100 MB 内存封顶
constexpr std::size_t kMaxHistoryPoints = 100000;

/// 释放快照携带的分布场内存：Eigen 动态向量 clear/resize(0) 不归还
/// 缓冲，须与空向量交换。历史中仅最新快照需要分布场（空间分布图只画
/// 最新时刻），c_s_surf_neg/pos 保留（CSV 导出列与表面浓度时序共用）
void slim_profiles(SolutionSnapshot& s) {
    Vector().swap(s.c_e);
    Vector().swap(s.phi_s);
    Vector().swap(s.phi_e);
    Vector().swap(s.c_s_neg);
    Vector().swap(s.c_s_pos);
    Vector().swap(s.c_surf_neg_x);
    Vector().swap(s.c_surf_pos_x);
    Vector().swap(s.i_surf_neg);
    Vector().swap(s.i_surf_pos);
    Vector().swap(s.eta_neg_x);
    Vector().swap(s.eta_pos_x);
    Vector().swap(s.u_neg_x);
    Vector().swap(s.u_pos_x);
    Vector().swap(s.i_e_face);
}

/// 历史点数超限时对旧点隔点抽稀：点距随历史长度对数变粗（1s→2s→4s…），
/// 曲线形状保持；圈末点（寿命曲线采样）与最新点永不被抽掉
void thin_history(std::vector<SolutionSnapshot>& snaps) {
    if (snaps.size() <= kMaxHistoryPoints) {
        return;
    }
    std::vector<SolutionSnapshot> thinned;
    thinned.reserve(snaps.size() / 2 + snaps.size() / 16);
    for (std::size_t i = 0; i < snaps.size(); ++i) {
        if (i % 2 == 0 || snaps[i].cycle_end || i + 1 == snaps.size()) {
            thinned.push_back(std::move(snaps[i]));
        }
    }
    snaps = std::move(thinned);
}

} // namespace

MainWindow::MainWindow()
    : plot_panel_(std::make_unique<PlotPanel>()),
      simulation_panel_(std::make_unique<SimulationPanel>(engine_, spm_params_)) {
    // 绘图面板只读共享快照历史（绑定于本成员向量，生命周期覆盖面板）
    plot_panel_->bind_history(&snapshots_);
    // 重新开始/重置仿真时清理上一次的曲线数据（面板绑定的同一向量，
    // clear 即对面板可见）
    simulation_panel_->set_restart_callback([this] {
        snapshots_.clear();
        sim_status_ = SimulationStatus::Idle;
    });
}

MainWindow::~MainWindow() = default;

void MainWindow::render() {
    // 采集后台仿真输出
    pump_simulation();

    render_menu_bar();
    render_toolbar();
    // 状态栏先于 dockspace 渲染：以视口侧栏（底部）形式自动收缩工作区，
    // 停靠面板不再延伸到状态栏下方
    if (show_status_bar_) {
        render_status_bar();
    }
    render_dockspace();

    // 渲染面板（查看菜单控制可见性；模型类型每帧同步驱动图表清单联动）
    if (plot_panel_) {
        if (simulation_panel_) {
            plot_panel_->set_model_type(simulation_panel_->current_model());
        }
        for (int i = 0; i < PlotPanel::kPanelCount; ++i) {
            if (show_plot_panels_[i]) {
                plot_panel_->render_panel(i);
            }
        }
    }
    if (simulation_panel_ && show_simulation_panel_) {
        simulation_panel_->render();
    }

    // 导出对话框轮询 + 结果反馈 / 关于弹窗
    render_export_csv_dialog();
    render_notify_popup();
    render_about_popup();
}

void MainWindow::render_dockspace() {
    // 创建 DockSpace 占据整个视口
    ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::SetNextWindowPos(viewport->WorkPos);
    ImGui::SetNextWindowSize(viewport->WorkSize);
    ImGui::SetNextWindowViewport(viewport->ID);

    ImGuiWindowFlags flags = ImGuiWindowFlags_NoDocking
                           | ImGuiWindowFlags_NoTitleBar
                           | ImGuiWindowFlags_NoCollapse
                           | ImGuiWindowFlags_NoResize
                           | ImGuiWindowFlags_NoMove
                           | ImGuiWindowFlags_NoBringToFrontOnFocus
                           | ImGuiWindowFlags_NoNavFocus;

    ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
    ImGui::Begin("DockSpace Window", nullptr, flags);
    ImGui::PopStyleVar(3);

    // dockspace 名带 _v3 版本号：修正默认布局方向语义时作废旧 ini 记录
    ImGuiID dockspace_id = ImGui::GetID("MainDockSpace_v3");

    // 首次运行或用户请求恢复：左侧 30% 仿真控制 + 右侧 70% 2×2 图表网格；
    // 此后用户拖拽布局由 imgui.ini 持久化，不再重建
    if (rebuild_layout_ || ImGui::DockBuilderGetNode(dockspace_id) == nullptr) {
        ImGui::DockBuilderRemoveNode(dockspace_id);
        ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_PassthruCentralNode);
        ImGui::DockBuilderSetNodeSize(dockspace_id, viewport->WorkSize);

        // SplitNode 语义：返回 split_dir 方向上占 size_ratio 的节点，
        // out 反方向参数为剩余部分（此前返回值与剩余部分弄反，致仿真控制被
        // dock 到右侧且随 ini 持久化）
        ImGuiID dock_sim = 0, dock_charts = 0;
        ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Left, 0.30f,
                                    &dock_sim, &dock_charts);
        ImGuiID dock_top = 0, dock_bottom = 0;
        ImGui::DockBuilderSplitNode(dock_charts, ImGuiDir_Up, 0.5f,
                                    &dock_top, &dock_bottom);
        ImGuiID dock_tl = 0, dock_tr = 0, dock_bl = 0, dock_br = 0;
        ImGui::DockBuilderSplitNode(dock_top, ImGuiDir_Left, 0.5f, &dock_tl, &dock_tr);
        ImGui::DockBuilderSplitNode(dock_bottom, ImGuiDir_Left, 0.5f, &dock_bl, &dock_br);

        ImGui::DockBuilderDockWindow("仿真控制", dock_sim);
        ImGui::DockBuilderDockWindow("图表 1", dock_tl);
        ImGui::DockBuilderDockWindow("图表 2", dock_tr);
        ImGui::DockBuilderDockWindow("图表 3", dock_bl);
        ImGui::DockBuilderDockWindow("图表 4", dock_br);
        ImGui::DockBuilderFinish(dockspace_id);

        rebuild_layout_ = false;
        show_simulation_panel_ = true;
        std::fill(std::begin(show_plot_panels_), std::end(show_plot_panels_), true);
        show_status_bar_ = true;
        show_toolbar_ = true;
    }

    ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f),
                     ImGuiDockNodeFlags_PassthruCentralNode);

    ImGui::End();
}

void MainWindow::pump_simulation() {
    auto fresh = engine_.collect_snapshots();
    if (!fresh.empty()) {
        // 旧「最新」快照让出分布场：空间分布图只画最新时刻，历史快照
        // 的 14 个分布向量（DFN ~22 KB/步）是长循环内存的大头——
        // 不释放的话 200 圈工况可达数十 GB（内存占用 99% 的根因）
        if (!snapshots_.empty()) {
            slim_profiles(snapshots_.back());
        }
        snapshots_.insert(snapshots_.end(),
                          std::make_move_iterator(fresh.begin()),
                          std::make_move_iterator(fresh.end()));
        thin_history(snapshots_);
    }
    sim_status_ = engine_.status();
}

void MainWindow::render_menu_bar() {
    // Ctrl+O：参数文件选择对话框（与菜单项等效）
    const ImGuiIO& io = ImGui::GetIO();
    if (io.KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_O, false)
        && simulation_panel_) {
        simulation_panel_->open_file_dialog();
    }

    // Ctrl+S：参数保存对话框（与菜单项等效）
    if (io.KeyCtrl && ImGui::IsKeyPressed(ImGuiKey_S, false)
        && simulation_panel_) {
        simulation_panel_->open_save_file_dialog();
    }

    // Ctrl+Shift+O / Ctrl+Shift+S：工况载入/保存对话框（与菜单项等效；
    // 无 Shift 修饰的 O/S 留给参数文件）
    if (io.KeyCtrl && io.KeyShift && ImGui::IsKeyPressed(ImGuiKey_O, false)
        && simulation_panel_) {
        simulation_panel_->open_protocol_dialog();
    }
    if (io.KeyCtrl && io.KeyShift && ImGui::IsKeyPressed(ImGuiKey_S, false)
        && simulation_panel_) {
        simulation_panel_->open_save_protocol_dialog();
    }

    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("文件")) {
            if (ImGui::MenuItem("加载参数...", "Ctrl+O")) {
                if (simulation_panel_) {
                    simulation_panel_->open_file_dialog();
                }
            }
            if (ImGui::MenuItem("保存参数...", "Ctrl+S")) {
                if (simulation_panel_) {
                    simulation_panel_->open_save_file_dialog();
                }
            }
            ImGui::Separator();
            if (ImGui::MenuItem("载入工况...", "Ctrl+Shift+O")) {
                if (simulation_panel_) {
                    simulation_panel_->open_protocol_dialog();
                }
            }
            if (ImGui::MenuItem("保存工况...", "Ctrl+Shift+S")) {
                if (simulation_panel_) {
                    simulation_panel_->open_save_protocol_dialog();
                }
            }
            ImGui::Separator();
            // 无仿真结果时禁用导出
            if (ImGui::MenuItem("导出结果 (CSV)...", nullptr, false,
                                !snapshots_.empty())) {
                open_export_csv_dialog();
            }
            ImGui::Separator();
            if (ImGui::MenuItem("退出", "Alt+F4")) {
                if (on_exit_) {
                    on_exit_();
                }
            }
            ImGui::EndMenu();
        }

        if (ImGui::BeginMenu("仿真")) {
            const bool running = sim_status_ == SimulationStatus::Running
                              || sim_status_ == SimulationStatus::Paused;
            if (ImGui::MenuItem("运行", "F5", false, !running) && simulation_panel_) {
                simulation_panel_->start();
            }
            const bool can_step = sim_status_ == SimulationStatus::Idle
                               || sim_status_ == SimulationStatus::Paused;
            if (ImGui::MenuItem("单步", "F10", false, can_step) && simulation_panel_) {
                simulation_panel_->step_once();
            }
            if (ImGui::MenuItem("暂停/继续", "F6", false, running) && simulation_panel_) {
                simulation_panel_->pause_resume();
            }
            if (ImGui::MenuItem("停止", "Shift+F5", false, running) && simulation_panel_) {
                simulation_panel_->stop();
            }
            ImGui::Separator();
            if (ImGui::MenuItem("重置")) {
                // reset() 内部通过重启回调清理曲线
                if (simulation_panel_) simulation_panel_->reset();
                sim_status_ = SimulationStatus::Idle;
            }
            ImGui::EndMenu();
        }

        render_view_menu();

        if (ImGui::BeginMenu("帮助")) {
            if (ImGui::MenuItem("关于 FluxCell")) {
                show_about_ = true;
            }
            ImGui::EndMenu();
        }

        ImGui::EndMainMenuBar();
    }
}

void MainWindow::render_toolbar() {
    if (!show_toolbar_) {
        return;
    }
    ImGuiViewport* viewport = ImGui::GetMainViewport();
    const ImGuiStyle& style = ImGui::GetStyle();
    // 一行控件高度 + 窗口上下内边距（与状态栏同口径，避免文字垂直裁剪）
    const float toolbar_height = ImGui::GetFrameHeight()
                               + style.WindowPadding.y * 2.0f;

    // 视口侧栏：自动叠加进工作区（菜单栏下方、dockspace 上方占位），
    // 内部基于 GetBuildWorkRect() 定位，同帧内先菜单栏后工具栏依次排列
    ImGuiWindowFlags flags = ImGuiWindowFlags_NoScrollbar
                           | ImGuiWindowFlags_NoSavedSettings;
    if (ImGui::BeginViewportSideBar("##MainToolbar", viewport, ImGuiDir_Up,
                                    toolbar_height, flags)) {
        // 按钮在一行高度内垂直居中
        ImGui::SetCursorPosY((toolbar_height - ImGui::GetFrameHeight()) * 0.5f);

        if (simulation_panel_) {
            const bool running = sim_status_ == SimulationStatus::Running;
            const bool paused = sim_status_ == SimulationStatus::Paused;
            const bool active = running || paused;

            // 运行：空闲/暂停/已结束可用
            ImGui::BeginDisabled(running);
            if (ImGui::Button("运行")) {
                simulation_panel_->start();
            }
            ImGui::EndDisabled();
            ImGui::SetItemTooltip("启动仿真 (F5)");

            // 暂停/继续：运行中显示"暂停"，已暂停显示"继续"
            ImGui::SameLine();
            ImGui::BeginDisabled(!active);
            if (ImGui::Button(paused ? "继续" : "暂停")) {
                simulation_panel_->pause_resume();
            }
            ImGui::EndDisabled();
            ImGui::SetItemTooltip("暂停/继续仿真 (F6)");

            // 停止
            ImGui::SameLine();
            ImGui::BeginDisabled(!active);
            if (ImGui::Button("停止")) {
                simulation_panel_->stop();
            }
            ImGui::EndDisabled();
            ImGui::SetItemTooltip("停止仿真 (Shift+F5)");

            // 单步：空闲/暂停可用
            ImGui::SameLine();
            ImGui::BeginDisabled(!(paused || sim_status_ == SimulationStatus::Idle));
            if (ImGui::Button("单步")) {
                simulation_panel_->step_once();
            }
            ImGui::EndDisabled();
            ImGui::SetItemTooltip("推进一个时间步 (F10)");

            // 分组分隔（侧栏窗口非菜单栏，Separator() 是水平线，需用内部竖线）
            ImGui::SameLine();
            ImGui::SeparatorEx(ImGuiSeparatorFlags_Vertical);

            // 重置：清空曲线回到初始状态
            ImGui::SameLine();
            if (ImGui::Button("重置")) {
                simulation_panel_->reset();
                sim_status_ = SimulationStatus::Idle;
            }
            ImGui::SetItemTooltip("清空结果并复位仿真");
        }
    }
    ImGui::End();
}

void MainWindow::render_view_menu() {
    if (ImGui::BeginMenu("查看")) {
        ImGui::MenuItem("仿真控制面板", nullptr, &show_simulation_panel_);
        for (int i = 0; i < PlotPanel::kPanelCount; ++i) {
            char label[32];
            std::snprintf(label, sizeof(label), "图表面板 %d", i + 1);
            ImGui::MenuItem(label, nullptr, &show_plot_panels_[i]);
        }
        ImGui::Separator();
        ImGui::MenuItem("工具栏", nullptr, &show_toolbar_);
        ImGui::MenuItem("状态栏", nullptr, &show_status_bar_);
        ImGui::Separator();
        if (ImGui::MenuItem("恢复默认布局")) {
            rebuild_layout_ = true;
        }
        ImGui::EndMenu();
    }
}

void MainWindow::render_status_bar() {
    ImGuiViewport* viewport = ImGui::GetMainViewport();
    const ImGuiStyle& style = ImGui::GetStyle();
    // 窗口高度 = 一行文字 + 窗口上下内边距。
    // 注意不能用 GetFrameHeight()：它不含 WindowPadding，文字会被垂直裁剪。
    const float status_height = ImGui::GetTextLineHeight() + style.WindowPadding.y * 2.0f;

    // 视口侧栏（底部）：与工具栏同走 BeginViewportSideBar，自动收缩工作区，
    // 停靠面板不会延伸到状态栏下方；此前 SetNextWindowPos 覆盖式定位
    // 仅叠加显示，会遮挡面板底部内容
    ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoSavedSettings
                           | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoInputs
                           | ImGuiWindowFlags_NoFocusOnAppearing;

    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.93f, 0.93f, 0.95f, 1.0f));
    // Begin 返回值可安全忽略（false 时 ImGui 内部跳过内容绘制），End 无条件配对
    ImGui::BeginViewportSideBar("##StatusBar", viewport, ImGuiDir_Down,
                                status_height, flags);

    // 状态指示（颜色在亮色背景上保持对比度）
    const char* status_text = "空闲";
    ImVec4 status_color(0.45f, 0.45f, 0.45f, 1.0f);
    switch (sim_status_) {
        case SimulationStatus::Running:
            status_text = "运行中"; status_color = ImVec4(0.05f, 0.60f, 0.10f, 1.0f); break;
        case SimulationStatus::Paused:
            status_text = "已暂停"; status_color = ImVec4(0.72f, 0.60f, 0.00f, 1.0f); break;
        case SimulationStatus::Completed:
            status_text = "已完成"; status_color = ImVec4(0.10f, 0.40f, 0.85f, 1.0f); break;
        case SimulationStatus::Failed:
            status_text = "失败"; status_color = ImVec4(0.85f, 0.15f, 0.15f, 1.0f); break;
        default: break;
    }

    // 预先计算右侧版本号区域；左侧信息逐项检查剩余空间，
    // 避免换行到第二行（单行高度窗口会裁剪）或与版本号重叠
    const char* version_text = "FluxCell v0.1.0";
    const float version_w = ImGui::CalcTextSize(version_text).x;
    const float version_x = ImGui::GetWindowWidth() - version_w - style.WindowPadding.x;
    // 当前光标位置追加一段文本后是否仍位于右边界内
    auto fits_before_version = [&](const char* text) {
        return ImGui::GetCursorPosX() + style.ItemSpacing.x + ImGui::CalcTextSize(text).x
               <= version_x;
    };

    ImGui::TextColored(status_color, "[%s]", status_text);

    // 运行/暂停时：迷你进度条 + 百分比（手绘以适配单行状态栏高度）
    if (sim_status_ == SimulationStatus::Running
        || sim_status_ == SimulationStatus::Paused) {
        const real_t prog = std::clamp(engine_.progress(), 0.0, 1.0);
        char pct_buf[32];
        std::snprintf(pct_buf, sizeof(pct_buf), "%.1f%%", prog * 100.0);
        if (fits_before_version(pct_buf)) {
            ImGui::SameLine();
            const float bar_w = 80.0f;
            const float bar_h = ImGui::GetTextLineHeight();
            const ImVec2 bar_p = ImGui::GetCursorScreenPos();
            ImDrawList* dl = ImGui::GetWindowDrawList();
            // 背景槽（半透明）+ 按状态色填充（运行绿 / 暂停黄）
            dl->AddRectFilled(bar_p, ImVec2(bar_p.x + bar_w, bar_p.y + bar_h),
                              IM_COL32(0, 0, 0, 30), 2.0f);
            dl->AddRectFilled(
                bar_p,
                ImVec2(bar_p.x + bar_w * static_cast<float>(prog),
                       bar_p.y + bar_h),
                ImGui::GetColorU32(status_color), 2.0f);
            ImGui::Dummy(ImVec2(bar_w, bar_h));
            ImGui::SameLine();
            ImGui::TextUnformatted(pct_buf);
        }
    }

    ImGui::SameLine();
    ImGui::Text("步数: %zu", snapshots_.size());

    // 当前电压 / 时间（空间不足时优先隐藏更次要的信息）
    if (!snapshots_.empty()) {
        char buf_v[64];
        char buf_t[64];
        std::snprintf(buf_v, sizeof(buf_v), "电压 = %.4f V", snapshots_.back().voltage);
        std::snprintf(buf_t, sizeof(buf_t), "时间 = %.0f s", snapshots_.back().time);
        if (fits_before_version(buf_v)) {
            ImGui::SameLine();
            ImGui::Text("%s", buf_v);
            if (fits_before_version(buf_t)) {
                ImGui::SameLine();
                ImGui::Text("%s", buf_t);
            }
        }
    }

    // 版本号右对齐；仅当左侧信息未侵入其区域时绘制
    if (ImGui::GetCursorPosX() + style.ItemSpacing.x <= version_x) {
        ImGui::SameLine(version_x);
        ImGui::TextDisabled("%s", version_text);
    }

    ImGui::End();
    ImGui::PopStyleColor();
}

void MainWindow::open_export_csv_dialog() {
    IGFD::FileDialogConfig config;
    config.path = "results";
    config.fileName = "results.csv";
    config.countSelectionMax = 1;
    config.flags = ImGuiFileDialogFlags_Modal
                 | ImGuiFileDialogFlags_ConfirmOverwrite;
    ImGuiFileDialog::Instance()->OpenDialog(kExportCsvDialogId,
                                            "导出结果 (CSV)",
                                            ".csv", config);
}

void MainWindow::render_export_csv_dialog() {
    if (!ImGuiFileDialog::Instance()->Display(kExportCsvDialogId)) {
        return;
    }
    if (ImGuiFileDialog::Instance()->IsOk() && !snapshots_.empty()) {
        try {
            // IGFD 返回本地代码页路径 → 转 UTF-8 显示；io 层再转回本地写文件
            const std::string picked =
                acp_to_utf8(ImGuiFileDialog::Instance()->GetFilePathName());
            io::export_csv(utf8_to_acp(picked), snapshots_, CellGeometry{});
            notify_message_ = "已导出 " + std::to_string(snapshots_.size())
                            + " 个快照\n" + picked;
        } catch (const std::exception& e) {
            notify_message_ = std::string("导出失败: ") + e.what();
        }
        notify_pending_ = true;
    }
    ImGuiFileDialog::Instance()->Close();
}

void MainWindow::render_notify_popup() {
    if (notify_pending_) {
        ImGui::OpenPopup("导出结果");
        notify_pending_ = false;
    }
    if (ImGui::BeginPopupModal("导出结果", nullptr,
                               ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::TextUnformatted(notify_message_.c_str());
        ImGui::Spacing();
        if (ImGui::Button("确定", ImVec2(120.0f, 0.0f))
            || ImGui::IsKeyPressed(ImGuiKey_Escape)) {
            notify_message_.clear();
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
}

void MainWindow::render_about_popup() {
    if (show_about_) {
        ImGui::OpenPopup("关于 FluxCell");
        show_about_ = false;
    }
    if (ImGui::BeginPopupModal("关于 FluxCell", nullptr,
                               ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::TextDisabled("FluxCell v0.1.0");
        ImGui::Separator();
        ImGui::Text("锂电池电化学仿真软件");
        ImGui::Spacing();
        ImGui::BulletText("模型：单颗粒 SPM、扩展 SPMe、全阶伪二维 DFN (P2D)");
        ImGui::BulletText("求解：Crank-Nicolson 固定步长 / SUNDIALS IDA 自适应");
        ImGui::BulletText("参数：Chen2020 (LG M50)，对标 PyBaMM 验证");
        ImGui::BulletText("输运：Nyman2008 浓度相关电解液 kappa(c)/D(c)");
        ImGui::Spacing();
        ImGui::TextDisabled("C++20 · Dear ImGui · ImPlot · SUNDIALS IDA · yaml-cpp");
        ImGui::Spacing();
        if (ImGui::Button("关闭", ImVec2(120.0f, 0.0f))
            || ImGui::IsKeyPressed(ImGuiKey_Escape)) {
            ImGui::CloseCurrentPopup();
        }
        ImGui::EndPopup();
    }
}

} // namespace fluxcell::gui
