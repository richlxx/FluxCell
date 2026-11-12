/// @file simulation_panel.cpp
/// @brief 仿真控制面板实现
#include <fluxcell/gui/simulation_panel.hpp>

#include <fluxcell/core/io/io.hpp>
#include <fluxcell/gui/plot_util.hpp>
#include <fluxcell/gui/text_encoding.hpp>

#include <imgui.h>
#include <implot.h>
#include <ImGuiFileDialog.h>

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <exception>
#include <vector>

#include <algorithm>

namespace fluxcell::gui {

namespace {

/// 参数文件对话框的唯一 ID（加载 / 保存各一）
constexpr const char* kParamsDialogId = "LoadParamsYamlDialog";
constexpr const char* kSaveParamsDialogId = "SaveParamsYamlDialog";

/// 工况文件对话框的唯一 ID（载入 / 保存各一）
constexpr const char* kProtocolDialogId = "LoadProtocolYamlDialog";
constexpr const char* kSaveProtocolDialogId = "SaveProtocolYamlDialog";

/// 电芯参数表行：名称列 + 只读值列（单值）
void cell_param_row(const char* label, const char* fmt, double v) {
    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::TextUnformatted(label);
    ImGui::TableNextColumn();
    ImGui::TextDisabled(fmt, v);
}

/// 电芯参数表行：名称列 + 只读值列（双值，如 H×W、x@0%/x@100%）
void cell_param_row2(const char* label, const char* fmt, double v1, double v2) {
    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::TextUnformatted(label);
    ImGui::TableNextColumn();
    ImGui::TextDisabled(fmt, v1, v2);
}

/// 电芯参数表行：名称列 + 只读值列（三值，如各层 ρcp）
void cell_param_row3(const char* label, const char* fmt,
                     double v1, double v2, double v3) {
    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::TextUnformatted(label);
    ImGui::TableNextColumn();
    ImGui::TextDisabled(fmt, v1, v2, v3);
}

/// InputText 的 std::string 扩容回调（ImGui 官方配方）：缓冲增长时同步
/// 底层 std::string，避免为失败信息预留固定 char[]
int input_text_resize_callback(ImGuiInputTextCallbackData* data) {
    if (data->EventFlag == ImGuiInputTextFlags_CallbackResize) {
        auto* s = static_cast<std::string*>(data->UserData);
        s->resize(data->BufTextLen);
        data->Buf = s->data();
    }
    return 0;
}

} // namespace

SimulationPanel::SimulationPanel(SimulationEngine& engine,
                                 model::SPMModel::Parameters& params)
    : engine_(engine), params_(params) {
    steps_.push_back(StepUI{}); // 默认：单段恒流放电
}

void SimulationPanel::set_restart_callback(std::function<void()> cb) {
    on_restart_ = std::move(cb);
}

void SimulationPanel::load_params() {
    try {
        // 统一用 SPMe 读取：固相键与 SPM 完全一致，电解液键缺失时取默认值
        // DFN 读取：继承 SPMe 结果，另取固相电导键（缺失取默认值）
        // Windows 窄字符 API 以本地代码页解析路径：GUI 内部统一 UTF-8，
        // 传给 std::ifstream 前转换；其余平台原样
        const std::string local = utf8_to_acp(params_path_);
        spme_params_ = io::load_spme_parameters(local);
        dfn_params_ = io::load_dfn_parameters(local);
        // 热模型选型与参数文件同步（DFN 专属；SPM/SPMe 忽略）
        thermal_type_ = dfn_params_.thermal_type == model::ThermalType::Lumped
                            ? 1 : 0;
        params_ = spme_params_.solid;
        load_error_.clear();
    } catch (const std::exception& e) {
        load_error_ = e.what();
    }
}

void SimulationPanel::save_params() {
    try {
        // 与 load_params 对称：GUI 内部 UTF-8 路径转本地代码页后交给 io 层；
        // 按 UI 当前模型类型写出对应层级的参数（文件口径，不含面板运行条件）
        const std::string local = utf8_to_acp(params_path_);
        switch (model_type_) {
            case 1:  io::save_spme_parameters(local, spme_params_); break;
            case 2:  io::save_dfn_parameters(local, dfn_params_); break;
            default: io::save_spm_parameters(local, params_); break;
        }
        save_error_.clear();
    } catch (const std::exception& e) {
        save_error_ = e.what();
    }
}

void SimulationPanel::load_protocol() {
    try {
        // GUI 内部 UTF-8 路径转本地代码页后交给 io 层（与参数文件同口径）
        const std::string local = utf8_to_acp(protocol_path_);
        const io::ProtocolFile pf = io::load_protocol(local);

        // 段序列：CCStep → StepUI（Mode/Control 枚举顺序两两对应，直接整型强转）
        steps_.clear();
        for (const auto& s : pf.protocol.steps) {
            StepUI ui;
            ui.mode = static_cast<int>(s.mode);
            ui.control = static_cast<int>(s.control);
            ui.current = static_cast<float>(s.current);
            ui.voltage = static_cast<float>(s.voltage);
            ui.power = static_cast<float>(s.power);
            ui.cutoff_time = static_cast<float>(s.cutoff_time);
            ui.cutoff_voltage = static_cast<float>(s.cutoff_voltage);
            ui.cutoff_current = static_cast<float>(s.cutoff_current);
            steps_.push_back(ui);
        }
        if (steps_.empty()) {
            steps_.push_back(StepUI{}); // 与构造函数同兑底：至少一个段
        }

        // 运行条件与模型/求解器选择
        cycles_ = pf.protocol.cycles;
        time_step_ = static_cast<float>(pf.protocol.time_step);
        real_time_factor_ = static_cast<float>(pf.protocol.real_time_factor);
        model_type_ = static_cast<int>(pf.protocol.model);
        solver_type_ = static_cast<int>(pf.protocol.solver);
        solver_rel_tol_ = static_cast<float>(pf.protocol.solver_rel_tol);
        initial_soc_pct_ = static_cast<float>(pf.initial_soc * 100.0);
        ambient_temp_c_ = static_cast<float>(pf.ambient_temp_c);
        // 热模型选型与文件同步（DFN 专属；与 load_params 的同步逻辑一致）
        thermal_type_ = pf.thermal_type == model::ThermalType::Lumped ? 1 : 0;
        dfn_params_.thermal_type = pf.thermal_type;

        // 载入的合法工况清掉旧的校验提示与载入错误
        protocol_error_.clear();
        protocol_load_error_.clear();
    } catch (const std::exception& e) {
        protocol_load_error_ = e.what();
    }
}

void SimulationPanel::save_protocol() {
    try {
        // 与 load_protocol 对称：UI 缓存 → ProtocolFile（含运行条件全量）
        const std::string local = utf8_to_acp(protocol_path_);
        io::ProtocolFile pf;
        pf.protocol = build_protocol();
        pf.initial_soc = static_cast<real_t>(initial_soc_pct_) / 100.0;
        pf.ambient_temp_c = static_cast<real_t>(ambient_temp_c_);
        pf.thermal_type = thermal_type_ == 1 ? model::ThermalType::Lumped
                                             : model::ThermalType::Isothermal;
        io::save_protocol(local, pf);
        protocol_save_error_.clear();
    } catch (const std::exception& e) {
        protocol_save_error_ = e.what();
    }
}

void SimulationPanel::open_file_dialog() {
    // 初始定位到当前路径的目录与文件名（空/无分隔符则落到 params/）
    const std::string current = params_path_;
    const std::string local = utf8_to_acp(current);
    const auto sep = local.find_last_of("/\\");

    IGFD::FileDialogConfig config;
    if (sep == std::string::npos) {
        config.path = "params";
        config.fileName = local;
    } else {
        config.path = local.substr(0, sep);
        config.fileName = local.substr(sep + 1);
    }
    config.countSelectionMax = 1;
    config.flags = ImGuiFileDialogFlags_Modal;
    ImGuiFileDialog::Instance()->OpenDialog(kParamsDialogId,
                                            "加载参数文件 (YAML)",
                                            ".yml,.yaml", config);
}

void SimulationPanel::open_save_file_dialog() {
    // 初始定位到当前路径的目录与文件名（空/无分隔符则落到 params/）
    const std::string current = params_path_;
    const std::string local = utf8_to_acp(current);
    const auto sep = local.find_last_of("/\\");

    IGFD::FileDialogConfig config;
    if (sep == std::string::npos) {
        config.path = "params";
        config.fileName = local;
    } else {
        config.path = local.substr(0, sep);
        config.fileName = local.substr(sep + 1);
    }
    config.countSelectionMax = 1;
    config.flags = ImGuiFileDialogFlags_Modal
                 | ImGuiFileDialogFlags_ConfirmOverwrite;
    ImGuiFileDialog::Instance()->OpenDialog(kSaveParamsDialogId,
                                            "保存参数文件 (YAML)",
                                            ".yml,.yaml", config);
}

void SimulationPanel::open_protocol_dialog() {
    // 初始定位到当前路径的目录与文件名（空/无分隔符则落到 params/）
    const std::string local = utf8_to_acp(protocol_path_);
    const auto sep = local.find_last_of("/\\");

    IGFD::FileDialogConfig config;
    if (sep == std::string::npos) {
        config.path = "params";
        config.fileName = local;
    } else {
        config.path = local.substr(0, sep);
        config.fileName = local.substr(sep + 1);
    }
    config.countSelectionMax = 1;
    config.flags = ImGuiFileDialogFlags_Modal;
    ImGuiFileDialog::Instance()->OpenDialog(kProtocolDialogId,
                                            "载入工况文件 (YAML)",
                                            ".yml,.yaml", config);
}

void SimulationPanel::open_save_protocol_dialog() {
    // 初始定位到当前路径的目录与文件名（空/无分隔符则落到 params/）
    const std::string local = utf8_to_acp(protocol_path_);
    const auto sep = local.find_last_of("/\\");

    IGFD::FileDialogConfig config;
    if (sep == std::string::npos) {
        config.path = "params";
        config.fileName = local;
    } else {
        config.path = local.substr(0, sep);
        config.fileName = local.substr(sep + 1);
    }
    config.countSelectionMax = 1;
    config.flags = ImGuiFileDialogFlags_Modal
                 | ImGuiFileDialogFlags_ConfirmOverwrite;
    ImGuiFileDialog::Instance()->OpenDialog(kSaveProtocolDialogId,
                                            "保存工况文件 (YAML)",
                                            ".yml,.yaml", config);
}

void SimulationPanel::render() {
    // 横向滚动条：工况行控件较多超宽时自动出现，避免控件被裁剪
    ImGui::Begin("仿真控制", nullptr, ImGuiWindowFlags_HorizontalScrollbar);

    render_protocol_section();
    // 电芯参数：独立折叠区（YAML 只读，三模型共享）
    render_cell_params();
    render_model_params_section();
    ImGui::Separator();
    // 控制按钮（运行/暂停/停止/单步/重置）已由主窗口工具栏承担
    render_progress_section();

    ImGui::End();

    // 参数文件选择对话框（模态；每帧调用直至关闭）
    if (ImGuiFileDialog::Instance()->Display(kParamsDialogId)) {
        if (ImGuiFileDialog::Instance()->IsOk()) {
            // 对话框返回本地代码页路径：转 UTF-8 供界面显示；
            // load_params 内部再转回本地代码页传给 std::ifstream
            const std::string picked =
                acp_to_utf8(ImGuiFileDialog::Instance()->GetFilePathName());
            std::snprintf(params_path_, sizeof(params_path_), "%s",
                          picked.c_str());
            load_params();
        }
        ImGuiFileDialog::Instance()->Close();
    }

    // 参数保存对话框（模态；每帧调用直至关闭）
    if (ImGuiFileDialog::Instance()->Display(kSaveParamsDialogId)) {
        if (ImGuiFileDialog::Instance()->IsOk()) {
            const std::string picked =
                acp_to_utf8(ImGuiFileDialog::Instance()->GetFilePathName());
            std::snprintf(params_path_, sizeof(params_path_), "%s",
                          picked.c_str());
            save_params();
        }
        ImGuiFileDialog::Instance()->Close();
    }

    // 工况载入对话框（模态；每帧调用直至关闭）
    if (ImGuiFileDialog::Instance()->Display(kProtocolDialogId)) {
        if (ImGuiFileDialog::Instance()->IsOk()) {
            // 对话框返回本地代码页路径：转 UTF-8 供界面显示；
            // load_protocol 内部再转回本地代码页传给 io 层
            const std::string picked =
                acp_to_utf8(ImGuiFileDialog::Instance()->GetFilePathName());
            std::snprintf(protocol_path_, sizeof(protocol_path_), "%s",
                          picked.c_str());
            load_protocol();
        }
        ImGuiFileDialog::Instance()->Close();
    }

    // 工况保存对话框（模态；每帧调用直至关闭）
    if (ImGuiFileDialog::Instance()->Display(kSaveProtocolDialogId)) {
        if (ImGuiFileDialog::Instance()->IsOk()) {
            const std::string picked =
                acp_to_utf8(ImGuiFileDialog::Instance()->GetFilePathName());
            std::snprintf(protocol_path_, sizeof(protocol_path_), "%s",
                          picked.c_str());
            save_protocol();
        }
        ImGuiFileDialog::Instance()->Close();
    }
}

void SimulationPanel::render_protocol_section() {
    // 纵向配置行的标签列起点（标签在左、控件对齐在右侧统一 x）
    const float label_col = ImGui::GetCursorPosX() + 100.0f;

    // 模型选择
    ImGui::TextUnformatted("模型");
    ImGui::SameLine(label_col);
    const char* model_names[] = {"SPM（单颗粒）", "SPMe（含电解液）", "DFN（全阶 P2D）"};
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::Combo("##model", &model_type_, model_names, IM_ARRAYSIZE(model_names));
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip(model_type_ == 1
            ? "SPMe：电解液 1D 浓度/电位分布，液相欧姆+浓差损失"
            : model_type_ == 2
                ? "DFN：每 x 位置独立颗粒，液相浓度/电位 1D 分布，固相欧姆损失，全阶 P2D"
                : "SPM：电解液均匀，无液相损失");
    }

    // 求解器选择（不同模型均可选；IDA 要求模型提供 DAE 视图）
    ImGui::TextUnformatted("求解器");
    ImGui::SameLine(label_col);
    const char* solver_names[] = {"固定步长 (Crank-Nicolson)", "自适应 (SUNDIALS IDA)"};
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::Combo("##solver", &solver_type_, solver_names, IM_ARRAYSIZE(solver_names));
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip(solver_type_ == 1
            ? "IDA 自动调整内部步长（BDF 多阶隐式），输出点由用户设定"
            : "各模型内置的 Crank-Nicolson 二阶精度固定步长积分");
    }
    if (solver_type_ == 1) {
        ImGui::TextUnformatted("相对容差");
        ImGui::SameLine(label_col);
        ImGui::SetNextItemWidth(-1.0f);
        ImGui::InputFloat("##reltol", &solver_rel_tol_, 0.0f, 0.0f, "%.1e");
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("IDA 局部误差相对容差（1e-6 ~ 1e-10）");
        }
    }

    // 参数文件：手输（回车确认）或点击“加载…”弹出文件选择对话框
    ImGui::TextUnformatted("参数文件");
    ImGui::SameLine(label_col);
    // 输入框宽度为右侧“加载…”按钮留位：负宽度按按钮实际宽 + 间距
    // 计算，避免按钮溢出面板右边界被裁剪
    const ImGuiStyle& style = ImGui::GetStyle();
    const float load_btn_w = ImGui::CalcTextSize("加载…").x
                           + style.FramePadding.x * 2.0f;
    ImGui::SetNextItemWidth(-load_btn_w - style.ItemSpacing.x);
    if (ImGui::InputText("##params", params_path_, sizeof(params_path_),
                         ImGuiInputTextFlags_EnterReturnsTrue)) {
        load_params();
    }
    ImGui::SameLine();
    if (ImGui::Button("加载…")) {
        open_file_dialog();
    }
    if (!load_error_.empty()) {
        ImGui::TextColored(ImVec4(0.9f, 0.3f, 0.3f, 1.0f), "加载失败: %s",
                           load_error_.c_str());
    }
    if (!save_error_.empty()) {
        ImGui::TextColored(ImVec4(0.9f, 0.3f, 0.3f, 1.0f), "保存失败: %s",
                           save_error_.c_str());
    }

    // 仿真条件
    ImGui::TextUnformatted("初始 SOC");
    ImGui::SameLine(label_col);
    ImGui::SetNextItemWidth(-1.0f);
    if (ImGui::InputFloat("##soc", &initial_soc_pct_, 1.0f, 10.0f, "%.0f %%")) {
        initial_soc_pct_ = std::clamp(initial_soc_pct_, 0.0f, 100.0f);
    }
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("0~100%%；决定初始嵌入率：石墨 x=%.4f→%.4f，NMC y=%.4f→%.4f",
                          static_cast<double>(params_.neg_stoich_at_soc0),
                          static_cast<double>(params_.neg_stoich_at_soc1),
                          static_cast<double>(params_.pos_stoich_at_soc1),
                          static_cast<double>(params_.pos_stoich_at_soc0));
    }
    ImGui::TextUnformatted("环境温度");
    ImGui::SameLine(label_col);
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::InputFloat("##temp", &ambient_temp_c_, 1.0f, 10.0f, "%.1f ℃");
    if (ImGui::IsItemHovered()) {
        if (model_type_ == 2 && thermal_type_ == 1) {
            ImGui::SetTooltip("集总热模型：Newton 冷却的环境温度 T_amb；"
                              "电芯初始温度也取该值");
        } else {
            ImGui::SetTooltip("等温假设：模型全程恒定温度（%g K）",
                              static_cast<double>(ambient_temp_c_) + 273.15);
        }
    }

    // 热模型选型（仅 DFN：SPM/SPMe 恒为等温，不提供耦合热选项）
    if (model_type_ == 2) {
        ImGui::TextUnformatted("热模型");
        ImGui::SameLine(label_col);
        const char* thermal_names[] = {"等温（恒温）", "集总热（lumped）"};
        ImGui::SetNextItemWidth(-1.0f);
        ImGui::Combo("##thermal", &thermal_type_, thermal_names,
                     IM_ARRAYSIZE(thermal_names));
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip(
                thermal_type_ == 1
                ? "单节点温度 ODE 与电化学双向耦合（参考 PyBaMM "
                  "options {\"thermal\": \"lumped\"}）：产热升温、"
                  "Newton 冷却降温，温度经 Arrhenius/BV 系数反馈动力学"
                : "等温假设：电芯温度恒为环境温度");
        }
    }

    // IDA 模式下时间步长为输出采样间隔，内部步长自适应
    ImGui::TextUnformatted(solver_type_ == 1 ? "输出步长" : "时间步长");
    ImGui::SameLine(label_col);
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::InputFloat("##dt", &time_step_, 0.1f, 1.0f, "%.2f s");
    ImGui::TextUnformatted("实时倍率");
    ImGui::SameLine(label_col);
    // 固定宽度输入框：为右侧提示留出可见空间（无需横向滚动）
    ImGui::SetNextItemWidth(100.0f);
    ImGui::InputFloat("##rtf", &real_time_factor_, 10.0f, 50.0f, "%.0fx");
    ImGui::SameLine();
    ImGui::TextDisabled("(0 = 不限速)");

    // 多段工况编辑（CC/CV/CP/静置；▲▼ 按钮调整段顺序）
    ImGui::Text("工况段");
    ImGui::SameLine();
    ImGui::TextDisabled("用 ▲▼ 调整段顺序；截止条件任一满足即进入下一段");

    int move_index = -1;
    int move_dir = 0;
    std::size_t remove_index = steps_.size();
    for (std::size_t i = 0; i < steps_.size(); ++i) {
        StepUI& step = steps_[i];
        ImGui::PushID(static_cast<int>(i));

        if (i > 0) {
            ImGui::Separator();
        }

        // 第一行：段号 + 模式/类型选择 + ▲▼ 上/下移 + 删除
        ImGui::Text("段 %d", static_cast<int>(i + 1));
        ImGui::SameLine();

        const char* mode_names[] = {"放电", "充电", "静置"};
        ImGui::SetNextItemWidth(80);
        ImGui::Combo("##mode", &step.mode, mode_names, IM_ARRAYSIZE(mode_names));

        if (step.mode == 2) {
            ImGui::SameLine();
            ImGui::TextDisabled("电流 = 0");
        } else {
            ImGui::SameLine();
            const char* control_names[] = {"恒流", "恒压", "恒功率"};
            ImGui::SetNextItemWidth(80);
            ImGui::Combo("##ctrl", &step.control, control_names,
                         IM_ARRAYSIZE(control_names));
        }
        ImGui::SameLine();
        ImGui::BeginDisabled(i == 0);
        if (ImGui::Button("\xE2\x96\xB2")) { // ▲ 上移
            move_index = static_cast<int>(i);
            move_dir = -1;
        }
        ImGui::EndDisabled();
        ImGui::SameLine();
        ImGui::BeginDisabled(i + 1 >= steps_.size());
        if (ImGui::Button("\xE2\x96\xBC")) { // ▼ 下移
            move_index = static_cast<int>(i);
            move_dir = 1;
        }
        ImGui::EndDisabled();
        ImGui::SameLine();
        if (ImGui::Button("删除")) {
            remove_index = i;
        }

        // 第二行：控制量与截止条件（随控制类型联动；label 在左、输入在右）
        if (step.mode == 2) {
            ImGui::TextUnformatted("截止时间");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(110);
            ImGui::InputFloat("##ct", &step.cutoff_time, 10.0f, 60.0f,
                              "%.0f s");
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("静置时长；仅按时间截止");
            }
        } else if (step.control == 1) {
            // 恒压：目标电压 + 时间/尾电流截止
            ImGui::TextUnformatted("目标电压");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(90);
            ImGui::InputFloat("##v", &step.voltage, 0.05f, 0.1f, "%.2f V");
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("本段端电压将跟踪该目标（充/放方向由电压差自然决定）");
            }
            ImGui::SameLine();
            ImGui::TextUnformatted("截止时间");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(100);
            ImGui::InputFloat("##ct", &step.cutoff_time, 10.0f, 60.0f,
                              "%.0f s");
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("本段最长时长；设为 0 表示不按时间截止");
            }
            ImGui::SameLine();
            ImGui::TextUnformatted("截止电流");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(100);
            ImGui::InputFloat("##ci", &step.cutoff_current, 0.05f, 0.1f,
                              "%.2f A");
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("尾电流：|I| 衰减到该值时段结束（CC-CV 充电惯例）；0 = 不启用");
            }
        } else {
            // 恒流 / 恒功率：控制量 + 时间/电压截止
            if (step.control == 0) {
                ImGui::TextUnformatted("电流[A]");
                ImGui::SameLine();
                ImGui::SetNextItemWidth(95);
                ImGui::InputFloat("##i", &step.current, 0.5f, 1.0f, "%.2f");
                ImGui::SameLine();
                ImGui::TextDisabled("%.2fC",
                                    params_.capacity > 0.0
                                        ? static_cast<double>(step.current)
                                              / params_.capacity
                                        : 0.0);
            } else {
                ImGui::TextUnformatted("功率[W]");
                ImGui::SameLine();
                ImGui::SetNextItemWidth(95);
                ImGui::InputFloat("##p", &step.power, 1.0f, 5.0f, "%.1f");
                ImGui::SameLine();
                ImGui::TextDisabled("\xE2\x89\x88 %.2fC",
                                    params_.capacity > 0.0
                                        ? static_cast<double>(step.power)
                                              / (3.7 * params_.capacity)
                                        : 0.0);
            }
            ImGui::SameLine();
            ImGui::TextUnformatted("截止时间");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(100);
            ImGui::InputFloat("##ct", &step.cutoff_time, 10.0f, 60.0f,
                              "%.0f s");
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("本段最长时长；设为 0 表示不按时间截止");
            }
            ImGui::SameLine();
            ImGui::TextUnformatted("截止电压");
            ImGui::SameLine();
            ImGui::SetNextItemWidth(100);
            ImGui::InputFloat("##cv", &step.cutoff_voltage, 0.05f, 0.1f,
                              "%.2f V");
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip(step.mode == 1 ? "充电上限电压；0 = 不启用"
                                                 : "放电下限电压；0 = 不启用");
            }
            // CC 段可选：电压截止后转恒压（CC-CV 一键工况，如 1C 充至
            // 4.4 V 后恒压至尾电流）；恒压目标 = 本段截止电压
            if (step.control == 0) {
                ImGui::Checkbox("截止后转恒压", &step.cccv);
                if (ImGui::IsItemHovered()) {
                    ImGui::SetTooltip(
                        "达到截止电压后自动转入恒压阶段（CC-CV）：\n"
                        "恒压目标 = 本段截止电压，尾电流/时间截止");
                }
                if (step.cccv) {
                    ImGui::SameLine();
                    ImGui::TextUnformatted("尾电流");
                    ImGui::SameLine();
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputFloat("##cccvi", &step.cutoff_current,
                                      0.05f, 0.1f, "%.2f A");
                    if (ImGui::IsItemHovered()) {
                        ImGui::SetTooltip(
                            "恒压段 |I| 衰减到该值时段结束；0 = 不按电流截止");
                    }
                }
            }
        }

        ImGui::PopID();
    }
    if (remove_index < steps_.size()) {
        steps_.erase(steps_.begin() + static_cast<std::ptrdiff_t>(remove_index));
    }
    // ▲▼ 重排在循环外执行，避免迭代中修改容器
    if (move_index >= 0) {
        const std::ptrdiff_t target =
            static_cast<std::ptrdiff_t>(move_index) + move_dir;
        if (target >= 0 &&
            target < static_cast<std::ptrdiff_t>(steps_.size())) {
            std::swap(steps_[static_cast<std::size_t>(move_index)],
                      steps_[static_cast<std::size_t>(target)]);
        }
    }

    if (ImGui::Button("+ 添加段")) {
        StepUI step; // 默认：静置 600 s，无电压截止
        step.mode = 2;
        step.current = 0.0f;
        step.cutoff_time = 600.0f;
        step.cutoff_voltage = 0.0f;
        steps_.push_back(step);
    }
    ImGui::SameLine();
    if (ImGui::Button("保存工况…")) {
        open_save_protocol_dialog();
    }
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip(
            "将当前段序列、循环圈数、模型/求解器选择与\n"
            "初始 SOC/环境温度存为 YAML 工况文件");
    }
    ImGui::SameLine();
    if (ImGui::Button("载入工况…")) {
        open_protocol_dialog();
    }
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("从 YAML 工况文件恢复完整仿真配置");
    }
    if (!protocol_load_error_.empty()) {
        ImGui::TextColored(ImVec4(0.9f, 0.3f, 0.3f, 1.0f),
                           "工况载入失败: %s",
                           protocol_load_error_.c_str());
    }
    if (!protocol_save_error_.empty()) {
        ImGui::TextColored(ImVec4(0.9f, 0.3f, 0.3f, 1.0f),
                           "工况保存失败: %s",
                           protocol_save_error_.c_str());
    }
    if (steps_.empty()) {
        ImGui::SameLine();
        ImGui::TextColored(ImVec4(0.9f, 0.3f, 0.3f, 1.0f), "至少需要一个工况段");
    }

    // 循环圈数（机理寿命仿真：段序列整体重复 N 圈）
    ImGui::TextUnformatted("循环圈数");
    ImGui::SameLine(label_col);
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::InputInt("##cycles", &cycles_);
    if (cycles_ < 1) {
        cycles_ = 1;
    }
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip(
            "上方段序列作为一个循环重复的圈数；>1 时每圈末快照\n"
            "（cycle_end）为寿命曲线采样点，降解累积反映容量衰减");
    }

    if (!protocol_error_.empty()) {
        ImGui::TextColored(ImVec4(0.9f, 0.3f, 0.3f, 1.0f), "工况错误: %s",
                           protocol_error_.c_str());
    }

    ImGui::TextDisabled("模型: SPM | 1C = %.1f A", params_.capacity);
}

void SimulationPanel::render_model_params_section() {
    if (!ImGui::CollapsingHeader("模型参数", ImGuiTreeNodeFlags_DefaultOpen)) {
        return;
    }
    ImGui::TextDisabled("（修改后需重新运行生效）");

    // SPMe/DFN：电解液传输参数（编辑统一写入 spme_params_，构建协议时同步到 DFN）
    if (model_type_ >= 1) {
        ImGui::SeparatorText("电解液");
        const float param_col = ImGui::GetCursorPosX() + 150.0f;
        float v;
        v = static_cast<float>(spme_params_.electrolyte_diffusivity);
        ImGui::TextUnformatted("D_e [m^2/s]");
        ImGui::SameLine(param_col);
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::InputFloat("##de", &v, 0.0f, 0.0f, "%.2e")) {
            spme_params_.electrolyte_diffusivity = std::max(1e-14f, v);
        }
        v = static_cast<float>(spme_params_.transference_number);
        ImGui::TextUnformatted("迁移数 t+");
        ImGui::SameLine(param_col);
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::InputFloat("##tn", &v, 0.01f, 0.1f, "%.3f")) {
            spme_params_.transference_number = std::clamp(v, 0.0f, 1.0f);
        }
        v = static_cast<float>(spme_params_.porosity_neg);
        ImGui::TextUnformatted("孔隙率（负极）");
        ImGui::SameLine(param_col);
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::InputFloat("##pn", &v, 0.01f, 0.1f, "%.3f")) {
            spme_params_.porosity_neg = std::clamp(v, 0.01f, 1.0f);
        }
        v = static_cast<float>(spme_params_.porosity_sep);
        ImGui::TextUnformatted("孔隙率（隔膜）");
        ImGui::SameLine(param_col);
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::InputFloat("##ps", &v, 0.01f, 0.1f, "%.3f")) {
            spme_params_.porosity_sep = std::clamp(v, 0.01f, 1.0f);
        }
        v = static_cast<float>(spme_params_.porosity_pos);
        ImGui::TextUnformatted("孔隙率（正极）");
        ImGui::SameLine(param_col);
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::InputFloat("##pp", &v, 0.01f, 0.1f, "%.3f")) {
            spme_params_.porosity_pos = std::clamp(v, 0.01f, 1.0f);
        }
        v = static_cast<float>(spme_params_.bruggeman_exponent);
        ImGui::TextUnformatted("Bruggeman 指数");
        ImGui::SameLine(param_col);
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::InputFloat("##br", &v, 0.05f, 0.5f, "%.2f")) {
            spme_params_.bruggeman_exponent = std::clamp(v, 0.0f, 5.0f);
        }
        v = static_cast<float>(spme_params_.electrolyte_conductivity);
        ImGui::TextUnformatted("电导率 κ [S/m]");
        ImGui::SameLine(param_col);
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::InputFloat("##k", &v, 0.01f, 0.1f, "%.3f")) {
            spme_params_.electrolyte_conductivity = std::max(0.001f, v);
        }
        v = static_cast<float>(spme_params_.c_e_init);
        ImGui::TextUnformatted("c_e 初始 [mol/m^3]");
        ImGui::SameLine(param_col);
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::InputFloat("##ce0", &v, 10.0f, 100.0f, "%.0f")) {
            spme_params_.c_e_init = std::max(1.0f, v);
        }
        ImGui::TextUnformatted("液相格数");
        ImGui::SameLine(param_col);
        ImGui::SetNextItemWidth(-1.0f);
        ImGui::InputInt("##cells", &spme_params_.electrolyte_cells);
        if (spme_params_.electrolyte_cells < 10) {
            spme_params_.electrolyte_cells = 10;
        }
    }

    if (model_type_ == 2) {
        ImGui::SeparatorText("固相电导 / 迭代控制");
        const float param_col = ImGui::GetCursorPosX() + 150.0f;
        float v;
        v = static_cast<float>(dfn_params_.neg_conductivity);
        ImGui::TextUnformatted("σ_n [S/m]");
        ImGui::SameLine(param_col);
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::InputFloat("##sn", &v, 0.0f, 0.0f, "%.3g")) {
            dfn_params_.neg_conductivity = std::max(1e-4f, v);
        }
        v = static_cast<float>(dfn_params_.pos_conductivity);
        ImGui::TextUnformatted("σ_p [S/m]");
        ImGui::SameLine(param_col);
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::InputFloat("##sp", &v, 0.0f, 0.0f, "%.3g")) {
            dfn_params_.pos_conductivity = std::max(1e-4f, v);
        }
        ImGui::TextUnformatted("最大迭代次数");
        ImGui::SameLine(param_col);
        ImGui::SetNextItemWidth(-1.0f);
        ImGui::InputInt("##maxit", &dfn_params_.max_iterations);
        if (dfn_params_.max_iterations < 1) {
            dfn_params_.max_iterations = 1;
        }
        v = static_cast<float>(dfn_params_.tolerance);
        ImGui::TextUnformatted("迭代容差");
        ImGui::SameLine(param_col);
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::InputFloat("##tol", &v, 0.0f, 0.0f, "%.1e")) {
            dfn_params_.tolerance = std::clamp(v, 1e-12f, 1e-2f);
        }
        v = static_cast<float>(dfn_params_.damping);
        ImGui::TextUnformatted("阻尼因子");
        ImGui::SameLine(param_col);
        ImGui::SetNextItemWidth(-1.0f);
        if (ImGui::InputFloat("##damp", &v, 0.05f, 0.1f, "%.2f")) {
            dfn_params_.damping = std::clamp(v, 0.01f, 1.0f);
        }

        // 热模型参数（仅 lumped 时生效；等温忽略）
        if (thermal_type_ == 1) {
            ImGui::SeparatorText("热模型（lumped）");
            const float th_col = ImGui::GetCursorPosX() + 150.0f;
            v = static_cast<float>(dfn_params_.thermal_params.h_total);
            ImGui::TextUnformatted("换热系数 h [W/(m²·K)]");
            ImGui::SameLine(th_col);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::InputFloat("##h", &v, 0.5f, 5.0f, "%.2f")) {
                dfn_params_.thermal_params.h_total = std::max(0.0f, v);
            }
            v = static_cast<float>(dfn_params_.thermal_params.cooling_area);
            ImGui::TextUnformatted("冷却面积 [m²]");
            ImGui::SameLine(th_col);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::InputFloat("##acool", &v, 0.0f, 0.0f, "%.5f")) {
                dfn_params_.thermal_params.cooling_area = std::max(0.0f, v);
            }
            v = static_cast<float>(dfn_params_.thermal_params.cell_volume);
            ImGui::TextUnformatted("电芯体积 V_cell [m³]");
            ImGui::SameLine(th_col);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::InputFloat("##vcell", &v, 0.0f, 0.0f, "%.3e")) {
                dfn_params_.thermal_params.cell_volume = std::max(1e-9f, v);
            }
            ImGui::TextDisabled("各层 ρcp / 集流体厚度 / 熵变系数在 YAML thermal 段配置");
        }

        // SEI 降解（机理寿命仿真；参数细节在 YAML degradation 段配置）
        ImGui::SeparatorText("SEI 降解（机理寿命）");
        ImGui::Checkbox("启用负极 SEI 生长",
                        &dfn_params_.degradation.enabled);
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip(
                "ec reaction limited SEI 生长 + 分布式膜电阻\n"
                "（PyBaMM options {\"SEI\": \"ec reaction limited\",\n"
                "\"SEI film resistance\": \"distributed\"} 口径）");
        }
        if (dfn_params_.degradation.enabled) {
            // 11 个 SEI 动力学参数（PyBaMM OKane2022 口径）：改动写入
            // dfn_params_ 缓存，start/单步立即生效，保存参数时写入
            // YAML degradation 段（round-trip 无损）
            auto& deg = dfn_params_.degradation;
            const float deg_col = ImGui::GetCursorPosX() + 170.0f;
            float v;

            v = static_cast<float>(deg.l_sei_0) * 1e9f;
            ImGui::TextUnformatted("初始厚度 [nm]");
            ImGui::SameLine(deg_col);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::InputFloat("##lsei0", &v, 0.0f, 0.0f, "%.3g")) {
                deg.l_sei_0 = std::max(0.0f, v) * 1e-9f;
            }
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip(
                    "初始 SEI 厚度 L_sei0（OKane2022 默认 5 nm）；"
                    "满电静置时也会持续生长");
            }

            v = static_cast<float>(deg.v_bar_sei);
            ImGui::TextUnformatted("偏摩尔体积 [m³/mol]");
            ImGui::SameLine(deg_col);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::InputFloat("##vbar", &v, 0.0f, 0.0f, "%.3e")) {
                deg.v_bar_sei = std::max(1e-10f, v);
            }
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip(
                    "SEI 偏摩尔体积：锂损失 = z/V_bar·a·ΔL·F 换算容量");
            }

            v = static_cast<float>(deg.z_sei);
            ImGui::TextUnformatted("锂摩尔比 z");
            ImGui::SameLine(deg_col);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::InputFloat("##zsei", &v, 0.0f, 0.0f, "%.2f")) {
                deg.z_sei = std::max(0.0f, v);
            }
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip(
                    "每生成 1 mol SEI 消耗的锂摩尔数（OKane2022 = 1）");
            }

            v = static_cast<float>(deg.u_sei);
            ImGui::TextUnformatted("平衡电位 U_sei [V]");
            ImGui::SameLine(deg_col);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::InputFloat("##usei", &v, 0.05f, 0.1f, "%.3f")) {
                deg.u_sei = std::clamp(v, 0.0f, 2.0f);
            }
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip(
                    "SEI 开路电位，决定 SEI 反应过电位 η_SEI = φ_s − φ_e − U_sei");
            }

            v = static_cast<float>(deg.r_sei);
            ImGui::TextUnformatted("电阻率 [Ω·m]");
            ImGui::SameLine(deg_col);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::InputFloat("##rsei", &v, 0.0f, 0.0f, "%.3e")) {
                deg.r_sei = std::max(0.0f, v);
            }
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip(
                    "SEI 电阻率 ρ：膜电阻随厚度增长（分布式），"
                    "压降 = i_surf·L·ρ；为 0 时仅生长不加阻抗");
            }

            v = static_cast<float>(deg.k_sei);
            ImGui::TextUnformatted("速率常数 k_sei [m/s]");
            ImGui::SameLine(deg_col);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::InputFloat("##ksei", &v, 0.0f, 0.0f, "%.3e")) {
                deg.k_sei = std::max(1e-30f, v);
            }
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip(
                    "EC 还原反应速率常数，控制生长速率的主参数");
            }

            v = static_cast<float>(deg.c_ec_0);
            ImGui::TextUnformatted("EC 初始浓度 [mol/m³]");
            ImGui::SameLine(deg_col);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::InputFloat("##cec0", &v, 10.0f, 100.0f, "%.0f")) {
                deg.c_ec_0 = std::max(1.0f, v);
            }
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("电解液 EC（碳酸乙烯酯）初始浓度 c_ec,0");
            }

            v = static_cast<float>(deg.d_ec);
            ImGui::TextUnformatted("EC 扩散系数 [m²/s]");
            ImGui::SameLine(deg_col);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::InputFloat("##dec", &v, 0.0f, 0.0f, "%.3e")) {
                deg.d_ec = std::max(1e-30f, v);
            }
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip(
                    "EC 在 SEI 中的扩散系数：越小生长越早转入扩散受限（1/(1+qL)）");
            }

            v = static_cast<float>(deg.e_sei);
            ImGui::TextUnformatted("活化能 [J/mol]");
            ImGui::SameLine(deg_col);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::InputFloat("##esei", &v, 100.0f, 1000.0f, "%.0f")) {
                deg.e_sei = std::max(0.0f, v);
            }
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip(
                    "SEI 生长 Arrhenius 活化能（等温模型下无温度反馈）");
            }

            v = static_cast<float>(deg.t_ref);
            ImGui::TextUnformatted("参考温度 [K]");
            ImGui::SameLine(deg_col);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::InputFloat("##tref", &v, 1.0f, 10.0f, "%.2f")) {
                deg.t_ref = std::max(1.0f, v);
            }
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Arrhenius 参考温度（等温时取环境温度）");
            }

            v = static_cast<float>(deg.alpha_sei);
            ImGui::TextUnformatted("传递系数 α");
            ImGui::SameLine(deg_col);
            ImGui::SetNextItemWidth(-1.0f);
            if (ImGui::InputFloat("##asei", &v, 0.05f, 0.1f, "%.2f")) {
                deg.alpha_sei = std::clamp(v, 0.01f, 0.99f);
            }
            if (ImGui::IsItemHovered()) {
                ImGui::SetTooltip(
                    "SEI 反应传递系数（对称 BV，0 < α < 1）");
            }
        }
    }
}

void SimulationPanel::render_cell_params() {
    if (!ImGui::CollapsingHeader("电芯参数（YAML 只读）",
                                 ImGuiTreeNodeFlags_DefaultOpen)) {
        return;
    }
    ImGui::TextDisabled("来源: %s", params_path_);

    // 两列布局：参数名 | 值。展示 YAML 加载的固相/几何参数（三模型共享），
    // 厚度/半径以 µm 显示，扩散系数/交换电流用科学计数法
    if (!ImGui::BeginTable("cell_params", 2,
                           ImGuiTableFlags_SizingStretchProp)) {
        return;
    }

    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::SeparatorText("几何");
    cell_param_row("负极厚度 L_n", "%.1f µm",
                   static_cast<double>(params_.neg_thickness) * 1e6);
    cell_param_row("隔膜厚度", "%.1f µm",
                   static_cast<double>(params_.sep_thickness) * 1e6);
    cell_param_row("正极厚度 L_p", "%.1f µm",
                   static_cast<double>(params_.pos_thickness) * 1e6);
    cell_param_row2("极片尺寸 (H×W)", "%.3f × %.2f m",
                    static_cast<double>(params_.electrode_height),
                    static_cast<double>(params_.electrode_width));
    cell_param_row("负极颗粒半径 R_n", "%.2f µm",
                   static_cast<double>(params_.neg_particle_radius) * 1e6);
    cell_param_row("正极颗粒半径 R_p", "%.2f µm",
                   static_cast<double>(params_.pos_particle_radius) * 1e6);

    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::SeparatorText("负极（石墨）");
    cell_param_row("活性材料分数 ε_s,n", "%.3f",
                   static_cast<double>(params_.neg_active_fraction));
    cell_param_row("固相扩散系数 D_s,n", "%.2e m²/s",
                   static_cast<double>(params_.neg_diffusion));
    cell_param_row("最大锂浓度 c_max,n", "%.0f mol/m³",
                   static_cast<double>(params_.neg_max_conc));
    cell_param_row("交换电流参考 m_ref,n", "%.2e",
                   static_cast<double>(params_.neg_m_ref));
    cell_param_row("活化能 E_r,n", "%.0f J/mol",
                   static_cast<double>(params_.neg_activation_energy));

    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::SeparatorText("正极（NMC）");
    cell_param_row("活性材料分数 ε_s,p", "%.3f",
                   static_cast<double>(params_.pos_active_fraction));
    cell_param_row("固相扩散系数 D_s,p", "%.2e m²/s",
                   static_cast<double>(params_.pos_diffusion));
    cell_param_row("最大锂浓度 c_max,p", "%.0f mol/m³",
                   static_cast<double>(params_.pos_max_conc));
    cell_param_row("交换电流参考 m_ref,p", "%.2e",
                   static_cast<double>(params_.pos_m_ref));
    cell_param_row("活化能 E_r,p", "%.0f J/mol",
                   static_cast<double>(params_.pos_activation_energy));

    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::SeparatorText("容量 / 化学计量");
    cell_param_row("名义容量", "%.2f Ah",
                   static_cast<double>(params_.capacity));
    cell_param_row2("负极嵌入率 x @0%/@100%", "%.4f / %.4f",
                    static_cast<double>(params_.neg_stoich_at_soc0),
                    static_cast<double>(params_.neg_stoich_at_soc1));
    cell_param_row2("正极嵌入率 y @0%/@100%", "%.4f / %.4f",
                    static_cast<double>(params_.pos_stoich_at_soc0),
                    static_cast<double>(params_.pos_stoich_at_soc1));
    cell_param_row("颗粒径向格数", "%.0f",
                   static_cast<double>(params_.particle_cells));

    // 热模型参数（DFN 专属；YAML thermal 段加载值）
    if (model_type_ == 2) {
        const auto& th = dfn_params_.thermal_params;
        ImGui::TableNextRow();
        ImGui::TableNextColumn();
        ImGui::SeparatorText("热模型（DFN）");
        cell_param_row2("环境/初始温度", "%.2f / %.2f K",
                        static_cast<double>(th.t_ambient),
                        static_cast<double>(th.t_initial));
        cell_param_row("换热系数 h", "%.2f W/(m²·K)",
                       static_cast<double>(th.h_total));
        cell_param_row("冷却面积 A_cooling", "%.5f m²",
                       static_cast<double>(th.cooling_area));
        cell_param_row("电芯体积 V_cell", "%.3e m³",
                       static_cast<double>(th.cell_volume));
        cell_param_row3("ρcp 负极/隔膜/正极", "%.0f / %.0f / %.0f",
                        static_cast<double>(th.rho_c_p_neg),
                        static_cast<double>(th.rho_c_p_sep),
                        static_cast<double>(th.rho_c_p_pos));
        cell_param_row2("ρcp 集流体 负/正", "%.0f / %.0f",
                        static_cast<double>(th.rho_c_p_cc_neg),
                        static_cast<double>(th.rho_c_p_cc_pos));
        cell_param_row2("集流体厚度 负/正", "%.1f / %.1f µm",
                        static_cast<double>(th.cc_thickness_neg) * 1e6,
                        static_cast<double>(th.cc_thickness_pos) * 1e6);
        // 熵变系数 dU/dT：LUT 显示点数，否则显示常数 [V/K]
        if (!th.entropic_lut_neg.empty()) {
            cell_param_row("熵变系数 dU/dT 负极", "LUT（%.0f 点）",
                           static_cast<double>(th.entropic_lut_neg.size()));
        } else {
            cell_param_row("熵变系数 dU/dT 负极", "%.2e V/K",
                           static_cast<double>(th.entropic_change_neg));
        }
        if (!th.entropic_lut_pos.empty()) {
            cell_param_row("熵变系数 dU/dT 正极", "LUT（%.0f 点）",
                           static_cast<double>(th.entropic_lut_pos.size()));
        } else {
            cell_param_row("熵变系数 dU/dT 正极", "%.2e V/K",
                           static_cast<double>(th.entropic_change_pos));
        }
    }

    ImGui::EndTable();

    // ── 电芯 OCV 曲线（vs SOC）──
    // x(SOC)、y(SOC) 沿参数的化学计量端点线性映射，U_ocv = U_p − U_n；
    // OCP 数据优先取参数的 YAML "ocp:" 段，未提供时用内置 Chen2020 公式
    constexpr int kOcvSamples = 120;
    std::vector<real_t> soc_axis(kOcvSamples);
    std::vector<real_t> u_neg(kOcvSamples);
    std::vector<real_t> u_pos(kOcvSamples);
    std::vector<real_t> u_ocv(kOcvSamples);
    for (int i = 0; i < kOcvSamples; ++i) {
        const real_t s = static_cast<real_t>(i) / (kOcvSamples - 1);
        soc_axis[static_cast<std::size_t>(i)] = s * 100.0f;
        u_neg[static_cast<std::size_t>(i)] = params_.neg_ocp_at(
            params_.neg_stoich_at_soc0
            + (params_.neg_stoich_at_soc1 - params_.neg_stoich_at_soc0) * s);
        u_pos[static_cast<std::size_t>(i)] = params_.pos_ocp_at(
            params_.pos_stoich_at_soc0
            + (params_.pos_stoich_at_soc1 - params_.pos_stoich_at_soc0) * s);
        u_ocv[static_cast<std::size_t>(i)] =
            u_pos[static_cast<std::size_t>(i)]
            - u_neg[static_cast<std::size_t>(i)];
    }
    if (ImPlot::BeginPlot("##cell_ocv", ImVec2(-1, 200))) {
        ImPlot::SetupAxes("SOC [%]", "U [V]",
                          ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_None);
        fit_y_padded(u_neg, u_pos, u_ocv);
        ImPlot::SetupFinish();
        ImPlot::PlotLine("负极 U_n(x)", soc_axis.data(), u_neg.data(),
                         kOcvSamples);
        ImPlot::PlotLine("正极 U_p(y)", soc_axis.data(), u_pos.data(),
                         kOcvSamples);
        ImPlot::SetNextLineStyle(IMPLOT_AUTO_COL, 2.5f);
        ImPlot::PlotLine("电芯 OCV", soc_axis.data(), u_ocv.data(),
                         kOcvSamples);
        ImPlot::EndPlot();
    }
    ImGui::TextDisabled(
        "OCP 来源: 负极 %s / 正极 %s",
        params_.neg_ocp_stoich.empty() ? "内置 Chen2020 公式" : "YAML ocp 段查表",
        params_.pos_ocp_stoich.empty() ? "内置 Chen2020 公式" : "YAML ocp 段查表");
}

void SimulationPanel::render_progress_section() {
    // 进度已移至主窗口状态栏（迷你进度条 + 百分比）；
    // 失败原因较长，仍留在面板内完整显示
    if (engine_.status() == SimulationStatus::Failed) {
        ImGui::TextColored(ImVec4(0.9f, 0.3f, 0.3f, 1.0f), "失败:");
        // 只读多行文本框呈现：文本可鼠标选中、Ctrl+C 复制，便于把完整
        // 错误（含失败时刻状态）贴进工单/日志；新失败发生时同步缓冲，
        // 其余帧指针稳定（ImGui std::string 配方的帧间要求）
        const std::string& msg = engine_.error_message();
        if (error_select_buf_ != msg) {
            error_select_buf_ = msg;
        }
        ImGui::InputTextMultiline(
            "##failure_msg", error_select_buf_.data(),
            error_select_buf_.capacity() + 1,
            ImVec2(0.0f, ImGui::GetTextLineHeight() * 4.0f),
            ImGuiInputTextFlags_ReadOnly | ImGuiInputTextFlags_CallbackResize,
            input_text_resize_callback, &error_select_buf_);
    }
}

void SimulationPanel::start() {
    // 先本地校验：非法工况不应终止正在进行的旧仿真
    const MultiStepProtocol protocol = build_protocol();
    try {
        protocol.validate();
    } catch (const std::exception& e) {
        protocol_error_ = e.what();
        return;
    }
    protocol_error_.clear();

    // 重新计算前通知主窗口清理上一次的曲线数据
    if (on_restart_) {
        on_restart_();
    }
    engine_.start(effective_params(), protocol);
}

void SimulationPanel::step_once() {
    const MultiStepProtocol protocol = build_protocol();
    try {
        protocol.validate();
    } catch (const std::exception& e) {
        protocol_error_ = e.what();
        return;
    }
    protocol_error_.clear();

    // 从空闲首次单步：引擎会重新初始化，同样需要清理旧曲线
    if (engine_.status() == SimulationStatus::Idle && on_restart_) {
        on_restart_();
    }
    engine_.single_step(effective_params(), protocol);
}

model::SPMModel::Parameters SimulationPanel::effective_params() const {
    // 参数文件值为基础，叠加面板仿真条件（SOC 百分比 → 0~1，摄氏 → 开尔文）
    auto params = params_;
    params.initial_soc = static_cast<real_t>(initial_soc_pct_) / 100.0;
    params.temperature = static_cast<real_t>(ambient_temp_c_) + 273.15;
    return params;
}

MultiStepProtocol SimulationPanel::build_protocol() const {
    MultiStepProtocol proto;
    proto.time_step = static_cast<real_t>(time_step_);
    proto.real_time_factor = static_cast<real_t>(real_time_factor_);
    proto.model = static_cast<CellModelType>(model_type_);
    proto.solver = static_cast<SolverType>(solver_type_);
    proto.solver_rel_tol = static_cast<real_t>(solver_rel_tol_);
    // SPMe 参数：YAML 电解液部分 + 面板运行条件（SOC/温度）覆盖固相
    proto.spme_params = spme_params_;
    proto.spme_params.solid = effective_params();
    // DFN 参数：同 SPMe 固相覆盖，保留 YAML 固相电导与迭代控制；
    // 电解液字段与面板编辑值同步（两模型共享同一组电解液 UI）
    proto.dfn_params = dfn_params_;
    proto.dfn_params.solid = effective_params();
    proto.dfn_params.electrolyte_diffusivity = proto.spme_params.electrolyte_diffusivity;
    proto.dfn_params.transference_number = proto.spme_params.transference_number;
    proto.dfn_params.porosity_neg = proto.spme_params.porosity_neg;
    proto.dfn_params.porosity_sep = proto.spme_params.porosity_sep;
    proto.dfn_params.porosity_pos = proto.spme_params.porosity_pos;
    proto.dfn_params.bruggeman_exponent = proto.spme_params.bruggeman_exponent;
    proto.dfn_params.electrolyte_conductivity = proto.spme_params.electrolyte_conductivity;
    proto.dfn_params.c_e_init = proto.spme_params.c_e_init;
    proto.dfn_params.electrolyte_cells = proto.spme_params.electrolyte_cells;
    // 热模型（DFN 专属）：UI 选型 + 面板环境温度联动（等温模式不使用热参数）
    proto.dfn_params.thermal_type =
        thermal_type_ == 1 ? model::ThermalType::Lumped
                           : model::ThermalType::Isothermal;
    const real_t t_amb = static_cast<real_t>(ambient_temp_c_) + 273.15;
    proto.dfn_params.thermal_params.t_ambient = t_amb;
    // GUI 未提供单独的初始电芯温度输入：lumped 从环境温度起步
    //（与等温模式“温度 = 面板环境温度”的语义一致）
    proto.dfn_params.thermal_params.t_initial = t_amb;
    // 循环圈数：段序列整体重复 N 次（机理寿命仿真工况）
    proto.cycles = cycles_;
    for (const StepUI& ui : steps_) {
        CCStep step;
        step.mode = static_cast<CCStep::Mode>(ui.mode);
        // 静置段固定为电流控制（电流 0），CV/CP 无意义
        step.control = ui.mode == 2
                         ? CCStep::Control::Current
                         : static_cast<CCStep::Control>(ui.control);
        step.current = static_cast<real_t>(ui.current);
        step.voltage = static_cast<real_t>(ui.voltage);
        step.power = static_cast<real_t>(ui.power);
        step.cutoff_time = static_cast<real_t>(ui.cutoff_time);
        step.cutoff_voltage = static_cast<real_t>(ui.cutoff_voltage);
        step.cutoff_current = static_cast<real_t>(ui.cutoff_current);
        proto.steps.push_back(step);
        // CC-CV 一键工况：电压截止后转恒压（截止电压即恒压目标），
        // 展开为 CC + CV 两段；恒压段截止时间沿用本段上限
        if (ui.cccv && ui.mode != 2 && ui.control == 0
            && ui.cutoff_voltage > 0.0f) {
            CCStep cv = step;
            cv.control = CCStep::Control::Voltage;
            cv.voltage = step.cutoff_voltage;
            cv.cutoff_voltage = 0.0f;
            cv.cutoff_current = static_cast<real_t>(ui.cutoff_current);
            proto.steps.push_back(cv);
        }
    }
    return proto;
}

void SimulationPanel::pause_resume() {
    if (engine_.status() == SimulationStatus::Paused) {
        engine_.resume();
    } else {
        engine_.pause();
    }
}

void SimulationPanel::stop() {
    engine_.stop();
}

void SimulationPanel::reset() {
    // 重置同样清理曲线，与 start() 一致
    if (on_restart_) {
        on_restart_();
    }
    engine_.stop();
}

} // namespace fluxcell::gui
