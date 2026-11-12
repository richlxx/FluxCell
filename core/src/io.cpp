/// @file io.cpp
/// @brief YAML 参数读取 + HDF5/CSV 结果写入
#include <fluxcell/core/io/io.hpp>

#include <yaml-cpp/yaml.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fluxcell::io {

// ══════════════════════════════════════════════════════════════
//  YAML 参数加载
// ══════════════════════════════════════════════════════════════

/// 严格解析整个字符串为实数（stod 会接受数字前缀，需校验完整消耗）
bool parse_fully(const std::string& text, real_t& out) {
    if (text.empty()) return false;
    try {
        std::size_t consumed = 0;
        out = std::stod(text, &consumed);
        return consumed == text.size();
    } catch (const std::exception&) {
        return false;
    }
}

model::EntropyLut load_entropy_lut_csv(const std::filesystem::path& csv_path) {
    std::ifstream file(csv_path);
    if (!file) {
        throw std::runtime_error("无法打开熵变系数文件: " + csv_path.string());
    }
    model::EntropyLut lut;
    std::string line;
    while (std::getline(file, line)) {
        // 两列 CSV：嵌锂量 sto、dU/dT [mV/K]；表头/空行/解析失败行跳过
        const auto comma = line.find(',');
        if (comma == std::string::npos) continue;
        real_t sto = 0.0;
        real_t dudt_mv = 0.0;
        if (!parse_fully(line.substr(0, comma), sto)
                || !parse_fully(line.substr(comma + 1), dudt_mv)) {
            continue;
        }
        lut.sto.push_back(sto);
        lut.dudt.push_back(dudt_mv * 1e-3); // mV/K → V/K
    }
    if (lut.size() < 2) {
        throw std::runtime_error("熵变系数文件有效数据点不足（需 ≥ 2 点）: "
                                 + csv_path.string());
    }
    // 按嵌锂量升序排列（插值节点必须严格递增）
    std::vector<std::size_t> order(lut.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(),
              [&](std::size_t a, std::size_t b) { return lut.sto[a] < lut.sto[b]; });
    std::vector<real_t> sto(lut.size());
    std::vector<real_t> dudt(lut.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        sto[i] = lut.sto[order[i]];
        dudt[i] = lut.dudt[order[i]];
    }
    lut.sto = std::move(sto);
    lut.dudt = std::move(dudt);
    return lut;
}

/// 解析 thermal.entropic_change_* 节点：标量常数 [V/K]、CSV 文件路径
/// （相对 yml 所在目录解析）或内联映射 {sto: [...], dudt: [...]}
void parse_entropy_parameter(const YAML::Node& node,
                             const std::filesystem::path& base_dir,
                             real_t& scalar_out, model::EntropyLut& lut_out) {
    if (node.IsMap()) {
        auto sto = node["sto"].as<std::vector<real_t>>();
        auto dudt = node["dudt"].as<std::vector<real_t>>();
        if (sto.size() != dudt.size() || sto.size() < 2) {
            throw std::invalid_argument(
                "entropic_change: sto/dudt 序列长度需一致且 ≥ 2 点");
        }
        lut_out.sto = std::move(sto);
        lut_out.dudt = std::move(dudt);
        return;
    }
    if (!node.IsScalar()) {
        throw std::invalid_argument(
            "entropic_change: 仅支持标量常数 / CSV 路径字符串 / 内联映射");
    }
    const std::string text = node.as<std::string>();
    real_t value = 0.0;
    if (parse_fully(text, value)) {
        scalar_out = value; // 数值标量 → 常数 dU/dT
        return;
    }
    std::filesystem::path csv_path(text); // 非数值 → CSV 文件路径
    if (csv_path.is_relative()) {
        csv_path = base_dir / csv_path;
    }
    lut_out = load_entropy_lut_csv(csv_path);
}

model::DFNModel::Parameters load_dfn_parameters(const std::filesystem::path& yaml_path) {
    model::DFNModel::Parameters p;
    // 固相 + 电解液部分复用 SPMe 读取（DFN 参数继承 SPMe 参数）
    static_cast<model::SPMeModel::Parameters&>(p) = load_spme_parameters(yaml_path);

    YAML::Node root = YAML::LoadFile(yaml_path.string());

    // DFN 特有键：固相电导率
    if (auto neg = root["negative_electrode"]) {
        if (neg["conductivity"])
            p.neg_conductivity = neg["conductivity"].as<real_t>();
    }
    if (auto pos = root["positive_electrode"]) {
        if (pos["conductivity"])
            p.pos_conductivity = pos["conductivity"].as<real_t>();
    }

    // 热模型（DFN 专属；缺省等温，参数缺省值即 PyBaMM Chen2020 口径）
    if (auto th = root["thermal"]) {
        if (th["type"]) {
            const std::string type = th["type"].as<std::string>();
            if (type == "lumped") {
                p.thermal_type = model::ThermalType::Lumped;
            } else if (type == "isothermal") {
                p.thermal_type = model::ThermalType::Isothermal;
            } else {
                throw std::invalid_argument(
                    "thermal.type: 未知热模型类型 '" + type
                    + "'（可选 isothermal/lumped）");
            }
        }
        if (th["ambient_temperature"])
            p.thermal_params.t_ambient = th["ambient_temperature"].as<real_t>();
        if (th["initial_temperature"])
            p.thermal_params.t_initial = th["initial_temperature"].as<real_t>();
        if (th["heat_transfer_coefficient"])
            p.thermal_params.h_total
                = th["heat_transfer_coefficient"].as<real_t>();
        if (th["cooling_area"])
            p.thermal_params.cooling_area = th["cooling_area"].as<real_t>();
        if (th["cell_volume"])
            p.thermal_params.cell_volume = th["cell_volume"].as<real_t>();
        // 各层体积比热容 ρ·c_p [J/(m³·K)]（cc = current collector 集流体）
        if (th["rho_cp_negative"])
            p.thermal_params.rho_c_p_neg = th["rho_cp_negative"].as<real_t>();
        if (th["rho_cp_separator"])
            p.thermal_params.rho_c_p_sep = th["rho_cp_separator"].as<real_t>();
        if (th["rho_cp_positive"])
            p.thermal_params.rho_c_p_pos = th["rho_cp_positive"].as<real_t>();
        if (th["rho_cp_cc_negative"])
            p.thermal_params.rho_c_p_cc_neg
                = th["rho_cp_cc_negative"].as<real_t>();
        if (th["rho_cp_cc_positive"])
            p.thermal_params.rho_c_p_cc_pos
                = th["rho_cp_cc_positive"].as<real_t>();
        if (th["cc_thickness_negative"])
            p.thermal_params.cc_thickness_neg
                = th["cc_thickness_negative"].as<real_t>();
        if (th["cc_thickness_positive"])
            p.thermal_params.cc_thickness_pos
                = th["cc_thickness_positive"].as<real_t>();
        // 熵变系数：标量常数 [V/K]、CSV 路径（相对本文件）或内联映射
        if (th["entropic_change_negative"]) {
            parse_entropy_parameter(th["entropic_change_negative"],
                                    yaml_path.parent_path(),
                                    p.thermal_params.entropic_change_neg,
                                    p.thermal_params.entropic_lut_neg);
        }
        if (th["entropic_change_positive"]) {
            parse_entropy_parameter(th["entropic_change_positive"],
                                    yaml_path.parent_path(),
                                    p.thermal_params.entropic_change_pos,
                                    p.thermal_params.entropic_lut_pos);
        }
    }

    // SEI 降解（机理寿命仿真；缺省关闭，参数缺省值 = PyBaMM OKane2022）
    if (auto dg = root["degradation"]) {
        if (dg["enabled"])
            p.degradation.enabled = dg["enabled"].as<bool>();
        if (auto sei = dg["sei"]) {
            if (sei["initial_thickness"])
                p.degradation.l_sei_0 = sei["initial_thickness"].as<real_t>();
            if (sei["molar_volume"])
                p.degradation.v_bar_sei = sei["molar_volume"].as<real_t>();
            if (sei["lithium_moles_ratio"])
                p.degradation.z_sei = sei["lithium_moles_ratio"].as<real_t>();
            if (sei["u_sei"])
                p.degradation.u_sei = sei["u_sei"].as<real_t>();
            if (sei["resistivity"])
                p.degradation.r_sei = sei["resistivity"].as<real_t>();
            if (sei["rate_constant"])
                p.degradation.k_sei = sei["rate_constant"].as<real_t>();
            if (sei["c_ec_0"])
                p.degradation.c_ec_0 = sei["c_ec_0"].as<real_t>();
            if (sei["diffusivity"])
                p.degradation.d_ec = sei["diffusivity"].as<real_t>();
            if (sei["activation_energy"])
                p.degradation.e_sei = sei["activation_energy"].as<real_t>();
            if (sei["reference_temperature"])
                p.degradation.t_ref = sei["reference_temperature"].as<real_t>();
            if (sei["alpha"])
                p.degradation.alpha_sei = sei["alpha"].as<real_t>();
        }
    }

    return p;
}

CellGeometry load_cell_geometry(const std::filesystem::path& yaml_path) {
    YAML::Node root = YAML::LoadFile(yaml_path.string());
    CellGeometry geom;

    if (auto g = root["geometry"]) {
        if (g["negative_thickness"]) geom.negative_thickness = g["negative_thickness"].as<real_t>();
        if (g["separator_thickness"]) geom.separator_thickness = g["separator_thickness"].as<real_t>();
        if (g["positive_thickness"]) geom.positive_thickness = g["positive_thickness"].as<real_t>();
        if (g["electrode_height"]) geom.electrode_height = g["electrode_height"].as<real_t>();
        if (g["electrode_width"]) geom.electrode_width = g["electrode_width"].as<real_t>();
        if (g["negative_particle_radius"]) geom.negative_particle_radius = g["negative_particle_radius"].as<real_t>();
        if (g["positive_particle_radius"]) geom.positive_particle_radius = g["positive_particle_radius"].as<real_t>();
        if (g["negative_mesh_cells"]) geom.negative_mesh_cells = g["negative_mesh_cells"].as<int>();
        if (g["separator_mesh_cells"]) geom.separator_mesh_cells = g["separator_mesh_cells"].as<int>();
        if (g["positive_mesh_cells"]) geom.positive_mesh_cells = g["positive_mesh_cells"].as<int>();
        if (g["particle_mesh_cells"]) geom.particle_mesh_cells = g["particle_mesh_cells"].as<int>();
    }

    return geom;
}

model::SPMModel::Parameters load_spm_parameters(const std::filesystem::path& yaml_path) {
    YAML::Node root = YAML::LoadFile(yaml_path.string());
    model::SPMModel::Parameters p;

    if (auto neg = root["negative_electrode"]) {
        if (neg["thickness"]) p.neg_thickness = neg["thickness"].as<real_t>();
        if (neg["particle_radius"]) p.neg_particle_radius = neg["particle_radius"].as<real_t>();
        if (neg["active_fraction"]) p.neg_active_fraction = neg["active_fraction"].as<real_t>();
        if (neg["diffusion_coefficient"]) p.neg_diffusion = neg["diffusion_coefficient"].as<real_t>();
        if (neg["max_concentration"]) p.neg_max_conc = neg["max_concentration"].as<real_t>();
        if (neg["exchange_current_density"]) p.neg_m_ref = neg["exchange_current_density"].as<real_t>();
        if (neg["activation_energy"]) p.neg_activation_energy = neg["activation_energy"].as<real_t>();
    }

    if (auto pos = root["positive_electrode"]) {
        if (pos["thickness"]) p.pos_thickness = pos["thickness"].as<real_t>();
        if (pos["particle_radius"]) p.pos_particle_radius = pos["particle_radius"].as<real_t>();
        if (pos["active_fraction"]) p.pos_active_fraction = pos["active_fraction"].as<real_t>();
        if (pos["diffusion_coefficient"]) p.pos_diffusion = pos["diffusion_coefficient"].as<real_t>();
        if (pos["max_concentration"]) p.pos_max_conc = pos["max_concentration"].as<real_t>();
        if (pos["exchange_current_density"]) p.pos_m_ref = pos["exchange_current_density"].as<real_t>();
        if (pos["activation_energy"]) p.pos_activation_energy = pos["activation_energy"].as<real_t>();
    }

    if (auto sep = root["separator"]) {
        if (sep["thickness"]) p.sep_thickness = sep["thickness"].as<real_t>();
    }

    if (auto ely = root["electrolyte"]) {
        if (ely["concentration"]) p.electrolyte_conc = ely["concentration"].as<real_t>();
    }

    if (auto g = root["geometry"]) {
        if (g["electrode_height"]) p.electrode_height = g["electrode_height"].as<real_t>();
        if (g["electrode_width"]) p.electrode_width = g["electrode_width"].as<real_t>();
        if (g["particle_mesh_cells"]) p.particle_cells = g["particle_mesh_cells"].as<int>();
    }

    if (auto kin = root["kinetics"]) {
        if (kin["temperature"]) p.temperature = kin["temperature"].as<real_t>();
    }

    if (auto cell = root["cell"]) {
        if (cell["capacity"]) p.capacity = cell["capacity"].as<real_t>();
        if (cell["initial_soc"]) p.initial_soc = cell["initial_soc"].as<real_t>();
        // SOC → 嵌入率线性映射端点（缺省用内置 Chen2020 值）
        if (cell["neg_stoich_at_soc0"])
            p.neg_stoich_at_soc0 = cell["neg_stoich_at_soc0"].as<real_t>();
        if (cell["neg_stoich_at_soc1"])
            p.neg_stoich_at_soc1 = cell["neg_stoich_at_soc1"].as<real_t>();
        if (cell["pos_stoich_at_soc0"])
            p.pos_stoich_at_soc0 = cell["pos_stoich_at_soc0"].as<real_t>();
        if (cell["pos_stoich_at_soc1"])
            p.pos_stoich_at_soc1 = cell["pos_stoich_at_soc1"].as<real_t>();
    }

    // OCP 曲线（可选）：stoich/voltage 等长采样序列（升序），线性插值；
    // 未提供时模型回退到内置 Chen2020 公式
    if (auto ocp = root["ocp"]) {
        const auto read_ocp_lut = [](const YAML::Node& node,
                                     const char* name,
                                     std::vector<real_t>& stoich,
                                     std::vector<real_t>& voltage) {
            if (node["stoich"])
                stoich = node["stoich"].as<std::vector<real_t>>();
            if (node["voltage"])
                voltage = node["voltage"].as<std::vector<real_t>>();
            if (stoich.size() != voltage.size()) {
                throw std::invalid_argument(
                    std::string("ocp.") + name
                    + ": stoich/voltage 长度不一致");
            }
        };
        if (auto neg = ocp["negative_electrode"]) {
            read_ocp_lut(neg, "negative_electrode",
                         p.neg_ocp_stoich, p.neg_ocp_voltage);
        }
        if (auto pos = ocp["positive_electrode"]) {
            read_ocp_lut(pos, "positive_electrode",
                         p.pos_ocp_stoich, p.pos_ocp_voltage);
        }
    }

    return p;
}

model::SPMeModel::Parameters load_spme_parameters(const std::filesystem::path& yaml_path) {
    model::SPMeModel::Parameters p;
    p.solid = load_spm_parameters(yaml_path); // 固相部分复用 SPM 读取

    YAML::Node root = YAML::LoadFile(yaml_path.string());

    if (auto ely = root["electrolyte"]) {
        if (ely["concentration"])
            p.solid.electrolyte_conc = ely["concentration"].as<real_t>();
        if (ely["initial_concentration"])
            p.c_e_init = ely["initial_concentration"].as<real_t>();
        if (ely["diffusivity"])
            p.electrolyte_diffusivity = ely["diffusivity"].as<real_t>();
        if (ely["transference_number"])
            p.transference_number = ely["transference_number"].as<real_t>();
        if (ely["porosity_neg"])
            p.porosity_neg = ely["porosity_neg"].as<real_t>();
        if (ely["porosity_sep"])
            p.porosity_sep = ely["porosity_sep"].as<real_t>();
        if (ely["porosity_pos"])
            p.porosity_pos = ely["porosity_pos"].as<real_t>();
        if (ely["bruggeman_exponent"])
            p.bruggeman_exponent = ely["bruggeman_exponent"].as<real_t>();
        if (ely["conductivity"])
            p.electrolyte_conductivity = ely["conductivity"].as<real_t>();

        // 可选：浓度相关输运表（标定曲线如 Nyman2008；节点升序、线性插值、
        // 端点钳位），加载后模型按浓度逐格求值，替代单点常量近似
        if (auto tbl = ely["transport_table"]) {
            if (tbl["nodes"])
                p.electrolyte_conc_lut = tbl["nodes"].as<std::vector<real_t>>();
            if (tbl["conductivity"])
                p.conductivity_lut = tbl["conductivity"].as<std::vector<real_t>>();
            if (tbl["diffusivity"])
                p.diffusivity_lut = tbl["diffusivity"].as<std::vector<real_t>>();
            if (p.conductivity_lut.empty() && p.diffusivity_lut.empty()) {
                throw std::invalid_argument(
                    "electrolyte.transport_table: 需至少提供 conductivity 或 diffusivity 序列");
            }
            if ((!p.conductivity_lut.empty()
                 && p.conductivity_lut.size() != p.electrolyte_conc_lut.size())
                || (!p.diffusivity_lut.empty()
                    && p.diffusivity_lut.size() != p.electrolyte_conc_lut.size())) {
                throw std::invalid_argument(
                    "electrolyte.transport_table: nodes 与值序列长度不一致");
            }
        }
    }

    if (auto g = root["geometry"]) {
        if (g["electrolyte_mesh_cells"])
            p.electrolyte_cells = g["electrolyte_mesh_cells"].as<int>();
    }

    return p;
}

// ══════════════════════════════════════════════════════════════
//  YAML 参数保存（与上方读取键结构对称）
// ══════════════════════════════════════════════════════════════

namespace {

/// 以 DFN 参数为超集写出 YAML 映射体（不含顶层 BeginMap/EndMap）：
/// with_electrolyte 控制是否写 electrolyte 扩展键（SPMe/DFN），
/// with_conductivity 控制是否写固相电导键（DFN）；
/// SPM/SPMe 保存时先聚合到默认 DFN 参数上再调用，避免重复分支
void emit_parameters(YAML::Emitter& out,
                     const model::DFNModel::Parameters& p,
                     bool with_electrolyte,
                     bool with_conductivity) {
    const auto& s = p.solid;
    const auto emit_seq = [&out](const std::vector<real_t>& v) {
        out << YAML::BeginSeq;
        for (const real_t x : v) {
            out << x;
        }
        out << YAML::EndSeq;
    };

    out << YAML::BeginMap;

    out << YAML::Key << "cell" << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "capacity" << YAML::Value << s.capacity;
    out << YAML::Key << "initial_soc" << YAML::Value << s.initial_soc;
    out << YAML::Key << "neg_stoich_at_soc0" << YAML::Value << s.neg_stoich_at_soc0;
    out << YAML::Key << "neg_stoich_at_soc1" << YAML::Value << s.neg_stoich_at_soc1;
    out << YAML::Key << "pos_stoich_at_soc0" << YAML::Value << s.pos_stoich_at_soc0;
    out << YAML::Key << "pos_stoich_at_soc1" << YAML::Value << s.pos_stoich_at_soc1;
    out << YAML::EndMap;

    out << YAML::Key << "kinetics" << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "temperature" << YAML::Value << s.temperature;
    out << YAML::EndMap;

    out << YAML::Key << "negative_electrode" << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "thickness" << YAML::Value << s.neg_thickness;
    out << YAML::Key << "particle_radius" << YAML::Value << s.neg_particle_radius;
    out << YAML::Key << "active_fraction" << YAML::Value << s.neg_active_fraction;
    out << YAML::Key << "diffusion_coefficient" << YAML::Value << s.neg_diffusion;
    out << YAML::Key << "max_concentration" << YAML::Value << s.neg_max_conc;
    out << YAML::Key << "exchange_current_density" << YAML::Value << s.neg_m_ref;
    out << YAML::Key << "activation_energy" << YAML::Value << s.neg_activation_energy;
    if (with_conductivity) {
        out << YAML::Key << "conductivity" << YAML::Value << p.neg_conductivity;
    }
    out << YAML::EndMap;

    out << YAML::Key << "positive_electrode" << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "thickness" << YAML::Value << s.pos_thickness;
    out << YAML::Key << "particle_radius" << YAML::Value << s.pos_particle_radius;
    out << YAML::Key << "active_fraction" << YAML::Value << s.pos_active_fraction;
    out << YAML::Key << "diffusion_coefficient" << YAML::Value << s.pos_diffusion;
    out << YAML::Key << "max_concentration" << YAML::Value << s.pos_max_conc;
    out << YAML::Key << "exchange_current_density" << YAML::Value << s.pos_m_ref;
    out << YAML::Key << "activation_energy" << YAML::Value << s.pos_activation_energy;
    if (with_conductivity) {
        out << YAML::Key << "conductivity" << YAML::Value << p.pos_conductivity;
    }
    out << YAML::EndMap;

    out << YAML::Key << "separator" << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "thickness" << YAML::Value << s.sep_thickness;
    out << YAML::EndMap;

    out << YAML::Key << "electrolyte" << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "concentration" << YAML::Value << s.electrolyte_conc;
    if (with_electrolyte) {
        out << YAML::Key << "initial_concentration" << YAML::Value << p.c_e_init;
        out << YAML::Key << "diffusivity" << YAML::Value << p.electrolyte_diffusivity;
        out << YAML::Key << "transference_number" << YAML::Value << p.transference_number;
        out << YAML::Key << "porosity_neg" << YAML::Value << p.porosity_neg;
        out << YAML::Key << "porosity_sep" << YAML::Value << p.porosity_sep;
        out << YAML::Key << "porosity_pos" << YAML::Value << p.porosity_pos;
        out << YAML::Key << "bruggeman_exponent" << YAML::Value << p.bruggeman_exponent;
        out << YAML::Key << "conductivity" << YAML::Value << p.electrolyte_conductivity;
        // 浓度相关输运表（非空时写出，与读取段 transport_table 对称）
        if (!p.electrolyte_conc_lut.empty()
            && (!p.conductivity_lut.empty() || !p.diffusivity_lut.empty())) {
            out << YAML::Key << "transport_table" << YAML::Value << YAML::BeginMap;
            out << YAML::Key << "nodes" << YAML::Value;
            emit_seq(p.electrolyte_conc_lut);
            if (!p.conductivity_lut.empty()) {
                out << YAML::Key << "conductivity" << YAML::Value;
                emit_seq(p.conductivity_lut);
            }
            if (!p.diffusivity_lut.empty()) {
                out << YAML::Key << "diffusivity" << YAML::Value;
                emit_seq(p.diffusivity_lut);
            }
            out << YAML::EndMap;
        }
    }
    out << YAML::EndMap;

    out << YAML::Key << "geometry" << YAML::Value << YAML::BeginMap;
    out << YAML::Key << "electrode_height" << YAML::Value << s.electrode_height;
    out << YAML::Key << "electrode_width" << YAML::Value << s.electrode_width;
    out << YAML::Key << "particle_mesh_cells" << YAML::Value << s.particle_cells;
    if (with_electrolyte) {
        out << YAML::Key << "electrolyte_mesh_cells" << YAML::Value << p.electrolyte_cells;
    }
    out << YAML::EndMap;

    // OCP 曲线（可选；非空且等长时写出）
    if (!s.neg_ocp_stoich.empty() || !s.pos_ocp_stoich.empty()) {
        out << YAML::Key << "ocp" << YAML::Value << YAML::BeginMap;
        const auto emit_ocp = [&out, &emit_seq](
                const char* name, const std::vector<real_t>& stoich,
                const std::vector<real_t>& voltage) {
            if (stoich.empty() || stoich.size() != voltage.size()) {
                return;
            }
            out << YAML::Key << name << YAML::Value << YAML::BeginMap;
            out << YAML::Key << "stoich" << YAML::Value;
            emit_seq(stoich);
            out << YAML::Key << "voltage" << YAML::Value;
            emit_seq(voltage);
            out << YAML::EndMap;
        };
        emit_ocp("negative_electrode", s.neg_ocp_stoich, s.neg_ocp_voltage);
        emit_ocp("positive_electrode", s.pos_ocp_stoich, s.pos_ocp_voltage);
        out << YAML::EndMap;
    }

    // 热模型段（DFN 专属；与读取段 thermal 对称。等温时仍写出全部热参数，
    // 保证 round-trip 无损与用户配置保留）
    if (with_conductivity) {
        const auto& th = p.thermal_params;
        out << YAML::Key << "thermal" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "type" << YAML::Value
            << (p.thermal_type == model::ThermalType::Lumped
                    ? "lumped" : "isothermal");
        out << YAML::Key << "ambient_temperature"
            << YAML::Value << th.t_ambient;
        out << YAML::Key << "initial_temperature"
            << YAML::Value << th.t_initial;
        out << YAML::Key << "heat_transfer_coefficient"
            << YAML::Value << th.h_total;
        out << YAML::Key << "cooling_area" << YAML::Value << th.cooling_area;
        out << YAML::Key << "cell_volume" << YAML::Value << th.cell_volume;
        out << YAML::Key << "rho_cp_negative"
            << YAML::Value << th.rho_c_p_neg;
        out << YAML::Key << "rho_cp_separator"
            << YAML::Value << th.rho_c_p_sep;
        out << YAML::Key << "rho_cp_positive"
            << YAML::Value << th.rho_c_p_pos;
        out << YAML::Key << "rho_cp_cc_negative"
            << YAML::Value << th.rho_c_p_cc_neg;
        out << YAML::Key << "rho_cp_cc_positive"
            << YAML::Value << th.rho_c_p_cc_pos;
        out << YAML::Key << "cc_thickness_negative"
            << YAML::Value << th.cc_thickness_neg;
        out << YAML::Key << "cc_thickness_positive"
            << YAML::Value << th.cc_thickness_pos;
        // 熵变系数：LUT 非空时内联写出（sto/dudt，V/K），否则写常数
        if (!th.entropic_lut_neg.empty()) {
            out << YAML::Key << "entropic_change_negative" << YAML::Value;
            out << YAML::Flow << YAML::BeginMap;
            out << YAML::Key << "sto" << YAML::Value << th.entropic_lut_neg.sto;
            out << YAML::Key << "dudt" << YAML::Value << th.entropic_lut_neg.dudt;
            out << YAML::EndMap;
            out << YAML::Block;
        } else {
            out << YAML::Key << "entropic_change_negative"
                << YAML::Value << th.entropic_change_neg;
        }
        if (!th.entropic_lut_pos.empty()) {
            out << YAML::Key << "entropic_change_positive" << YAML::Value;
            out << YAML::Flow << YAML::BeginMap;
            out << YAML::Key << "sto" << YAML::Value << th.entropic_lut_pos.sto;
            out << YAML::Key << "dudt" << YAML::Value << th.entropic_lut_pos.dudt;
            out << YAML::EndMap;
            out << YAML::Block;
        } else {
            out << YAML::Key << "entropic_change_positive"
                << YAML::Value << th.entropic_change_pos;
        }
        out << YAML::EndMap;

        // SEI 降解段（DFN 专属；始终写出全部键，保证 round-trip 无损）
        const auto& deg = p.degradation;
        out << YAML::Key << "degradation" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "enabled" << YAML::Value << deg.enabled;
        out << YAML::Key << "sei" << YAML::Value << YAML::BeginMap;
        out << YAML::Key << "initial_thickness"
            << YAML::Value << deg.l_sei_0;
        out << YAML::Key << "molar_volume" << YAML::Value << deg.v_bar_sei;
        out << YAML::Key << "lithium_moles_ratio" << YAML::Value << deg.z_sei;
        out << YAML::Key << "u_sei" << YAML::Value << deg.u_sei;
        out << YAML::Key << "resistivity" << YAML::Value << deg.r_sei;
        out << YAML::Key << "rate_constant" << YAML::Value << deg.k_sei;
        out << YAML::Key << "c_ec_0" << YAML::Value << deg.c_ec_0;
        out << YAML::Key << "diffusivity" << YAML::Value << deg.d_ec;
        out << YAML::Key << "activation_energy" << YAML::Value << deg.e_sei;
        out << YAML::Key << "reference_temperature" << YAML::Value << deg.t_ref;
        out << YAML::Key << "alpha" << YAML::Value << deg.alpha_sei;
        out << YAML::EndMap;
        out << YAML::EndMap;
    }

    out << YAML::EndMap;
}

/// 写出 emitter 内容；打开/写入失败抛异常（GUI 侧统一捕获展示）
void write_yaml_file(const std::filesystem::path& yaml_path,
                     const YAML::Emitter& out) {
    std::ofstream ofs(yaml_path, std::ios::binary);
    if (!ofs.is_open()) {
        throw std::runtime_error("无法打开文件写入: " + yaml_path.string());
    }
    ofs << out.c_str() << "\n";
    if (!ofs) {
        throw std::runtime_error("文件写入失败: " + yaml_path.string());
    }
}

} // namespace

void save_spm_parameters(const std::filesystem::path& yaml_path,
                         const model::SPMModel::Parameters& params) {
    model::DFNModel::Parameters full; // 仅取固相部分，其余字段不写出
    full.solid = params;
    YAML::Emitter out;
    out << YAML::Comment("FluxCell SPM 模型参数");
    emit_parameters(out, full, /*with_electrolyte=*/false,
                    /*with_conductivity=*/false);
    write_yaml_file(yaml_path, out);
}

void save_spme_parameters(const std::filesystem::path& yaml_path,
                          const model::SPMeModel::Parameters& params) {
    model::DFNModel::Parameters full; // 仅取固相 + 电解液部分
    static_cast<model::SPMeModel::Parameters&>(full) = params;
    YAML::Emitter out;
    out << YAML::Comment("FluxCell SPMe 模型参数");
    emit_parameters(out, full, /*with_electrolyte=*/true,
                    /*with_conductivity=*/false);
    write_yaml_file(yaml_path, out);
}

void save_dfn_parameters(const std::filesystem::path& yaml_path,
                         const model::DFNModel::Parameters& params) {
    YAML::Emitter out;
    out << YAML::Comment("FluxCell DFN 模型参数");
    emit_parameters(out, params, /*with_electrolyte=*/true,
                    /*with_conductivity=*/true);
    write_yaml_file(yaml_path, out);
}

// ══════════════════════════════════════════════════════════════
//  工况文件（多段协议 + 运行条件）
// ══════════════════════════════════════════════════════════════

const char* protocol_mode_name(CCStep::Mode m) {
    switch (m) {
        case CCStep::Mode::Discharge: return "discharge";
        case CCStep::Mode::Charge:    return "charge";
        default:                      return "rest";
    }
}

const char* protocol_control_name(CCStep::Control c) {
    switch (c) {
        case CCStep::Control::Current: return "current";
        case CCStep::Control::Voltage: return "voltage";
        default:                       return "power";
    }
}

const char* protocol_model_name(CellModelType m) {
    switch (m) {
        case CellModelType::SPM:  return "SPM";
        case CellModelType::SPMe: return "SPMe";
        default:                  return "DFN";
    }
}

const char* protocol_solver_name(SolverType s) {
    return s == SolverType::SundialsIDA ? "SundialsIDA" : "EigenCN";
}

const char* protocol_thermal_name(model::ThermalType t) {
    return t == model::ThermalType::Lumped ? "lumped" : "isothermal";
}

CCStep::Mode parse_step_mode(const std::string& s, std::size_t index) {
    if (s == "discharge") return CCStep::Mode::Discharge;
    if (s == "charge")    return CCStep::Mode::Charge;
    if (s == "rest")      return CCStep::Mode::Rest;
    throw std::runtime_error("工况段 " + std::to_string(index + 1)
                             + " 的 mode 非法: " + s
                             + "（应为 discharge/charge/rest）");
}

CCStep::Control parse_step_control(const std::string& s, std::size_t index) {
    if (s == "current") return CCStep::Control::Current;
    if (s == "voltage") return CCStep::Control::Voltage;
    if (s == "power")   return CCStep::Control::Power;
    throw std::runtime_error("工况段 " + std::to_string(index + 1)
                             + " 的 control 非法: " + s
                             + "（应为 current/voltage/power）");
}

CellModelType parse_cell_model(const std::string& s) {
    if (s == "SPM")  return CellModelType::SPM;
    if (s == "SPMe") return CellModelType::SPMe;
    if (s == "DFN")  return CellModelType::DFN;
    throw std::runtime_error("工况文件 model 非法: " + s
                             + "（应为 SPM/SPMe/DFN）");
}

SolverType parse_solver(const std::string& s) {
    if (s == "EigenCN")     return SolverType::EigenCN;
    if (s == "SundialsIDA") return SolverType::SundialsIDA;
    throw std::runtime_error("工况文件 solver 非法: " + s
                             + "（应为 EigenCN/SundialsIDA）");
}

model::ThermalType parse_thermal_type(const std::string& s) {
    if (s == "isothermal") return model::ThermalType::Isothermal;
    if (s == "lumped")     return model::ThermalType::Lumped;
    throw std::runtime_error("工况文件 thermal_type 非法: " + s
                             + "（应为 isothermal/lumped）");
}

void save_protocol(const std::filesystem::path& yaml_path,
                   const ProtocolFile& pf) {
    YAML::Emitter out;
    out << YAML::Comment("FluxCell 工况文件：多段协议 + 运行条件（不含电芯参数）");
    out << YAML::BeginMap;
    out << YAML::Key << "model" << YAML::Value
        << protocol_model_name(pf.protocol.model);
    out << YAML::Key << "solver" << YAML::Value
        << protocol_solver_name(pf.protocol.solver);
    out << YAML::Key << "solver_rel_tol" << YAML::Value
        << pf.protocol.solver_rel_tol;
    out << YAML::Key << "thermal_type" << YAML::Value
        << protocol_thermal_name(pf.thermal_type);
    out << YAML::Key << "initial_soc" << YAML::Value << pf.initial_soc;
    out << YAML::Key << "ambient_temp_c" << YAML::Value << pf.ambient_temp_c;
    out << YAML::Key << "cycles" << YAML::Value << pf.protocol.cycles;
    out << YAML::Key << "time_step" << YAML::Value << pf.protocol.time_step;
    out << YAML::Key << "real_time_factor" << YAML::Value
        << pf.protocol.real_time_factor;
    out << YAML::Key << "steps" << YAML::Value << YAML::BeginSeq;
    for (const auto& s : pf.protocol.steps) {
        out << YAML::BeginMap;
        out << YAML::Key << "mode" << YAML::Value << protocol_mode_name(s.mode);
        out << YAML::Key << "control" << YAML::Value
            << protocol_control_name(s.control);
        out << YAML::Key << "current" << YAML::Value << s.current;
        out << YAML::Key << "voltage" << YAML::Value << s.voltage;
        out << YAML::Key << "power" << YAML::Value << s.power;
        out << YAML::Key << "cutoff_time" << YAML::Value << s.cutoff_time;
        out << YAML::Key << "cutoff_voltage" << YAML::Value << s.cutoff_voltage;
        out << YAML::Key << "cutoff_current" << YAML::Value << s.cutoff_current;
        out << YAML::EndMap;
    }
    out << YAML::EndSeq;
    out << YAML::EndMap;
    write_yaml_file(yaml_path, out);
}

ProtocolFile load_protocol(const std::filesystem::path& yaml_path) {
    YAML::Node root = YAML::LoadFile(yaml_path.string());
    ProtocolFile pf;

    if (auto n = root["model"])
        pf.protocol.model = parse_cell_model(n.as<std::string>());
    if (auto n = root["solver"])
        pf.protocol.solver = parse_solver(n.as<std::string>());
    if (auto n = root["solver_rel_tol"])
        pf.protocol.solver_rel_tol = n.as<real_t>();
    if (auto n = root["thermal_type"])
        pf.thermal_type = parse_thermal_type(n.as<std::string>());
    if (auto n = root["initial_soc"])
        pf.initial_soc = n.as<real_t>();
    if (auto n = root["ambient_temp_c"])
        pf.ambient_temp_c = n.as<real_t>();
    if (auto n = root["cycles"])
        pf.protocol.cycles = n.as<int>();
    if (auto n = root["time_step"])
        pf.protocol.time_step = n.as<real_t>();
    if (auto n = root["real_time_factor"])
        pf.protocol.real_time_factor = n.as<real_t>();

    if (auto steps = root["steps"]) {
        pf.protocol.steps.clear();
        for (const auto& n : steps) {
            const std::size_t index = pf.protocol.steps.size();
            CCStep s;
            if (auto v = n["mode"])
                s.mode = parse_step_mode(v.as<std::string>(), index);
            if (auto v = n["control"])
                s.control = parse_step_control(v.as<std::string>(), index);
            if (auto v = n["current"])   s.current = v.as<real_t>();
            if (auto v = n["voltage"])   s.voltage = v.as<real_t>();
            if (auto v = n["power"])     s.power = v.as<real_t>();
            if (auto v = n["cutoff_time"])    s.cutoff_time = v.as<real_t>();
            if (auto v = n["cutoff_voltage"]) s.cutoff_voltage = v.as<real_t>();
            if (auto v = n["cutoff_current"]) s.cutoff_current = v.as<real_t>();
            pf.protocol.steps.push_back(s);
        }
    }
    return pf;
}

// ══════════════════════════════════════════════════════════════
//  HDF5 写入（占位实现）
// ══════════════════════════════════════════════════════════════

struct HDF5Writer::Impl {
    std::filesystem::path filepath;
    bool initialized = false;
    int step_count = 0;
};

HDF5Writer::HDF5Writer(const std::filesystem::path& filepath)
    : impl_(std::make_unique<Impl>()) {
    impl_->filepath = filepath;
}

HDF5Writer::~HDF5Writer() {
    close();
}

void HDF5Writer::initialize(int n_c_e, int n_phi, int n_particle) {
    // TODO: 创建 HDF5 文件和数据集
    // - /time (unlimited dimension)
    // - /voltage
    // - /c_e (n_c_e × unlimited)
    // - /phi_e (n_phi × unlimited)
    // - /c_s_neg (n_particle × n_neg × unlimited)
    // - /c_s_pos (n_particle × n_pos × unlimited)
    (void)n_c_e;
    (void)n_phi;
    (void)n_particle;
    impl_->initialized = true;
    impl_->step_count = 0;
}

void HDF5Writer::append(const SolutionSnapshot& snapshot) {
    if (!impl_->initialized) return;
    // TODO: 追加数据到 HDF5 数据集
    (void)snapshot;
    impl_->step_count++;
}

void HDF5Writer::close() {
    impl_->initialized = false;
}

// ══════════════════════════════════════════════════════════════
//  CSV 导出
// ══════════════════════════════════════════════════════════════

void export_csv(const std::filesystem::path& filepath,
                const std::vector<SolutionSnapshot>& snapshots,
                const CellGeometry& geometry) {
    (void)geometry;
    std::ofstream ofs(filepath);
    if (!ofs.is_open()) {
        throw std::runtime_error("Cannot open file: " + filepath.string());
    }

    // 若快照包含颗粒表面浓度（SPM），则追加两列
    const bool has_surf = !snapshots.empty()
                       && snapshots.front().c_s_surf_neg.size() > 0;

    ofs << std::scientific << std::setprecision(8);
    if (has_surf) {
        ofs << "time,voltage,current,temperature,c_s_surf_neg,c_s_surf_pos\n";
    } else {
        ofs << "time,voltage,current,temperature\n";
    }

    for (const auto& snap : snapshots) {
        ofs << snap.time << "," << snap.voltage << ","
            << snap.current << "," << snap.temperature;
        if (has_surf) {
            ofs << "," << snap.c_s_surf_neg(0) << "," << snap.c_s_surf_pos(0);
        }
        ofs << "\n";
    }
}

} // namespace fluxcell::io
