/// @file io.hpp
/// @brief 参数读取/保存（YAML）与结果写入（HDF5 / CSV）
#pragma once

#include <fluxcell/core/model/dfn_model.hpp>
#include <fluxcell/core/model/spm_model.hpp>
#include <fluxcell/core/model/spme_model.hpp>
#include <fluxcell/core/simulation/simulation.hpp>
#include <fluxcell/core/types.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace fluxcell::io {

// ══════════════════════════════════════════════════════════════
//  YAML 参数读取
// ══════════════════════════════════════════════════════════════

/// 从 CSV 加载熵变系数 dU/dT 查找表（两列：嵌锂量 sto、dU/dT [mV/K]；
/// 表头/空行/解析失败行自动跳过，值换算为 V/K，有效点需 ≥ 2）
model::EntropyLut load_entropy_lut_csv(const std::filesystem::path& csv_path);

/// 从 YAML 文件加载 DFN 模型参数（SPMe 键 + negative/positive_electrode
/// 的 conductivity 固相电导键；缺失键取默认值）
model::DFNModel::Parameters load_dfn_parameters(const std::filesystem::path& yaml_path);

/// 从 YAML 文件加载电池几何参数
CellGeometry load_cell_geometry(const std::filesystem::path& yaml_path);

/// 从 YAML 文件加载 SPM 模型参数（Chen2020 格式）
model::SPMModel::Parameters load_spm_parameters(const std::filesystem::path& yaml_path);

/// 从 YAML 文件加载 SPMe 模型参数（SPM 键 + electrolyte 扩展键；
/// 缺失的电解液键取默认值，固相部分与 load_spm_parameters 完全一致）
model::SPMeModel::Parameters load_spme_parameters(const std::filesystem::path& yaml_path);

/// 将 SPM 模型参数保存为 YAML 文件（键结构与 load_spm_parameters 对称，
/// 含非空且等长的 OCP 曲线表；可直接被 load_spm_parameters 读回）
void save_spm_parameters(const std::filesystem::path& yaml_path,
                         const model::SPMModel::Parameters& params);

/// 将 SPMe 模型参数保存为 YAML 文件（SPM 部分 + electrolyte 扩展段与
/// 可选 transport_table；键结构与 load_spme_parameters 对称）
void save_spme_parameters(const std::filesystem::path& yaml_path,
                          const model::SPMeModel::Parameters& params);

/// 将 DFN 模型参数保存为 YAML 文件（在 SPMe 基础上追加
/// negative/positive_electrode 的 conductivity 固相电导键）
void save_dfn_parameters(const std::filesystem::path& yaml_path,
                         const model::DFNModel::Parameters& params);

// ══════════════════════════════════════════════════════════════
//  工况文件（多段协议 + 运行条件）
// ══════════════════════════════════════════════════════════════

/// 工况文件内容：多段协议 + 面板运行条件（不含电芯化学参数——
/// 参数仍由独立参数文件提供，两层文件组合成完整仿真配置）
struct ProtocolFile {
    MultiStepProtocol protocol;    ///< 多段协议（spme_params/dfn_params 不参与存取）
    real_t initial_soc = 0.0;     ///< 初始 SOC（0~1，与参数文件 cell.initial_soc 同口径）
    real_t ambient_temp_c = 25.0; ///< 环境温度 [℃]（等温环境/集总热 T_amb）
    model::ThermalType thermal_type = model::ThermalType::Isothermal; ///< DFN 热模型
};

/// 将工况保存为 YAML 文件（枚举写成字符串：discharge/charge/rest、
/// current/voltage/power、SPM/SPMe/DFN、EigenCN/SundialsIDA、
/// isothermal/lumped；可直接被 load_protocol 读回）
void save_protocol(const std::filesystem::path& yaml_path,
                   const ProtocolFile& pf);

/// 从 YAML 文件加载工况（缺失键取默认值；枚举字符串非法时抛出）
ProtocolFile load_protocol(const std::filesystem::path& yaml_path);

// ══════════════════════════════════════════════════════════════
//  HDF5 结果写入
// ══════════════════════════════════════════════════════════════

/// HDF5 结果写入器
///
/// 支持增量写入（每个时间步追加），适合长时间仿真。
class HDF5Writer {
public:
    explicit HDF5Writer(const std::filesystem::path& filepath);
    ~HDF5Writer();

    /// 初始化数据集（在第一次写入前调用）
    void initialize(int n_c_e, int n_phi, int n_particle);

    /// 追加一个时间步的结果
    void append(const SolutionSnapshot& snapshot);

    /// 关闭文件
    void close();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// ══════════════════════════════════════════════════════════════
//  CSV 导出
// ══════════════════════════════════════════════════════════════

/// 将一系列快照导出为 CSV 文件
void export_csv(const std::filesystem::path& filepath,
                const std::vector<SolutionSnapshot>& snapshots,
                const CellGeometry& geometry);

} // namespace fluxcell::io
