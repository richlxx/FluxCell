/// @file types.hpp
/// @brief FluxCell 核心类型定义
#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cstdint>
#include <string>
#include <vector>

namespace fluxcell {

/// 浮点精度（默认 double）
using real_t = double;

/// 动态大小向量
using Vector = Eigen::Matrix<real_t, Eigen::Dynamic, 1>;

/// 稀疏矩阵（列主序）
using SparseMatrix = Eigen::SparseMatrix<real_t, Eigen::ColMajor>;

/// 三元组列表（用于稀疏矩阵组装）
using TripletList = std::vector<Eigen::Triplet<real_t>>;

/// 电池域标识
enum class Domain : uint8_t {
    NegativeElectrode = 0,
    Separator = 1,
    PositiveElectrode = 2,
    Count
};

/// 仿真工况类型
enum class ProtocolType : uint8_t {
    ConstantCurrent,   ///< 恒流充放电 (CC)
    ConstantVoltage,   ///< 恒压 (CV)
    CCCV,              ///< 恒流-恒压
    Pulse,             ///< 脉冲工况
    Custom             ///< 自定义电流曲线
};

/// 仿真状态
enum class SimulationStatus : uint8_t {
    Idle,
    Running,
    Paused,
    Completed,
    Failed
};

/// 电池几何参数
struct CellGeometry {
    real_t negative_thickness = 80e-6;   ///< 负极厚度 [m]
    real_t separator_thickness = 25e-6;  ///< 隔膜厚度 [m]
    real_t positive_thickness = 70e-6;   ///< 正极厚度 [m]
    real_t electrode_height = 0.208;     ///< 电极高度 [m]
    real_t electrode_width = 0.137;      ///< 电极宽度 [m]
    real_t negative_particle_radius = 10e-6; ///< 负极颗粒半径 [m]
    real_t positive_particle_radius = 10e-6; ///< 正极颗粒半径 [m]
    int negative_mesh_cells = 20;
    int separator_mesh_cells = 10;
    int positive_mesh_cells = 20;
    int particle_mesh_cells = 15;        ///< 颗粒径向网格数
};

/// 仿真结果快照（单时间步）
struct SolutionSnapshot {
    real_t time = 0.0;
    real_t voltage = 0.0;
    real_t current = 0.0;
    real_t temperature = 298.15;
    real_t soc = 0.0;     ///< 荷电状态 0~1（库仑计数，由引擎填充）
    real_t internal_step = 0.0; ///< 求解器最近内部步长 [s]（IDA 诊断；CN = dt）
    Vector c_e;           ///< 电解液浓度分布（SPM 恒定，SPMe/DFN 有效）
    Vector phi_s;         ///< 固相电位分布（SPMe/DFN 有效）
    Vector phi_e;         ///< 电解液电位分布（SPMe/DFN 有效）
    Vector c_s_surf_neg;  ///< 负极表面固相浓度
    Vector c_s_surf_pos;  ///< 正极表面固相浓度
    Vector c_s_neg;       ///< 负极颗粒径向浓度分布 c_s(r)（控制体平均）
    Vector c_s_pos;       ///< 正极颗粒径向浓度分布 c_s(r)
    real_t x_surf = 0.0;  ///< 负极表面嵌入率 x = c_surf / c_max
    real_t y_surf = 0.0;  ///< 正极表面嵌入率 y
    real_t eta_neg = 0.0; ///< 负极活化过电位 η_n [V]
    real_t eta_pos = 0.0; ///< 正极活化过电位 η_p [V]
    real_t eta_e = 0.0;   ///< 液相过电位（欧姆 + 浓差，SPM 恒为 0）[V]
    real_t u_neg = 0.0;   ///< 负极开路电位 U_n(x_surf) [V]
    real_t u_pos = 0.0;   ///< 正极开路电位 U_p(y_surf) [V]

    // ── 产热诊断（DFN；PyBaMM base_thermal 三项口径）[W] ──
    real_t q_ohm = 0.0;   ///< 欧姆热（固相 i²/σ^eff + 液相 −i_e·∂φ_e/∂x）
    real_t q_rxn = 0.0;   ///< 不可逆界面热 a·i_surf·η
    real_t q_rev = 0.0;   ///< 可逆熵热 a·i_surf·T·dU/dT
    real_t q_total = 0.0; ///< 总产热（q_ohm + q_rxn + q_rev）

    // ── DFN 分布输出（SPM/SPMe 留空）──
    int sep_start = -1;   ///< 隔膜区起始格索引（SPMe/DFN 填充；<0 无区间信息）
    int pos_start = -1;   ///< 正极区起始格索引（SPMe/DFN 填充）
    Vector c_surf_neg_x;  ///< 负极各 x 格颗粒表面浓度 [mol/m³]
    Vector c_surf_pos_x;  ///< 正极各 x 格颗粒表面浓度 [mol/m³]
    Vector i_surf_neg;    ///< 负极孔壁电流密度分布 [A/m²]
    Vector i_surf_pos;    ///< 正极孔壁电流密度分布 [A/m²]
    Vector eta_neg_x;     ///< 负极活化过电位分布 [V]
    Vector eta_pos_x;     ///< 正极活化过电位分布 [V]
    Vector u_neg_x;       ///< 负极 OCP 分布 [V]
    Vector u_pos_x;       ///< 正极 OCP 分布 [V]
    Vector i_e_face;      ///< 液相面电流密度（格界面，m+1）[A/m²]

    // ── 机理寿命诊断（DFN 降解启用时非零；engine 填循环圈数）──
    int cycle_number = 0;          ///< 当前循环圈数（1-based）
    bool cycle_end = false;        ///< 本快照为某圈最后一段的结束时刻（寿命曲线采样点）
    real_t li_loss_ah = 0.0;       ///< SEI 累计锂损失 [Ah]（PyBaMM Loss of capacity）
    real_t degraded_capacity_ah = 0.0; ///< 含降解估计容量 [Ah]
    real_t sei_thickness = 0.0;    ///< 负极 SEI 平均厚度 [m]
};

} // namespace fluxcell
