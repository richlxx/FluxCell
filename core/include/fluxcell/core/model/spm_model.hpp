/// @file spm_model.hpp
/// @brief 单颗粒模型 (Single Particle Model, SPM)
///
/// 参考 PyBaMM `pybamm.lithium_ion.SPM` + Chen2020 参数集：
/// - 每个电极一个代表性球形颗粒，径向 Fick 扩散（球坐标 FVM）
/// - 电解液浓度恒定（SPM 标准简化），无液相/固相欧姆损失
/// - 对称 Butler-Volmer 动力学，arcsinh 解析反解过电位
/// - 端电压 V = U_p(y_surf) - U_n(x_surf) + η_p - η_n
#pragma once

#include <fluxcell/core/model/cell_model.hpp>
#include <fluxcell/core/types.hpp>

#include <Eigen/SparseLU>
#include <memory>
#include <vector>

namespace fluxcell::model {

/// 单颗粒模型
class SPMModel : public CellModel, public solver::DaeProblem {
public:
    /// SPM 参数（默认值 = PyBaMM Chen2020 / LG M50）
    struct Parameters {
        // ── 几何 ──
        real_t neg_thickness = 85.2e-6;        ///< 负极厚度 L_n [m]
        real_t sep_thickness = 12.0e-6;        ///< 隔膜厚度 [m]（SPM 不参与计算）
        real_t pos_thickness = 75.6e-6;        ///< 正极厚度 L_p [m]
        real_t electrode_height = 0.065;       ///< 电极高度 [m]
        real_t electrode_width = 1.58;         ///< 电极宽度 [m]
        real_t neg_particle_radius = 5.86e-6;  ///< 负极颗粒半径 R_n [m]
        real_t pos_particle_radius = 5.22e-6;  ///< 正极颗粒半径 R_p [m]

        // ── 负极（石墨）──
        real_t neg_active_fraction = 0.75;     ///< 活性物质体积分数 ε_s,n
        real_t neg_diffusion = 3.3e-14;        ///< D_s,n [m²/s]
        real_t neg_max_conc = 33133.0;         ///< c_s_max,n [mol/m³]
        real_t neg_m_ref = 6.48e-7;            ///< 交换电流密度参考值 [(A/m²)(m³/mol)^1.5]
        real_t neg_activation_energy = 35000.0;///< E_r,n [J/mol]

        // ── 正极（NMC）──
        real_t pos_active_fraction = 0.665;    ///< ε_s,p
        real_t pos_diffusion = 4.0e-15;        ///< D_s,p [m²/s]
        real_t pos_max_conc = 63104.0;         ///< c_s_max,p [mol/m³]
        real_t pos_m_ref = 3.42e-6;            ///< m_ref,p
        real_t pos_activation_energy = 17800.0;///< E_r,p [J/mol]

        // ── 电解液 / 环境 ──
        real_t electrolyte_conc = 1000.0;      ///< c_e（恒定）[mol/m³]
        real_t temperature = 298.15;           ///< 等温环境温度 T [K]
        real_t capacity = 5.0;                 ///< 名义容量 [Ah]

        // ── 初始条件：初始 SOC 线性映射化学计量端点 ──
        // 端点可由 YAML cell: 段覆盖（neg/pos_stoich_at_soc0/soc1）
        // Chen2020 (LG M50)：PyBaMM initial_soc 反推值（由 0%/100% SOC 平衡
        // 电压 2.5V/4.2V 与电极锂守恒 C_n(x₁−x₀)=C_p(y₀−y₁) 联立求解）。
        // 论文表 4 端点 (0.0279, 0.9014, 0.9084, 0.2661) 与电极几何容量比
        // 不自洽 —— 直抄会使 0% SOC 正极嵌入率偏高约 0.054，同 SOC 下 OCV
        // 偏低 60~100mV，CC-CV 充电多充约 9% SOC（CV 末端 SOC 虚高至 110%）
        real_t initial_soc = 1.0;              ///< 初始荷电状态 0~1
        real_t neg_stoich_at_soc0 = 0.026346;  ///< 负极 0% SOC 嵌入率 x_n
        real_t neg_stoich_at_soc1 = 0.910618;  ///< 负极 100% SOC 嵌入率 x_n
        real_t pos_stoich_at_soc0 = 0.853975;  ///< 正极 0% SOC 嵌入率 y_p
        real_t pos_stoich_at_soc1 = 0.263845;  ///< 正极 100% SOC 嵌入率 y_p

        // ── 数值 ──
        int particle_cells = 20;               ///< 颗粒径向控制体数

        // ── OCP 曲线（可选；YAML "ocp:" 段加载，空 = 内置 Chen2020 公式）──
        std::vector<real_t> neg_ocp_stoich;    ///< 负极 U_n(x) 嵌入率采样点（升序）
        std::vector<real_t> neg_ocp_voltage;   ///< 负极 U_n 电压采样点 [V]
        std::vector<real_t> pos_ocp_stoich;    ///< 正极 U_p(y) 嵌入率采样点（升序）
        std::vector<real_t> pos_ocp_voltage;   ///< 正极 U_p 电压采样点 [V]

        /// 负极开路电位 U_n(x)：优先查表（线性插值、端点钳位），空表用内置公式
        real_t neg_ocp_at(real_t x) const;
        /// 正极开路电位 U_p(y)：优先查表（线性插值、端点钳位），空表用内置公式
        real_t pos_ocp_at(real_t y) const;
    };

    /// @param params 模型参数
    /// @param initial_current 初始工作电流 [A]（正 = 放电；用于计算 t=0 的过电位，
    ///                        与 PyBaMM 阶跃工况对齐；0 = 纯 OCV）
    explicit SPMModel(const Parameters& params, real_t initial_current = 0.0);
    ~SPMModel();

    /// 以恒定电流推进一个时间步（Crank-Nicolson 扩散 + BV 动力学）
    /// @param current  外加电流 [A]（正 = 放电）
    /// @param dt       时间步长 [s]
    void step(real_t current, real_t dt) override;

    /// 当前端电压 [V]（基于最近一次 step 的表面浓度/过电位）
    real_t voltage() const override { return voltage_; }
    real_t time() const override { return time_; }

    /// 负极颗粒表面浓度 [mol/m³]（FVM 边界值：由最外两个控制体二阶外推，
    /// 与 PyBaMM 的 edge 值约定对齐；控制体平均用于锂守恒）
    real_t neg_surface_conc() const;
    /// 正极颗粒表面浓度 [mol/m³]（同上）
    real_t pos_surface_conc() const;

    /// 负极嵌入率 x = c_s_surf / c_s_max
    real_t neg_stoichiometry() const { return neg_surface_conc() / params_.neg_max_conc; }
    /// 正极嵌入率 y
    real_t pos_stoichiometry() const { return pos_surface_conc() / params_.pos_max_conc; }

    /// 最近一步的过电位 [V]
    real_t neg_overpotential() const { return eta_neg_; }
    real_t pos_overpotential() const { return eta_pos_; }

    /// 最近一步的电极开路电位 U_n(x_surf) / U_p(y_surf) [V]
    real_t neg_ocp() const { return u_neg_; }
    real_t pos_ocp() const { return u_pos_; }

    /// 电极比表面积 a = 3·ε_s/R [m²/m³]
    real_t neg_specific_area() const { return a_neg_; }
    real_t pos_specific_area() const { return a_pos_; }

    /// 颗粒内锂总量（mol，两电极合计×各自颗粒数归一）— 用于守恒性测试
    real_t total_lithium_mol() const;

    /// 径向浓度分布（调试/可视化）
    const Vector& neg_concentration() const { return c_s_neg_; }
    const Vector& pos_concentration() const { return c_s_pos_; }

    /// 生成当前状态的仿真快照
    SolutionSnapshot snapshot(real_t current) const override;

    const Parameters& parameters() const { return params_; }

    real_t nominal_capacity() const override { return params_.capacity; }
    real_t initial_soc() const override { return params_.initial_soc; }

    // ── solver::DaeProblem 接口（SPM 为纯 ODE，无代数行）──
    int dae_size() const override;
    void dae_mask(Vector& id) const override;
    void dae_atol(Vector& atol) const override;
    void dae_pack_state(Vector& y, Vector& yp) const override;
    void dae_residual(real_t t, const Eigen::Ref<const Vector>& y,
                      const Eigen::Ref<const Vector>& yp,
                      Eigen::Ref<Vector> res) override;
    void dae_load_state(real_t t,
                                const Eigen::Ref<const Vector>& y) override;
    void dae_set_current(real_t current) override { dae_current_ = current; }
    solver::DaeProblem* dae_problem() override { return this; }

private:
    /// 由表面浓度计算过电位与端电压
    void update_voltage(real_t current);

    /// 构建/重建 Crank-Nicolson 算子（dt 变化时调用）
    void ensure_cn(real_t dt);

    Parameters params_;
    int n_ = 0;                ///< 颗粒径向控制体数
    real_t a_neg_ = 0.0;       ///< 负极比表面积 [m²/m³]
    real_t a_pos_ = 0.0;
    real_t area_ = 0.0;        ///< 电极面积 A [m²]

    Vector c_s_neg_;           ///< 负极颗粒径向浓度 [mol/m³]
    Vector c_s_pos_;           ///< 正极颗粒径向浓度 [mol/m³]

    SparseMatrix D_neg_;       ///< 负极颗粒球坐标扩散算子
    SparseMatrix D_pos_;       ///< 正极颗粒球坐标扩散算子
    std::vector<real_t> vol_neg_; ///< 负极颗粒控制体体积 [m³]
    std::vector<real_t> vol_pos_; ///< 正极颗粒控制体体积 [m³]

    // Crank-Nicolson 预分解矩阵（dt 变化时重新分解）
    struct CNOperators;
    std::unique_ptr<CNOperators> cn_;
    real_t cn_dt_ = -1.0;

    real_t voltage_ = 0.0;
    real_t eta_neg_ = 0.0;
    real_t eta_pos_ = 0.0;
    real_t u_neg_ = 0.0;
    real_t u_pos_ = 0.0;
    real_t time_ = 0.0;
    real_t dae_current_ = 0.0; ///< DAE 残差的端口电流 [A]（段切换时更新）
};

} // namespace fluxcell::model
