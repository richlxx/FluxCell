/// @file spme_model.hpp
/// @brief 扩展单颗粒模型 (Single Particle Model with Electrolyte, SPMe)
///
/// 在 SPM 基础上引入电解液相 1D 沿厚度方向模型（参考 PyBaMM `pybamm.lithium_ion.SPMe`）：
/// - 液相质量守恒：ε_e ∂c_e/∂t = ∂/∂x(D_e^eff ∂c_e/∂x) + (1−t₊)/F · ∂i_e/∂x
///   （源项经 ∂i_e/∂x = a·F·j 展开；D_e^eff = D_e(c_e)·ε_e^bruggeman 浓度相关，
///   单点常量回退）
/// - 液相过电位（逐格解析积分，κ^eff = κ(c_e)·ε^bruggeman）：
///   η_e = φ_e(L) − φ_e(0) = −∫ i_e/κ^eff dx + (2RT/F)(1−t₊)·ln(c_e(L)/c_e(0))
/// - 端电压：V = U_p − U_n + η_p − η_n + η_e（放电时 η_e < 0）
/// - BV 动力学的交换电流密度使用各电极区域的 x 平均电解液浓度
/// - 固相颗粒部分与 SPM 相同：球坐标 FVM + Crank-Nicolson + 对称 BV arcsinh 反解
#pragma once

#include <fluxcell/core/model/cell_model.hpp>
#include <fluxcell/core/model/spm_model.hpp>
#include <fluxcell/core/types.hpp>

#include <Eigen/SparseLU>
#include <memory>

namespace fluxcell::model {

/// 扩展单颗粒模型
class SPMeModel : public CellModel, public solver::DaeProblem {
public:
    /// SPMe 参数：固相部分复用 SPM 参数 + 电解液相扩展
    struct Parameters {
        SPMModel::Parameters solid;            ///< 固相/几何/动力学/初始条件

        // ── 电解液相 ──
        real_t electrolyte_diffusivity = 3.0e-10; ///< D_e [m²/s]（Chen2020）
        real_t transference_number = 0.3634;      ///< 阳离子迁移数 t₊
        real_t porosity_neg = 0.25;               ///< 负极孔隙率 ε_e（= 1 − ε_s,n）
        real_t porosity_sep = 1.0;                ///< 隔膜孔隙率
        real_t porosity_pos = 0.335;              ///< 正极孔隙率（= 1 − ε_s,p）
        real_t bruggeman_exponent = 1.5;          ///< Bruggeman 指数
        real_t electrolyte_conductivity = 1.0;    ///< κ [S/m]（无表时回退常量）
        real_t c_e_init = 1000.0;                 ///< 电解液初始浓度 [mol/m³]
        int electrolyte_cells = 100;              ///< 厚度方向控制体数

        // ── 浓度相关输运表（可选；YAML electrolyte.transport_table 加载）──
        // 非空时 κ/D 按浓度节点线性插值（升序、端点钳位），替代单点常量近似
        //（如 Chen2020 的 Nyman2008 曲线）。两值序列共享同一浓度节点。
        std::vector<real_t> electrolyte_conc_lut; ///< 浓度节点 [mol/m³]（升序）
        std::vector<real_t> conductivity_lut;     ///< κ 采样 [S/m]（空 = 常量）
        std::vector<real_t> diffusivity_lut;      ///< D 采样 [m²/s]（空 = 常量）

        /// κ(c_e)：有表线性插值，无表回退常量 electrolyte_conductivity
        real_t conductivity_at(real_t c_e) const;
        /// D(c_e)：有表线性插值，无表回退常量 electrolyte_diffusivity
        real_t diffusivity_at(real_t c_e) const;
        /// 校验表自洽（值序列与节点等长 ≥2、节点严格升序），不合法抛异常
        void validate_transport_table() const;
    };

    /// @param params 模型参数
    /// @param initial_current 初始工作电流 [A]（正 = 放电；用于 t = 0 过电位）
    explicit SPMeModel(const Parameters& params, real_t initial_current = 0.0);
    ~SPMeModel() override;

    void step(real_t current, real_t dt) override;
    real_t voltage() const override { return voltage_; }
    real_t time() const override { return time_; }
    SolutionSnapshot snapshot(real_t current) const override;

    /// 负极/正极颗粒表面浓度 [mol/m³]（FVM 边界二阶外推，同 SPM）
    real_t neg_surface_conc() const;
    real_t pos_surface_conc() const;

    /// 负极/正极嵌入率
    real_t neg_stoichiometry() const { return neg_surface_conc() / params_.solid.neg_max_conc; }
    real_t pos_stoichiometry() const { return pos_surface_conc() / params_.solid.pos_max_conc; }

    /// 最近一步的过电位 [V]
    real_t neg_overpotential() const { return eta_neg_; }
    real_t pos_overpotential() const { return eta_pos_; }
    real_t electrolyte_overpotential() const { return eta_e_; } ///< 液相欧姆+浓差 [V]

    /// 最近一步的电极开路电位 [V]
    real_t neg_ocp() const { return u_neg_; }
    real_t pos_ocp() const { return u_pos_; }

    /// 电解液浓度分布 / 液相电位分布（相对 φ_e(0) 基准）[mol/m³] / [V]
    const Vector& electrolyte_conc() const { return c_e_; }
    const Vector& electrolyte_potential() const { return phi_e_; }

    /// 颗粒内锂总量 + 液相锂总量（mol）— 用于守恒性测试
    real_t total_lithium_mol() const;

    const Parameters& parameters() const { return params_; }

    real_t nominal_capacity() const override { return params_.solid.capacity; }
    real_t initial_soc() const override { return params_.solid.initial_soc; }

    // ── solver::DaeProblem 接口（SPMe 为纯 ODE，无代数行）──
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
    /// 由表面浓度与电解液分布计算过电位与端电压
    void update_voltage(real_t current);

    /// 构建/重建 Crank-Nicolson 算子（dt 变化时调用）
    void ensure_cn(real_t dt);

    /// 液相守恒通量散度 L(c)·c：D^eff = D(c_e)·ε^bruggeman 逐格、
    /// 界面算术平均，两端零通量；浓度相关 D 的 CN/DAE 路径共用
    Vector electrolyte_div_flux(const Vector& c) const;

    /// 液相 i 格有效电导率 κ(c_e(i))·ε^bruggeman [S/m]
    ///（κ 取数值下限，防盐耗尽 c_e → 0 时欧姆项除零）
    real_t kappa_eff_at(int i) const;

    /// 液相区域恒源向量（解析电流分布，CN 与 DAE 路径共用）
    Vector dae_electrolyte_source() const;

    Parameters params_;
    real_t dae_current_ = 0.0; ///< DAE 残差的端口电流 [A]（段切换时更新）

    // ── 颗粒（固相，同 SPM）──
    int n_ = 0;
    real_t a_neg_ = 0.0;
    real_t a_pos_ = 0.0;
    real_t area_ = 0.0;
    Vector c_s_neg_;
    Vector c_s_pos_;
    SparseMatrix D_neg_;
    SparseMatrix D_pos_;
    std::vector<real_t> vol_neg_;
    std::vector<real_t> vol_pos_;

    struct CNOperators;
    std::unique_ptr<CNOperators> cn_; ///< 颗粒径向 CN 预分解
    real_t cn_dt_ = -1.0;

    // ── 电解液（液相，x 方向 1D）──
    int m_ = 0;                  ///< 厚度方向控制体数
    real_t dx_ = 0.0;            ///< 厚度方向步长 [m]
    int sep_start_ = 0;          ///< 隔膜区起始格索引
    int pos_start_ = 0;          ///< 正极区起始格索引
    Vector c_e_;                 ///< 电解液浓度 [mol/m³]
    Vector phi_e_;               ///< 液相电位分布（φ_e(0) = 0 基准）[V]
    Vector eps_e_;               ///< 各格孔隙率
    SparseMatrix le_lhs_;        ///< 电解液 CN 左端矩阵 εI − dt/2·L(c_e)
    SparseMatrix le_plus_;       ///< εI + dt/2·L(c_e)
    Eigen::SparseLU<SparseMatrix> le_solve_;

    // ── 状态 ──
    real_t voltage_ = 0.0;
    real_t eta_neg_ = 0.0;
    real_t eta_pos_ = 0.0;
    real_t eta_e_ = 0.0;
    real_t u_neg_ = 0.0;
    real_t u_pos_ = 0.0;
    real_t time_ = 0.0;
};

} // namespace fluxcell::model
