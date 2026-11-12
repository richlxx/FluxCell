/// @file dfn_model.hpp
/// @brief 全阶伪二维模型 (Doyle-Fuller-Newman / P2D)
///
/// 在 SPMe 基础上取消"每电极单颗粒 + 电解液解析电流"假设（参考
/// PyBaMM `pybamm.lithium_ion.DFN`）：
/// - 固相：每个 x 位置有独立颗粒，c_s(x, r) 为二维场（球坐标 FVM + CN，
///   同 SPM 算子对每列独立作用）；表面通量由局部 BV 决定
/// - 液相浓度：ε_e ∂c_e/∂t = ∂/∂x(D_e^eff ∂c_e/∂x) + (1−t₊)/F·a·i_surf
///   （CN 三对角，同 SPMe，但源项逐格）
/// - 反应分布：各格孔壁电流密度 i_surf = 2·i0·sinh(F·η/2RT)，
///   η = φ_s − φ_e − U(c_s_surf)，由 Gummel 型定点迭代 + 两个集流体
///   电位基准（δ_n、δ_p）的标量 Newton（电流守恒残差）求解
/// - 电位：液相 ∂φ_e/∂x = −i_e/κ^eff + (2RT/F)(1−t₊)·∂ln c_e/∂x，
///   固相 ∂φ_s/∂x = −i_s/σ^eff（逐格离散积分，i_e/i_s 由电流守恒确定）
/// - 端电压：V = φ_s,p(集流体) − φ_s,n(集流体基准)
#pragma once

#include <fluxcell/core/model/cell_model.hpp>
#include <fluxcell/core/model/spme_model.hpp>
#include <fluxcell/core/types.hpp>

#include <Eigen/SparseLU>
#include <algorithm>
#include <memory>
#include <vector>

namespace fluxcell::model {

/// 热模型类型（DFN 专属选项；SPM/SPMe 恒为等温）
enum class ThermalType : uint8_t {
    Isothermal, ///< 等温：温度恒定（与热耦合前的行为一致）
    Lumped,     ///< 集总热模型：单节点温度 ODE 与电化学双向耦合
                ///  （参考 PyBaMM options {"thermal": "lumped"}，Timms2021）
};

/// 熵变系数 dU/dT 查找表：对颗粒化学计量比 sto 的分段线性插值。
/// 数据源自外部测量 CSV（两列：嵌锂量 sto、dU/dT [mV/K]，加载时换算
/// 为 V/K），如 params/neg_graphite_dudt.csv 与 pos_ncm_dudt.csv
struct EntropyLut {
    std::vector<real_t> sto;  ///< 嵌锂量（化学计量比）节点，严格升序
    std::vector<real_t> dudt; ///< dU/dT [V/K]（与 sto 等长）

    bool empty() const { return sto.empty(); }
    std::size_t size() const { return sto.size(); }

    /// 分段线性插值；范围外取端点值（与 numpy.interp 一致），空表恒 0
    real_t eval(real_t x) const {
        if (sto.empty()) return 0.0;
        if (x <= sto.front()) return dudt.front();
        if (x >= sto.back()) return dudt.back();
        const auto it = std::upper_bound(sto.begin(), sto.end(), x);
        const auto i = static_cast<std::size_t>(it - sto.begin());
        const real_t t = (x - sto[i - 1]) / (sto[i] - sto[i - 1]);
        return dudt[i - 1] + t * (dudt[i] - dudt[i - 1]);
    }
};

/// 集总热模型参数（PyBaMM Chen2020 默认值口径）
struct ThermalParameters {
    real_t t_ambient = 298.15;      ///< 环境温度 [K]（Newton 冷却基准）
    real_t t_initial = 298.15;      ///< 初始电芯温度 [K]
    real_t h_total = 10.0;          ///< 总换热系数 [W/(m²·K)]
    real_t cooling_area = 0.00531;  ///< 冷却面积 [m²]
    real_t cell_volume = 2.42e-05;  ///< 电芯体积（热容/散热归一）[m³]

    // ── 各层体积比热容 ρ·c_p [J/(m³·K)] 与集流体厚度 [m] ──
    // 有效体积比热容按厚度加权（PyBaMM rho_c_p_eff，含集流体层）：
    // (ρcp_ccn·L_ccn + ρcp_n·L_n + ρcp_s·L_s + ρcp_p·L_p + ρcp_ccp·L_ccp) / L_total
    real_t rho_c_p_neg = 1657.0 * 700.0;     ///< 负极（1657 kg/m³ × 700）
    real_t rho_c_p_sep = 397.0 * 700.0;      ///< 隔膜（397 × 700）
    real_t rho_c_p_pos = 3262.0 * 700.0;     ///< 正极（3262 × 700）
    real_t rho_c_p_cc_neg = 8960.0 * 385.0;  ///< 负集流体（8960 × 385）
    real_t rho_c_p_cc_pos = 2700.0 * 897.0;  ///< 正集流体（2700 × 897）
    real_t cc_thickness_neg = 1.2e-05;       ///< 负集流体厚度 [m]
    real_t cc_thickness_pos = 1.6e-05;       ///< 正集流体厚度 [m]

    real_t entropic_change_neg = 0.0; ///< dU_n/dT [V/K]（常数 fallback）
    real_t entropic_change_pos = 0.0; ///< dU_p/dT [V/K]

    // dU/dT(sto) 查找表 [V/K]（非空时优先于上方常数；Chen2020 未提供
    // 熵热系数，默认 0；用户测量数据经 load_entropy_lut_csv 加载）
    EntropyLut entropic_lut_neg;
    EntropyLut entropic_lut_pos;
};

/// SEI 生长降解参数（机理寿命仿真；PyBaMM options
/// {"SEI": "ec reaction limited", "SEI film resistance": "distributed"}
/// 口径，默认值取 OKane2022 参数集，电化学部分沿用 Chen2020）
struct DegradationParameters {
    bool enabled = false;          ///< 启用负极 SEI 生长降解
    real_t l_sei_0 = 5e-9;         ///< 初始 SEI 厚度 [m]
    real_t v_bar_sei = 9.585e-05;  ///< SEI 偏摩尔体积 [m³/mol]
    real_t z_sei = 1.0;            ///< 锂/SEI 摩尔比（OKane2022 = 1）
    real_t u_sei = 0.4;            ///< SEI 开路电位 [V]
    real_t r_sei = 2.0e5;          ///< SEI 电阻率 [Ohm·m]
    real_t k_sei = 1.0e-12;        ///< SEI 反应速率常数 [m/s]
    real_t c_ec_0 = 4541.0;        ///< 电解液 EC 初始浓度 [mol/m³]
    real_t d_ec = 2.0e-18;         ///< EC 在 SEI 中的扩散系数 [m²/s]
    real_t e_sei = 38000.0;        ///< SEI 生长活化能 [J/mol]
    real_t t_ref = 298.15;         ///< Arrhenius 参考温度 [K]
    real_t alpha_sei = 0.5;        ///< SEI 反应传递系数（对称）
};

/// 全阶 DFN (P2D) 模型
class DFNModel : public CellModel, public solver::DaeProblem {
public:
    /// DFN 参数 = SPMe 参数（固相 + 电解液）+ 固相电导/迭代控制 + 热模型
    struct Parameters : SPMeModel::Parameters {
        real_t neg_conductivity = 215.0;  ///< σ_s,n [S/m]（Chen2020 石墨）
        real_t pos_conductivity = 0.18;   ///< σ_s,p [S/m]（Chen2020 NMC）
        int max_iterations = 100;         ///< 反应分布定点迭代上限
        real_t tolerance = 1e-9;          ///< 电流守恒/自洽相对容差
        real_t damping = 0.5;             ///< 反应电流场阻尼因子（0~1）

        ThermalType thermal_type = ThermalType::Isothermal;
        ThermalParameters thermal_params; ///< 热模型参数（等温模式忽略）

        DegradationParameters degradation; ///< SEI 降解（enabled = false 时全路径零开销）
    };

    /// @param params 模型参数
    /// @param initial_current 初始工作电流 [A]（正 = 放电；t = 0 过电位）
    explicit DFNModel(const Parameters& params, real_t initial_current = 0.0);
    ~DFNModel() override;

    void step(real_t current, real_t dt) override;
    real_t voltage() const override { return voltage_; }
    real_t time() const override { return time_; }
    SolutionSnapshot snapshot(real_t current) const override;

    /// 各 x 格颗粒表面浓度 [mol/m³]（二阶外推，列向量）
    const Vector& neg_surface_conc_profile() const { return c_surf_n_; }
    const Vector& pos_surface_conc_profile() const { return c_surf_p_; }

    /// 各 x 格颗粒表面（孔壁）电流密度 [A/m²]（正 = 放电方向）
    const Vector& neg_reaction_current() const { return i_surf_n_; }
    const Vector& pos_reaction_current() const { return i_surf_p_; }

    /// 电解液浓度分布 [mol/m³] / 液相电位分布（φ_e(0) = 0 基准）[V]
    const Vector& electrolyte_conc() const { return c_e_; }
    const Vector& electrolyte_potential() const { return phi_e_; }
    const Vector& solid_potential() const { return phi_s_full_; }

    /// 液相过电位（欧姆 + 浓差）[V]
    real_t electrolyte_overpotential() const { return eta_e_; }

    /// 当前电芯温度 [K]（等温模式恒为参数温度；集总模式为状态变量）
    real_t temperature() const { return t_cell_; }

    /// 最近一次产热功率 [W]：q_total = q_ohm + q_rxn + q_rev（PyBaMM 口径）
    real_t total_heating() const { return q_total_; }
    real_t ohmic_heating() const { return q_ohm_; }
    real_t reaction_heating() const { return q_rxn_; }
    real_t reversible_heating() const { return q_rev_; }

    // ── 机理寿命诊断（degradation.enabled = false 时恒 0/false）──
    bool degradation_enabled() const { return params_.degradation.enabled; }

    /// 负极 SEI 平均厚度（x 方向算术平均）[m]
    real_t sei_thickness_avg() const;

    /// SEI 生长累计锂损失 [Ah]（PyBaMM "Loss of capacity to negative SEI"）
    real_t lithium_loss_ah() const;

    /// 含降解的估计可用容量 [Ah]（标称容量 − 累计锂损失；无 LAM 时的
    /// 锂库存口径，与 PyBaMM LLCI 单一机理一致）
    real_t degraded_capacity_ah() const;

    /// 固相（各列颗粒）+ 液相锂总量 [mol] — 守恒性测试用
    real_t total_lithium_mol() const;

    const Parameters& parameters() const { return params_; }

    real_t nominal_capacity() const override { return params_.solid.capacity; }
    real_t initial_soc() const override { return params_.solid.initial_soc; }

    // ── solver::DaeProblem 接口（指数-1 DAE，增广外电路端口）──
    // 状态布局：[c_s_neg (n×m_neg 列主), c_s_pos, c_e, i_surf_n, i_surf_p,
    //           δ_n, δ_p, (T), (L_sei...), I_cell]；
    // 微分变量 = 浓度场 + 温度/SEI（可选），代数变量 = i_surf、δ 与端口
    // 电流 I_cell（末位）。外电路方程为残差末行（参考
    // pybamm.external_circuit）：CC 锁 I_cell − I_target，CV 锁
    // V(y) − V_target——端口电流作为 DAE 代数变量由 Newton 一并求解，
    // CV 段无需引擎侧外层反解。
    // lumped 热模型在 δ 之后追加温度 t_cell_（微分变量，off_t = off_d + 2）；
    // SEI 降解在温度之后追加 m_neg 个 SEI 厚度状态（微分变量，
    // off_sei = off_d + 2 + (lumped ? 1 : 0)）；I_cell 恒在末位
    // （off_i = dae_size() − 1）。
    // 代数约束：逐格 BV 方程（含 SEI 膜电阻压降 −(j+j_sei)·L·R_sei；
    // φ 场由 i_surf/δ 解析积分）+ 两电极电流守恒（含 SEI 副反应分流，
    // i_app = I_cell/area）+ 外电路行。
    // BV 反应分布由 IDA 的 Newton 与浓度场一致求解，
    // 替代固定步长路径的 Gummel 定点迭代 + 显式子步细分。
    int dae_size() const override;
    void dae_mask(Vector& id) const override;
    void dae_atol(Vector& atol) const override;
    void dae_pack_state(Vector& y, Vector& yp) const override;
    void dae_residual(real_t t, const Eigen::Ref<const Vector>& y,
                      const Eigen::Ref<const Vector>& yp,
                      Eigen::Ref<Vector> res) override;
    void dae_load_state(real_t t,
                                const Eigen::Ref<const Vector>& y) override;
    /// 电流切换时除记录外，还用半解析 Gummel 迭代重估代数变量
    ///（i_surf 分布与 δ），为 IDA reinit 提供自洽初值
    void dae_set_current(real_t current) override;
    /// 切换恒压控制：外电路行 = V(y) − voltage；代数变量按当前端口
    /// 电流重估自洽初值（CC→CV 切换时端口电流继承 CC 段末值，端电压
    /// 已 ≈ 目标，CalcIC 最小修正收敛）
    void dae_set_voltage(real_t voltage) override;
    bool dae_port_augmented() const override { return true; }
    real_t dae_port_current() const override { return i_port_; }
    real_t dae_current_target() const override { return dae_current_; }
    /// 推进失败时的状态摘要（stderr）：各场区 min/max，定位 corrector 失败区
    void dae_dump_state(real_t t,
                        const Eigen::Ref<const Vector>& y) const override;
    /// SEI 厚度列裁剪为仅 SEI 行对角（移除 j_sei ∝ 1/L 链差分巨元）
    void dae_jac_col_prune(int j, std::vector<int>& rows) const override;
    solver::DaeProblem* dae_problem() override { return this; }

private:
    /// 构建/重建 Crank-Nicolson 算子（颗粒每列独立 + 液相，dt 变化时重建）
    void ensure_cn(real_t dt);

    /// 有效体积比热容 ρc_p^eff [J/(m³·K)]（各层按厚度加权，含集流体）
    real_t rho_c_p_eff() const;

    /// 集总热容 C_th = ρc_p^eff · V_cell [J/K]
    real_t thermal_capacitance() const;

    /// 当前状态瞬时产热功率（PyBaMM base_thermal 三项口径）：
    /// - Q_ohm：固相 i_s²/σ^eff + 液相 −i_e·∂φ_e（含浓差扩散项）
    /// - Q_rxn：不可逆界面热 a·i_surf·η
    /// - Q_rev：可逆熵热 a·i_surf·T·dU/dT
    /// 成员场（i_e_face_/φ/η/i_surf）必须已是自洽解；结果缓存到 q_*_ 成员
    /// @param i_app 端口面电流密度 [A/m²]（固相电流 = i_app − i_e_face）
    /// @return 总产热 [W]
    real_t refresh_heat(real_t i_app);

    /// 液相守恒通量散度 L(c)·c：D^eff = D(c_e)·ε^bruggeman 逐格、
    /// 界面算术平均，两端零通量；浓度相关 D 的 CN/DAE 路径共用
    Vector electrolyte_div_flux(const Vector& c) const;

    /// 液相 i 格有效电导率 κ(c_e(i))·ε^bruggeman [S/m]
    ///（κ 取数值下限，防盐耗尽 c_e → 0 时欧姆项除零）
    real_t kappa_eff_at(const Eigen::Ref<const Vector>& c_e, int i) const;

    /// 从 c_s 各列提取表面浓度、嵌入率对应的 OCP（迭代与快照共用）
    void refresh_surface_quantities();

    /// 求解代数反应分布：Gummel 定点迭代（i_surf ↔ φ 场）+
    /// δ_n/δ_p 标量 Newton（电流守恒）。更新 i_surf/φ 场/eta 分布与 voltage_。
    /// @param i_app 端口面电流密度 [A/m²]（正 = 放电）
    void solve_reaction_distribution(real_t i_app);

    /// 由当前 i_surf 场积分液相/固相电流与电位分布（i_e_face/φ_e/φ_s/φ_s_full）
    void integrate_fields(real_t i_app);

    /// 自洽收尾：积分 φ 场并更新过电位分布/液相过电位/端电压
    ///（Gummel 收敛后与 DAE load_state 共用）
    void finalize_fields(real_t i_app);

    /// 估计表面浓度–BV 局部回路的电化学弛豫率 λ（i_bv 经 i0 与 U 依赖
    /// c_surf，CN 表面源为显式处理），按子步长 ≤ 0.4/λ 返回 dt 的子步数
    int plan_substeps(real_t dt) const;

    /// 单个子步推进（step 的主体：解代数分布 → CN 推进浓度 → 重解分布）
    void step_once(real_t current, real_t dt);

    // ── SEI 降解（degradation.enabled 时生效）──

    /// 负极 SEI 副反应电流（ec reaction limited，PyBaMM sei_growth.py
    /// Yang2017 口径）：j_sei = −F·c_ec0·K·arrh/(1 + L·K/D_ec)，
    /// K = k_sei·exp(−α·F·eta_SEI/RT)，eta_SEI = (φ_s−φ_e) − U_sei − j·L·R_sei，
    /// arrh 乘整体（与 PyBaMM 一致，分母不含）。恒负 = 锂消耗。
    /// 同时输出 dL/dt = g/(1 + q·L) 的 g [m/s] 与 q [1/m]（CN 隐式推进用）
    /// @param k 负极格索引（L 取 l_sei_(k)）
    /// @param L 本格 SEI 厚度 [m]
    /// @param eta_bare φ_s − φ_e [V]（不含 U_sei）
    real_t sei_kinetics_neg(int k, real_t L, real_t eta_bare,
                            real_t i_surf_k, real_t T,
                            real_t& g, real_t& q) const;

    /// 用当前成员状态刷新 j_sei_n_/eta_sei_n_（构造、dae_load_state、
    /// SEI 推进后同步副反应电流与 SEI 过电位诊断）
    void refresh_sei_fields();

    /// SEI 厚度推进（step_once 步末调用）：φ/T/i_surf 用步末自洽解冻结，
    /// 逐格完全隐式一元二次根（q→0 退化为显式），随后刷新 j_sei_n_
    void advance_sei_cn(real_t dt);

    Parameters params_;

    // ── x 方向网格（同 SPMe）──
    int m_ = 0;            ///< 厚度方向控制体数
    int n_ = 0;            ///< 颗粒径向控制体数
    real_t dx_ = 0.0;      ///< 厚度方向步长 [m]
    int sep_start_ = 0;    ///< 隔膜区起始格索引
    int pos_start_ = 0;    ///< 正极区起始格索引
    int m_neg_ = 0;        ///< 负极格数 = sep_start_
    int m_pos_ = 0;        ///< 正极格数 = m_ - pos_start_
    Vector eps_e_;         ///< 各格孔隙率
    real_t sigma_n_eff_ = 0.0; ///< 负极固相有效电导率 [S/m]
    real_t sigma_p_eff_ = 0.0; ///< 正极固相有效电导率 [S/m]
    real_t a_n_ = 0.0;     ///< 负极比表面积 a = 3·ε_s/R [m²/m³]
    real_t a_p_ = 0.0;     ///< 正极比表面积 [m²/m³]
    real_t area_ = 0.0;    ///< 电极面积 A [m²]

    // ── 固相颗粒（2D：列 = 各 x 位置的颗粒）──
    Eigen::MatrixXd c_s_neg_;  ///< n_ × m_neg_ 径向浓度 [mol/m³]
    Eigen::MatrixXd c_s_pos_;  ///< n_ × m_pos_
    SparseMatrix D_neg_;       ///< 球坐标径向扩散算子（同 SPM）
    SparseMatrix D_pos_;
    SparseMatrix cn_plus_neg_; ///< I + dt/2·D
    SparseMatrix cn_plus_pos_;
    Eigen::SparseLU<SparseMatrix> cn_solve_neg_; ///< (I − dt/2·D)⁻¹
    Eigen::SparseLU<SparseMatrix> cn_solve_pos_;
    real_t surf_area_n_ = 0.0; ///< 颗粒外表面面积 [m²]
    real_t surf_vol_n_ = 0.0;  ///< 表面控制体体积 [m³]
    real_t surf_area_p_ = 0.0;
    real_t surf_vol_p_ = 0.0;
    std::vector<real_t> vol_neg_; ///< 径向控制体体积（守恒用）
    std::vector<real_t> vol_pos_;
    real_t cn_dt_ = -1.0;

    // ── 液相 CN（同 SPMe，D_e 浓度相关每步重建）──
    SparseMatrix le_lhs_;  ///< εI − dt/2·L(c_e)
    SparseMatrix le_plus_; ///< εI + dt/2·L(c_e)
    Eigen::SparseLU<SparseMatrix> le_solve_;

    // ── 场变量 ──
    Vector c_e_;        ///< 电解液浓度 [mol/m³]
    Vector phi_e_;      ///< 液相电位（φ_e(0) = 0 基准）[V]
    Vector phi_s_n_;    ///< 负极固相电位（基准 φ_s_n(0) = δ_n）[V]
    Vector phi_s_p_;    ///< 正极固相电位（基准 φ_s_p(0) = δ_p）[V]
    Vector phi_s_full_; ///< 全厚度固相电位（隔膜区延伸，快照用）[V]
    Vector i_e_face_;   ///< 液相面电流密度（格界面，m_+1）[A/m²]
    Vector i_surf_n_;   ///< 负极孔壁电流密度 [A/m²]
    Vector i_surf_p_;   ///< 正极孔壁电流密度 [A/m²]
    Vector c_surf_n_;   ///< 负极各格表面浓度 [mol/m³]
    Vector c_surf_p_;   ///< 正极各格表面浓度 [mol/m³]
    Vector u_n_vec_;    ///< 负极各格 OCP [V]
    Vector u_p_vec_;    ///< 正极各格 OCP [V]
    Vector eta_n_;      ///< 负极过电位分布 [V]
    Vector eta_p_;      ///< 正极过电位分布 [V]

    // ── 标量状态 ──
    real_t delta_n_ = 0.0;   ///< 负极集流体固相电位（迭代变量）[V]
    real_t delta_p_ = 0.0;   ///< 正极/隔膜界面固相电位（迭代变量）[V]
    real_t voltage_ = 0.0;
    real_t eta_e_ = 0.0;
    real_t time_ = 0.0;
    real_t dae_current_ = 0.0; ///< 恒流外电路行的端口电流目标 [A]（段切换
                               ///  时更新；Gummel 重估与短路比较基准）
    real_t i_port_ = 0.0;      ///< 最近一致化的端口电流 [A]（状态中的
                               ///  I_cell，dae_load_state 同步；dae_port_current）
    real_t v_target_ = 0.0;    ///< 恒压外电路行的端电压目标 [V]
    solver::PortControl port_control_ = solver::PortControl::Current;
                               ///< 外电路控制模式（CC/CV；切换后须 reinit）

    // ── 热状态（集总热模型；等温模式 t_cell_ 恒为参数温度）──
    real_t t_cell_ = 298.15;     ///< 电芯温度 [K]（lumped 模式为状态变量）
    real_t rho_c_p_eff_ = 0.0;   ///< 有效体积比热容缓存 [J/(m³·K)]
    real_t q_total_ = 0.0;       ///< 最近产热缓存 [W]（诊断/快照）
    real_t q_ohm_ = 0.0;
    real_t q_rxn_ = 0.0;
    real_t q_rev_ = 0.0;

    // ── SEI 降解状态（degradation.enabled = false 时为空）──
    Vector l_sei_;     ///< 负极各格 SEI 厚度 [m]
    Vector j_sei_n_;   ///< 负极 SEI 副反应电流密度 [A/m²]（恒负 = 锂消耗）
    Vector eta_sei_n_; ///< 负极 SEI 反应过电位分布 [V]（产热诊断）
};

} // namespace fluxcell::model
