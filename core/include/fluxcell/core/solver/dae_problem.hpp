/// @file dae_problem.hpp
/// @brief DAE 问题抽象：模型向自适应求解器（SUNDIALS IDA）描述微分-代数方程
///
/// 纯 C++ 接口，不依赖 SUNDIALS：微分方程部分 y' = f(t, y)，代数约束部分
/// 0 = g(t, y)。状态向量布局由具体模型定义，微分/代数变量通过掩码区分。
/// SPM / SPMe 为纯 ODE（无代数行）；DFN 将逐格 BV 动力学与电流守恒
/// 作为代数约束（指数-1 DAE），由求解器 Newton 一致求解，替代固定步长
/// 路径的显式子步细分。
#pragma once

#include <vector>

#include <fluxcell/core/types.hpp>

namespace fluxcell::solver {

/// 外电路控制模式（增广端口的代数方程选择）
enum class PortControl : uint8_t {
    Current, ///< 恒流：外电路行 = I_cell − I_target
    Voltage, ///< 恒压：外电路行 = V(y) − V_target（PyBaMM
             ///  external_circuit.ConstantVoltage 同构）
};

/// DAE 问题接口（由支持自适应求解的模型实现）
class DaeProblem {
public:
    virtual ~DaeProblem() = default;

    /// 状态向量维度（微分 + 代数变量总数）
    virtual int dae_size() const = 0;

    /// 变量类型掩码：1.0 = 微分变量（参与误差控制），0.0 = 代数变量
    virtual void dae_mask(Vector& id) const = 0;

    /// 绝对误差容限向量（按各变量物理尺度给出）
    virtual void dae_atol(Vector& atol) const = 0;

    /// 由模型当前状态打包 y 与 y'（微分变量的导数 = f(y)，代数变量导数 = 0）
    virtual void dae_pack_state(Vector& y, Vector& yp) const = 0;

    /// 残差 F(t, y, y')：微分行 = y' − f(t, y)，代数行 = g(t, y)
    ///（Ref 签名：接受连续存储的 Vector 或零拷贝 Map 视图）
    virtual void dae_residual(real_t t, const Eigen::Ref<const Vector>& y,
                              const Eigen::Ref<const Vector>& yp,
                              Eigen::Ref<Vector> res) = 0;

    /// 接受求解器状态：y → 模型内部场，并重算电压/过电位等导出量
    virtual void dae_load_state(real_t t,
                                const Eigen::Ref<const Vector>& y) = 0;

    /// 更新残差依赖的端口电流 [A]（正 = 放电）。
    /// 电流突变后残差不连续，调用方必须随后对求解器 reinit。
    /// 含代数变量的实现（如 DFN）应在切换时将代数变量重估为
    /// 与新电流自洽的值，为重初始化的一致化提供良好初值。
    virtual void dae_set_current(real_t current) = 0;

    // ── 增广外电路端口（可选；PyBaMM external_circuit 同构）──
    // 模型将端口电流 I_cell 作为状态末位的代数变量，外电路方程
    // 作为残差末行由求解器 Newton 一并求解：CC 锁电流、CV 锁端
    // 电压——CV 段无需引擎侧外层反解。默认关闭（纯 ODE 模型
    // /未增广实现），CV 由引擎侧反馈反解承担。

    /// 模型是否已增广端口电流状态（支持 CC/CV 代数外电路行）
    virtual bool dae_port_augmented() const { return false; }

    /// 切换恒压控制：外电路行 = V(y) − voltage = 0；代数变量按
    /// 当前端口电流重估自洽初值。切换后须 reinit 生效
    virtual void dae_set_voltage(real_t voltage) {
        (void)voltage;
    }

    /// 最近一致化的端口电流 [A]（增广模型 = 状态中的 I_cell）
    virtual real_t dae_port_current() const { return 0.0; }

    /// 恒流外电路行的电流目标 [A]（dae_set_current 最近设置的值；
    /// 供求解器探测外电路两种行形态后恢复原模式用）
    virtual real_t dae_current_target() const { return 0.0; }

    /// 推进失败时的状态摘要诊断（定位 corrector 失败的状态区）；
    /// 默认无输出，由具体模型按自身布局覆写
    virtual void dae_dump_state(real_t t,
                                const Eigen::Ref<const Vector>& y) const {
        (void)t;
        (void)y;
    }

    /// 差分 Jacobian 列支持集裁剪（可选，探测完成后、着色前调用）：
    /// 把该列探测到的行支持集替换为更小的结构集（默认不裁剪）。
    /// 用于移除「经中间变量链差分出的病态巨元」：如 SEI 厚度 L 在
    /// K 饱和区 j_sei ∝ 1/L，经守恒行/BV 行/热行反应热差分出 ~1e4-3e5
    /// 元素，而热行对角在 h 大时仅 ~1e-3 → 行严重对角劣势 → KLU
    /// pivot 失效（corrector 全灭，lumped+SEI 长循环 CV 段 flag -4
    /// 实证）。被删反馈对 Newton 的实际影响 ~∂j_sei/∂L×ΔL
    /// （dz_L ~ 1e-11）≈ 高阶小量，且 j_sei↔φ 回路已在残差内不动点
    /// 精确解出，Jacobian 以裁剪后的 quasi-Newton 近似收敛
    virtual void dae_jac_col_prune(int j, std::vector<int>& rows) const {
        (void)j;
        (void)rows;
    }
};

} // namespace fluxcell::solver
