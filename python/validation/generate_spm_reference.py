"""生成 PyBaMM SPM 基准数据。

使用 PyBaMM 的 SPM (Single Particle Model)，在多个 C-rate 下运行
恒流放电仿真，将时间/电压/电流数据导出为 CSV，作为 FluxCell SPM
求解器的验证基准。

参数集：Chen2020（石墨/NMC532）

输出文件（validation/reference_data/）：
    chen2020_spm_discharge_{crate}C.csv
        列: time_s, voltage_V, current_A

用法：
    uv run python validation/generate_spm_reference.py
"""

from pathlib import Path

import pybamm

# 输出目录
OUTPUT_DIR = Path(__file__).parent / "reference_data"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 待验证的 C-rate 列表
C_RATES = [0.5, 1.0, 2.0]


def run_discharge(c_rate: float, parameter_set: str = "Chen2020") -> Path:
    """运行单个 C-rate 的 SPM 恒流放电仿真并导出 CSV。

    Args:
        c_rate: 放电倍率（如 1.0 表示 1C）
        parameter_set: PyBaMM 参数集名称

    Returns:
        导出的 CSV 文件路径
    """
    # ── 模型与参数 ──
    model = pybamm.lithium_ion.SPM(name="SPM reference")
    params = pybamm.ParameterValues(parameter_set)

    # ── 工况：恒流放电至 2.5V（period=1s 保证初始瞬态段有足够采样密度，
    #    否则与 FluxCell 等间隔结果对比时会产生插值伪影）──
    experiment = pybamm.Experiment(
        [f"Discharge at {c_rate}C until 2.5 V"],
        period="1s",
    )

    # ── 求解 ──
    sim = pybamm.Simulation(
        model,
        parameter_values=params,
        experiment=experiment,
        solver=pybamm.CasadiSolver(mode="safe"),
    )
    solution = sim.solve()

    # ── 提取数据（附关键内部变量用于精细对齐）──
    def flat(name: str):
        """提取变量并展平为 1D 数组（部分变量形状为 (1, N)）。"""
        import numpy as np
        return np.ravel(solution[name].entries)
    
    time_s = flat("Time [s]")
    voltage_v = flat("Terminal voltage [V]")
    current_a = flat("Current [A]")
    u_neg = flat("Negative electrode open-circuit potential [V]")
    u_pos = flat("Positive electrode open-circuit potential [V]")
    eta_neg = flat("Negative electrode reaction overpotential [V]")
    eta_pos = flat("Positive electrode reaction overpotential [V]")

    # ── 导出 CSV ──
    csv_path = OUTPUT_DIR / f"chen2020_spm_discharge_{c_rate:g}C.csv"
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("time_s,voltage_V,current_A,u_neg_V,u_pos_V,eta_neg_V,eta_pos_V\n")
        for t, v, i, un, up, en, ep in zip(
            time_s, voltage_v, current_a, u_neg, u_pos, eta_neg, eta_pos
        ):
            f.write(f"{t:.6e},{v:.9e},{i:.9e},{un:.9e},{up:.9e},{en:.9e},{ep:.9e}\n")

    print(
        f"[OK] {c_rate:g}C SPM 放电: {len(time_s)} 个时间点, "
        f"终止电压 {voltage_v[-1]:.4f} V, 耗时 {time_s[-1]:.1f} s"
    )
    return csv_path


def main() -> None:
    print(f"PyBaMM version: {pybamm.__version__}")
    print(f"Output directory: {OUTPUT_DIR}\n")

    for c_rate in C_RATES:
        run_discharge(c_rate)

    print("\nAll SPM reference data generated successfully.")


if __name__ == "__main__":
    main()
