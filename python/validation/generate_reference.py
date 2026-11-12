"""生成 PyBaMM DFN 基准数据。

使用 PyBaMM 的 DFN (Doyle-Fuller-Newman) 模型，在多个 C-rate 下运行
恒流放电仿真，将时间/电压/电流数据导出为 CSV，作为 FluxCell 求解器的
验证基准。

参数集：Chen2020（石墨/NMC532，与 FluxCell 默认参数最接近的公开参数集）

输出文件（validation/reference_data/）：
    chen2020_dfn_discharge_{crate}C.csv
        列: time_s, voltage_V, current_A
"""

from pathlib import Path

import pybamm

# 输出目录
OUTPUT_DIR = Path(__file__).parent / "reference_data"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 待验证的 C-rate 列表
C_RATES = [0.5, 1.0, 2.0]


def run_discharge(c_rate: float, parameter_set: str = "Chen2020") -> Path:
    """运行单个 C-rate 的 DFN 恒流放电仿真并导出 CSV。

    Args:
        c_rate: 放电倍率（如 1.0 表示 1C）
        parameter_set: PyBaMM 参数集名称

    Returns:
        导出的 CSV 文件路径
    """
    # ── 模型与参数 ──
    model = pybamm.lithium_ion.DFN(name="DFN reference")
    params = pybamm.ParameterValues(parameter_set)

    # ── 工况：恒流放电至 2.5V ──
    experiment = pybamm.Experiment(
        [f"Discharge at {c_rate}C until 2.5 V"],
    )

    # ── 求解 ──
    sim = pybamm.Simulation(
        model,
        parameter_values=params,
        experiment=experiment,
        solver=pybamm.CasadiSolver(mode="safe"),
    )
    solution = sim.solve()

    # ── 提取数据 ──
    time_s = solution["Time [s]"].entries
    voltage_v = solution["Terminal voltage [V]"].entries
    current_a = solution["Current [A]"].entries

    # ── 导出 CSV ──
    csv_path = OUTPUT_DIR / f"chen2020_dfn_discharge_{c_rate:g}C.csv"
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("time_s,voltage_V,current_A\n")
        for t, v, i in zip(time_s, voltage_v, current_a):
            f.write(f"{t:.6e},{v:.9e},{i:.9e}\n")

    print(
        f"[OK] {c_rate:g}C 放电: {len(time_s)} 个时间点, "
        f"终止电压 {voltage_v[-1]:.4f} V, 耗时 {time_s[-1]:.1f} s"
    )
    return csv_path


def main() -> None:
    print(f"PyBaMM version: {pybamm.__version__}")
    print(f"Output directory: {OUTPUT_DIR}\n")

    for c_rate in C_RATES:
        run_discharge(c_rate)

    print("\nAll reference data generated successfully.")


if __name__ == "__main__":
    main()
