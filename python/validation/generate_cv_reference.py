"""生成 PyBaMM CC-CV 充电基准数据（对标 FluxCell 的 DFN(IDA) CC-CV 工况）。

工况（复现 GUI 报告场景）：
    1. Charge at 1C until 4.2 V   （恒流充电，4.2 V 截止）
    2. Hold at 4.2 V until C/50   （恒压充电，尾电流 C/50 = 0.1 A 截止）

参数集：Chen2020（LG M50, 5 Ah），初始 SOC = 0%（用户场景：从空电完整充电）。

输出文件（validation/reference_data/）：
    chen2020_dfn_cccv_1C.csv
        列: time_s, voltage_V, current_A, soc（库仑计数，soc = 1 + ∫I dt/(3600·Q)）

用法：
    uv run python validation/generate_cv_reference.py
"""

from pathlib import Path

import numpy as np
import pybamm

OUTPUT_DIR = Path(__file__).parent / "reference_data"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

INITIAL_SOC = 0.0  # 与 sim.solve(initial_soc=...) 保持一致


def main() -> None:
    model = pybamm.lithium_ion.DFN(name="DFN CC-CV reference")
    params = pybamm.ParameterValues("Chen2020")

    experiment = pybamm.Experiment(
        [
            "Charge at 1C until 4.2 V",
            "Hold at 4.2 V until C/50",
        ],
        period="1s",
    )

    sim = pybamm.Simulation(
        model,
        parameter_values=params,
        experiment=experiment,
        solver=pybamm.CasadiSolver(mode="safe"),
    )
    # initial_soc=0：由 0%/100% OCV 端点反推初始化学计量
    # （与 FluxCell 的 neg/pos_stoich_at_soc0 同口径）
    solution = sim.solve(initial_soc=0.0)

    t = solution["Time [s]"].entries
    v = solution["Terminal voltage [V]"].entries
    i = solution["Current [A]"].entries
    q_nom = params["Nominal cell capacity [A.h]"]

    # 库仑计数 SOC（与 FluxCell 引擎同口径）：放电 +I 下降、充电 -I 回升，
    # 起点为 INITIAL_SOC
    dt = np.diff(t, prepend=t[0])
    soc = INITIAL_SOC - np.cumsum(i * dt) / 3600.0 / q_nom

    csv_path = OUTPUT_DIR / "chen2020_dfn_cccv_1C.csv"
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("time_s,voltage_V,current_A,soc\n")
        for ti, vi, ii, si in zip(t, v, i, soc):
            f.write(f"{ti:.6e},{vi:.9e},{ii:.9e},{si:.9e}\n")

    print(f"steps: {len(t)}  t_end: {t[-1]:.0f} s  "
          f"V_end: {v[-1]:.4f} V  I_end: {i[-1]:.4f} A  "
          f"SOC_end: {soc[-1]*100:.2f}%")
    print(f"-> {csv_path}")

    # 关键节点摘要（与 FluxCell 诊断输出同粒度，每 120 s）
    print("\n---- 每 120 s 摘要 ----")
    for idx in range(0, len(t), 120):
        print(f"t={t[idx]:7.0f} s  V={v[idx]:.4f}  I={i[idx]:+7.3f} A  "
              f"SOC={soc[idx]*100:6.2f}%")


if __name__ == "__main__":
    main()
