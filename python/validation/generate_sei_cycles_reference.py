"""生成 PyBaMM DFN + SEI 降解长循环基准数据（圈数/工况可调）。

在 generate_degradation_reference.py 的 2 圈基准基础上扩展：
圈数 N 与 CCCV 工况参数化，用于长循环寿命对标（如 50 圈）与
FluxCell 循环圈数工况（含恒压阶段）的基准生成。

示例：
    python generate_sei_cycles_reference.py --cycles 50
    python generate_sei_cycles_reference.py --cycles 20 --cv-voltage 4.4

输出文件（validation/reference_data/）：
    okane2022_dfn_sei_cc{rate}c_{N}cycles.csv
    okane2022_dfn_sei_cccv{V}v_{rate}c_{N}cycles.csv
        列: time_s, voltage_V, li_loss_ah, sei_thickness_m, cycle_end
"""

import argparse
from pathlib import Path

import numpy as np
import pybamm

OUTPUT_DIR = Path(__file__).parent / "reference_data"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cycles", type=int, default=2, help="循环圈数")
    parser.add_argument("--rate", type=float, default=1.0, help="充放电倍率 [C]")
    parser.add_argument("--v-high", type=float, default=4.2,
                        help="充电截止/恒压电压 [V]")
    parser.add_argument("--v-low", type=float, default=2.5,
                        help="放电截止电压 [V]")
    parser.add_argument("--cv-voltage", type=float, default=None,
                        help="启用 CCCV：CC 充电到该电压后恒压（默认不启用）")
    parser.add_argument("--sample-dt", type=float, default=60.0,
                        help="采样间隔 [s]")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"PyBaMM version: {pybamm.__version__}")
    print(f"cycles={args.cycles} rate={args.rate}C "
          f"V:[{args.v_low}, {args.v_high}] cv={args.cv_voltage}")

    # ── 模型：DFN + ec reaction limited SEI + 分布式膜电阻 + 等温 ──
    model = pybamm.lithium_ion.DFN(
        options={
            "SEI": "ec reaction limited",
            "SEI film resistance": "distributed",
            "thermal": "isothermal",
        },
        name="DFN SEI degradation long cycles",
    )
    params = pybamm.ParameterValues("OKane2022")

    # ── 工况：1C 放电 → 1C 充电（→ 恒压）──
    discharge = f"Discharge at {args.rate}C until {args.v_low} V"
    charge = f"Charge at {args.rate}C until {args.v_high} V"
    steps = [discharge, charge]
    if args.cv_voltage is not None:
        steps.append(f"Hold at {args.cv_voltage} V until C/50")
    experiment = pybamm.Experiment([tuple(steps)] * args.cycles)

    sim = pybamm.Simulation(
        model,
        parameter_values=params,
        experiment=experiment,
        solver=pybamm.CasadiSolver(mode="safe"),
    )
    solution = sim.solve()

    # ── 段/圈末诊断 ──
    for i, cyc in enumerate(solution.cycles):
        ends = []
        for j, st in enumerate(cyc.steps):
            ends.append(
                f"st{j} t={st.t[-1] - st.t[0]:.0f}s "
                f"V={st['Terminal voltage [V]'].entries[-1]:.3f}"
            )
        print(f"  cycle {i}: " + " | ".join(ends))

    # ── 均匀采样 ──
    t_eval = np.arange(0.0, solution.t[-1] + args.sample_dt, args.sample_dt)
    voltage = solution["Terminal voltage [V]"](t_eval)
    li_loss = solution["Loss of capacity to negative SEI [A.h]"](t_eval)
    sei_thick = solution["X-averaged negative SEI thickness [m]"](t_eval)

    # 圈末绝对时刻（cycle.t 为绝对时间）
    cycle_end_times = [float(cyc.t[-1]) for cyc in solution.cycles]

    if args.cv_voltage is not None:
        csv_path = OUTPUT_DIR / (
            f"okane2022_dfn_sei_cccv{args.cv_voltage:g}v_"
            f"{args.rate:g}c_{args.cycles}cycles.csv"
        )
    else:
        csv_path = OUTPUT_DIR / (
            f"okane2022_dfn_sei_cc{args.rate:g}c_{args.cycles}cycles.csv"
        )
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("time_s,voltage_V,li_loss_ah,sei_thickness_m,cycle_end\n")
        for i, t in enumerate(t_eval):
            is_end = 1 if any(
                t <= te < t + args.sample_dt for te in cycle_end_times
            ) else 0
            f.write(
                f"{t:.6e},{voltage[i]:.9e},{li_loss[i]:.9e},"
                f"{sei_thick[i]:.9e},{is_end}\n"
            )

    print(f"\n[OK] {len(t_eval)} 个采样点, 总时长 {t_eval[-1]:.0f} s")
    for i, te in enumerate(cycle_end_times):
        print(
            f"     圈 {i + 1} 末: t={te:.0f} s  "
            f"V={solution['Terminal voltage [V]'](te):.4f}  "
            f"loss={solution['Loss of capacity to negative SEI [A.h]'](te):.5f} Ah  "
            f"L_sei={solution['X-averaged negative SEI thickness [m]'](te):.4e} m"
        )
    print(f"\nWritten: {csv_path}")


if __name__ == "__main__":
    main()
