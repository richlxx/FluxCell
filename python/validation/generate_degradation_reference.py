"""生成 PyBaMM DFN + SEI 降解（循环工况）基准数据。

使用 PyBaMM 的 DFN 模型（options {"SEI": "ec reaction limited",
"SEI film resistance": "distributed", "thermal": "isothermal"}，
OKane2022 参数集）运行 2 圈 1C 放电/充电循环，导出时间/电压/
SEI 锂损失/SEI 厚度数据，作为 FluxCell 机理寿命仿真的验证基准。

口径说明：
    - OKane2022 = Chen2020 电化学参数 + SEI/plating/LAM 降解参数；
      等温 T=298.15 时其 Arrhenius 输运退化为常数基准，固相
      Bruggeman 1.5 与 FluxCell σ_eff = σ·ε^1.5 口径一致，
      因此 FluxCell 侧可直接用 params/dfn_chen2020.yml（电化学
      Chen2020 + degradation 段 SEI 参数取 OKane2022）对标
    - 每圈末（cycle_end）快照为寿命曲线采样点

输出文件（validation/reference_data/）：
    okane2022_dfn_sei_cc_2cycles.csv
        列: time_s, voltage_V, li_loss_ah, sei_thickness_m, cycle_end
"""

from pathlib import Path

import numpy as np
import pybamm

# 输出目录
OUTPUT_DIR = Path(__file__).parent / "reference_data"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 循环圈数与采样间隔 [s]
N_CYCLES = 2
SAMPLE_DT = 30.0


def main() -> None:
    print(f"PyBaMM version: {pybamm.__version__}")
    print(f"Output directory: {OUTPUT_DIR}\n")

    # ── 模型：DFN + ec reaction limited SEI + 分布式膜电阻 + 等温 ──
    model = pybamm.lithium_ion.DFN(
        options={
            "SEI": "ec reaction limited",
            "SEI film resistance": "distributed",
            "thermal": "isothermal",
        },
        name="DFN SEI degradation reference",
    )
    params = pybamm.ParameterValues("OKane2022")

    # ── 工况：1C 放电至 2.5V + 1C 充电至 4.2V，重复 N_CYCLES 圈 ──
    experiment = pybamm.Experiment(
        [
            ("Discharge at 1C until 2.5 V", "Charge at 1C until 4.2 V"),
        ]
        * N_CYCLES,
    )

    sim = pybamm.Simulation(
        model,
        parameter_values=params,
        experiment=experiment,
        solver=pybamm.CasadiSolver(mode="safe"),
    )
    solution = sim.solve()

    # ── 圈末时刻（每圈最后一段的截止时刻，寿命曲线采样点）──
    # 先诊断每段时长与末态（确认段截止事件都正常触发）
    for i, cyc in enumerate(solution.cycles):
        for j, st in enumerate(cyc.steps):
            print(
                f"  cycle {i} step {j}: t_end={st.t[-1]:.2f}  "
                f"V={st['Terminal voltage [V]'].entries[-1]:.4f}  "
                f"I={st['Current [A]'].entries[-1]:.3f}"
            )
    for i, cyc in enumerate(solution.cycles):
        print(
            f"  cycle {i}: t[0]={cyc.t[0]:.2f}  t[-1]={cyc.t[-1]:.2f}"
        )
    # cycle.t 为绝对时间（t[0] = 圈首绝对时刻，见上方诊断输出）
    cycle_end_times = [float(cyc.t[-1]) for cyc in solution.cycles]
    print(f"Cycle end times (absolute): {cycle_end_times}")

    # ── 均匀采样 ──
    t_eval = np.arange(0.0, solution.t[-1] + SAMPLE_DT, SAMPLE_DT)
    voltage = solution["Terminal voltage [V]"](t_eval)
    li_loss = solution["Loss of capacity to negative SEI [A.h]"](t_eval)
    sei_thick = solution[
        "X-averaged negative SEI thickness [m]"
    ](t_eval)

    # ── 导出 CSV ──
    csv_path = OUTPUT_DIR / "okane2022_dfn_sei_cc_2cycles.csv"
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("time_s,voltage_V,li_loss_ah,sei_thickness_m,cycle_end\n")
        for i, t in enumerate(t_eval):
            # 圈末采样窗口行：圈末时刻非 30s 整数倍，标记其所在的
            # 采样窗口（圈末与窗口起点差 ≤ 采样间隔，对寿命量的
            # 影响远小于对标准差）
            is_end = 1 if any(
                t <= te < t + SAMPLE_DT for te in cycle_end_times
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
