"""生成 PyBaMM DFN 耦合热模型（lumped thermal）基准数据。

使用 PyBaMM 的 DFN 模型 + options {"thermal": "lumped"}（Timms2021 集总
热模型：单节点温度 ODE，产热三项 = 欧姆 + 界面 + 熵热，Newton 冷却），
在 1C 恒流放电工况下运行仿真，将时间/电压/温度数据导出为 CSV，
作为 FluxCell DFN 耦合热模型的验证基准。

模型与参数：
    - pybamm.lithium_ion.DFN(options={"thermal": "lumped"})
    - Chen2020 参数集（含热参数默认值：h_total=10 W/(m²·K)，
      A_cooling=0.00531 m²，V_cell=2.42e-5 m³，初始/环境温度 298.15 K）
    - 熵变系数 dU/dT：Chen2020 未提供（默认 0），此处按 params/ 下的
      用户测量数据（neg_graphite_dudt.csv / pos_ncm_dudt.csv，两列：
      嵌锂量 sto、dU/dT [mV/K]）以分段线性插值设置，与 FluxCell 的
      EntropyLut 逐点一致

输出文件（validation/reference_data/）：
    chen2020_dfn_thermal_entropic_discharge_1C.csv
        列: time_s, voltage_V, temperature_K
"""

from pathlib import Path

import numpy as np
import pybamm

# 仓库根（python/validation/ 上两级）与熵变系数 CSV（与 FluxCell 共用）
REPO_ROOT = Path(__file__).resolve().parents[2]
ENTROPY_CSV = {
    "Negative electrode OCP entropic change [V.K-1]":
        REPO_ROOT / "params" / "neg_graphite_dudt.csv",
    "Positive electrode OCP entropic change [V.K-1]":
        REPO_ROOT / "params" / "pos_ncm_dudt.csv",
}


def load_dudt_csv(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """加载熵变系数 CSV（嵌锂量 sto、dU/dT [mV/K]，换算为 V/K）。"""
    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    data = data[~np.isnan(data).any(axis=1)]  # 过滤空行/表头残留
    return data[:, 0], data[:, 1] * 1e-3


def _install_entropy_changes(params: pybamm.ParameterValues) -> None:
    """把用户测量 dU/dT 以 linear Interpolant 注入参数集（默认参数键是
    标量 0.0；替换为函数后 PyBaMM 会以各电极化学计量比 sto 求值）。
    节点扩展至 [0, 1] 端点常数延拓：初始化 Newton 会探索 sto 全域，
    且与 C++ EntropyLut 的 clamp 口径一致。"""
    for key, path in ENTROPY_CSV.items():
        sto, dudt = load_dudt_csv(path)
        sto = np.concatenate(([0.0], sto, [1.0]))
        dudt = np.concatenate(([dudt[0]], dudt, [dudt[-1]]))
        params[key] = (
            lambda s, _x=sto, _y=dudt: pybamm.Interpolant(
                _x, _y, s, interpolator="linear"
            )
        )

# 输出目录
OUTPUT_DIR = Path(__file__).parent / "reference_data"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# 待验证的 C-rate 列表
C_RATES = [1.0]


def run_thermal_discharge(c_rate: float, parameter_set: str = "Chen2020") -> Path:
    """运行 DFN + lumped thermal 恒流放电仿真并导出 CSV。

    Args:
        c_rate: 放电倍率（如 1.0 表示 1C）
        parameter_set: PyBaMM 参数集名称

    Returns:
        导出的 CSV 文件路径
    """
    # ── 模型与参数（lumped 热模型与 Chen2020 热参数）──
    model = pybamm.lithium_ion.DFN(
        name="DFN lumped thermal reference",
        options={"thermal": "lumped"},
    )
    params = pybamm.ParameterValues(parameter_set)
    _install_entropy_changes(params)  # 熵变系数 dU/dT（用户测量 LUT）

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
    temperature_k = solution["Volume-averaged cell temperature [K]"].entries

    # ── 导出 CSV ──
    csv_path = OUTPUT_DIR / f"chen2020_dfn_thermal_entropic_discharge_{c_rate:g}C.csv"
    with open(csv_path, "w", encoding="utf-8") as f:
        f.write("time_s,voltage_V,temperature_K\n")
        for t, v, temp in zip(time_s, voltage_v, temperature_k):
            f.write(f"{t:.6e},{v:.9e},{temp:.9e}\n")

    print(
        f"[OK] {c_rate:g}C 耦合热放电: {len(time_s)} 个时间点, "
        f"终止电压 {voltage_v[-1]:.4f} V, 终止温度 {temperature_k[-1]:.2f} K, "
        f"耗时 {time_s[-1]:.1f} s"
    )
    return csv_path


def main() -> None:
    print(f"PyBaMM version: {pybamm.__version__}")
    print(f"Output directory: {OUTPUT_DIR}\n")

    for c_rate in C_RATES:
        run_thermal_discharge(c_rate)

    print("\nThermal reference data generated successfully.")


if __name__ == "__main__":
    main()
