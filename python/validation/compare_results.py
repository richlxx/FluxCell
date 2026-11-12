"""对比 FluxCell 仿真结果与 PyBaMM 基准数据。

读取 FluxCell 导出的 CSV（time,voltage,current,temperature 格式），
或直接解析 gtest 诊断快照日志（如 *DFNDischarge1CDiagnostic* 输出的
"t=  ... s  V=..." 行），与 PyBaMM DFN 基准曲线对比，计算误差指标并绘制对比图。

用法：
    uv run python validation/compare_results.py --fluxcell result.csv --crate 1C
    fluxcell_tests.exe --gtest_filter='*DFNDischarge1CDiagnostic*' | \
        python validation/compare_results.py --gtest-log - --model dfn --crate 1C

误差指标：
    - 最大绝对误差 (max |ΔV|)
    - 平均绝对误差 (mean |ΔV|)
    - RMS 误差
"""

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np

REFERENCE_DIR = Path(__file__).parent / "reference_data"


def parse_gtest_snapshots(lines: list[str]) -> dict[str, np.ndarray]:
    """解析 gtest 诊断快照行（"t=  ... s  V=..." 格式）为列字典。"""
    pat = re.compile(r"t=\s*(\d+)\s+s\s+V=([\d.]+)")
    times: list[float] = []
    volts: list[float] = []
    for line in lines:
        m = pat.search(line)
        if m:
            times.append(float(m.group(1)))
            volts.append(float(m.group(2)))
    if not times:
        raise ValueError("gtest 日志中未解析到快照行（期望格式: t=  ... s  V=...）")
    return {"time": np.array(times), "voltage": np.array(volts)}


def load_csv_columns(path: Path) -> dict[str, np.ndarray]:
    """加载 CSV 文件为列字典。"""
    with open(path, encoding="utf-8") as f:
        reader = csv.reader(f)
        header = next(reader)
        columns = {name: [] for name in header}
        for row in reader:
            for name, value in zip(header, row):
                columns[name].append(float(value))
    return {name: np.array(values) for name, values in columns.items()}


def interpolate_reference(ref_time: np.ndarray, ref_voltage: np.ndarray,
                           flux_time: np.ndarray) -> np.ndarray:
    """将基准电压曲线插值到 FluxCell 的时间点上。"""
    return np.interp(flux_time, ref_time, ref_voltage)


def compute_metrics(flux_voltage: np.ndarray, ref_voltage: np.ndarray) -> dict[str, float]:
    """计算误差指标（单位：mV）。"""
    diff_v = flux_voltage - ref_voltage
    diff_mv = diff_v * 1000.0
    return {
        "max_abs_err_mv": float(np.max(np.abs(diff_mv))),
        "mean_abs_err_mv": float(np.mean(np.abs(diff_mv))),
        "rms_err_mv": float(np.sqrt(np.mean(diff_mv**2))),
        "n_points": len(diff_v),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Compare FluxCell results with PyBaMM reference")
    parser.add_argument("--fluxcell", default=None,
                        help="FluxCell 导出的 CSV 文件路径")
    parser.add_argument("--gtest-log", default=None,
                        help="gtest 诊断日志（快照行），'-' 表示读 stdin；与 --fluxcell 二选一")
    parser.add_argument("--crate", default="1C", help="基准数据对应的 C-rate（如 1C）")
    parser.add_argument("--model", default="dfn", choices=["dfn", "spm"],
                        help="基准模型类型（dfn 或 spm，决定基准文件名前缀）")
    parser.add_argument("--tol-max-mv", type=float, default=10.0, help="最大误差容限 [mV]")
    parser.add_argument("--tol-mean-mv", type=float, default=2.0, help="平均误差容限 [mV]")
    parser.add_argument("--plot", action="store_true", help="绘制对比图（需要 matplotlib）")
    args = parser.parse_args()
    if bool(args.gtest_log) == bool(args.fluxcell):
        parser.error("--gtest-log 与 --fluxcell 必须二选一")

    # ── 加载 FluxCell 结果 ──
    if args.gtest_log:
        if args.gtest_log == "-":
            lines = sys.stdin.read().splitlines()
        else:
            lines = Path(args.gtest_log).read_text(encoding="utf-8", errors="replace").splitlines()
        flux_data = parse_gtest_snapshots(lines)
        flux_label = "<gtest log>"
    else:
        flux_path = Path(args.fluxcell)
        if not flux_path.exists():
            raise FileNotFoundError(f"FluxCell result file not found: {flux_path}")
        flux_data = load_csv_columns(flux_path)
        flux_label = str(flux_path)

    # ── 加载 PyBaMM 基准 ──
    crate_str = args.crate.removesuffix("C")
    ref_path = REFERENCE_DIR / f"chen2020_{args.model}_discharge_{float(crate_str):g}C.csv"
    if not ref_path.exists():
        raise FileNotFoundError(
            f"Reference data not found: {ref_path}\n"
            "请先运行: uv run python validation/generate_reference.py"
        )
    ref_data = load_csv_columns(ref_path)

    # ── 截取共同时间范围并插值 ──
    t_common_end = min(flux_data["time"].max(), ref_data["time_s"].max())
    mask = flux_data["time"] <= t_common_end
    flux_time = flux_data["time"][mask]
    flux_voltage = flux_data["voltage"][mask]

    ref_voltage_interp = interpolate_reference(
        ref_data["time_s"], ref_data["voltage_V"], flux_time
    )

    # ── 计算误差 ──
    metrics = compute_metrics(flux_voltage, ref_voltage_interp)

    print("=" * 60)
    print(f"FluxCell 结果: {flux_label}")
    print(f"PyBaMM 基准:   {ref_path}")
    print(f"对比时间范围: 0 ~ {t_common_end:.1f} s")
    print("-" * 60)
    print(f"数据点数:         {metrics['n_points']}")
    print(f"最大绝对误差:     {metrics['max_abs_err_mv']:.3f} mV")
    print(f"平均绝对误差:     {metrics['mean_abs_err_mv']:.3f} mV")
    print(f"RMS 误差:         {metrics['rms_err_mv']:.3f} mV")
    print("-" * 60)

    # ── 判定 ──
    passed = (
        metrics["max_abs_err_mv"] < args.tol_max_mv
        and metrics["mean_abs_err_mv"] < args.tol_mean_mv
    )
    verdict = "PASS ✓" if passed else "FAIL ✗"
    print(f"验证结果: {verdict} "
          f"(容限: max < {args.tol_max_mv} mV, mean < {args.tol_mean_mv} mV)")
    print("=" * 60)

    # ── 绘图（可选）──
    if args.plot:
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        ax1.plot(ref_data["time_s"], ref_data["voltage_V"], "b-", label="PyBaMM (reference)")
        ax1.plot(flux_time, flux_voltage, "r--", label="FluxCell")
        ax1.set_ylabel("Voltage [V]")
        ax1.set_title(f"FluxCell vs PyBaMM {args.model.upper()} — {args.crate} discharge")
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        diff_mv = (flux_voltage - ref_voltage_interp) * 1000.0
        ax2.plot(flux_time, diff_mv, "g-")
        ax2.axhline(0, color="k", linewidth=0.5)
        ax2.axhline(args.tol_max_mv, color="r", linestyle="--", linewidth=0.8, label=f"+{args.tol_max_mv} mV")
        ax2.axhline(-args.tol_max_mv, color="r", linestyle="--", linewidth=0.8, label=f"-{args.tol_max_mv} mV")
        ax2.set_xlabel("Time [s]")
        ax2.set_ylabel("Voltage error [mV]")
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        out_png = flux_path.with_suffix(".comparison.png") if args.fluxcell else Path("dfn_comparison.png")
        fig.tight_layout()
        fig.savefig(out_png, dpi=150)
        print(f"\nComparison plot saved to: {out_png}")

    raise SystemExit(0 if passed else 1)


if __name__ == "__main__":
    main()
