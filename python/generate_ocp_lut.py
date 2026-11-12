"""生成/更新 params/*.yml 的 ocp: 段（内置 OCP 公式的等间隔采样点）。

OCP 数据以 stoich/voltage 采样序列存入 YAML，模型加载后线性插值使用；
未提供 ocp: 段的参数文件回退到内置 Chen2020 公式（见 core/src/physics.cpp）。

本脚本幂等：重复运行会先移除旧的 ocp 段再重新追加。
"""
import math
from pathlib import Path

# ── OCP 公式（与 core/src/physics.cpp 严格一致）──


def graphite_ocp(x: float) -> float:
    """石墨负极 OCP — Safari & Delacourt (2009) 简化拟合"""
    x = min(max(x, 0.001), 0.999)
    return (0.194 + 1.5 * math.exp(-120.0 * x)
            + 0.0351 * math.tanh((x - 0.286) / 0.083)
            - 0.0045 * math.tanh((x - 0.9) / 0.119)
            - 0.035 * math.tanh((x - 0.99) / 0.05)
            - 0.0147 * math.tanh((x - 0.5) / 0.034)
            - 0.102 * math.tanh((x - 0.194) / 0.142)
            - 0.022 * math.tanh((x - 0.9) / 0.0164)
            - 0.011 * math.tanh((x - 0.124) / 0.0226)
            + 0.0155 * math.tanh((x - 0.105) / 0.029))


def nmc111_ocp(y: float) -> float:
    """NMC111 正极 OCP — Lu et al. (2013) 拟合"""
    y = min(max(y, 0.001), 0.999)
    return (-0.8090 * y + 4.4875
            - 0.0428 * math.tanh(18.5138 * (y - 0.5542))
            - 17.7326 * math.tanh(15.7890 * (y - 0.3117))
            + 17.5842 * math.tanh(15.9308 * (y - 0.3120)))


def graphite_chen2020_ocp(x: float) -> float:
    """石墨负极 OCP — Chen et al. (2020), PyBaMM graphite_LGM50_ocp_Chen2020"""
    x = min(max(x, 1e-6), 1.0 - 1e-6)
    return (1.9793 * math.exp(-39.3631 * x) + 0.2482
            - 0.0909 * math.tanh(29.8538 * (x - 0.1234))
            - 0.04478 * math.tanh(14.9159 * (x - 0.2769))
            - 0.0205 * math.tanh(30.4444 * (x - 0.6103)))


def nmc_chen2020_ocp(y: float) -> float:
    """NMC 正极 OCP — Chen et al. (2020), PyBaMM nmc_LGM50_ocp_Chen2020"""
    y = min(max(y, 1e-6), 1.0 - 1e-6)
    return (-0.8090 * y + 4.4875
            - 0.0428 * math.tanh(18.5138 * (y - 0.5542))
            - 17.7326 * math.tanh(15.7890 * (y - 0.3117))
            + 17.5842 * math.tanh(15.9308 * (y - 0.3120)))


N = 201  # [0,1] 上 201 点（步长 0.005，线性插值误差 < 1 mV）

OCV_MARKER = "# ── 开路电位 (OCP) 曲线"


def fmt_flow_list(values: list[float]) -> str:
    """格式化为折行的 YAML flow 序列（[...] 包裹，每行最多 10 个数）"""
    lines = []
    for i in range(0, len(values), 10):
        chunk = ", ".join(f"{v:.6g}" for v in values[i:i + 10])
        prefix = "      [" if i == 0 else "      "
        sep = "," if i + 10 < len(values) else ""
        lines.append(prefix + chunk + sep)
    lines[-1] += "]"
    return "\n".join(lines)


def ocp_block(neg_fn, pos_fn, neg_label: str, pos_label: str) -> str:
    stoich = [i / (N - 1) for i in range(N)]
    return (
        f"{OCV_MARKER}：stoich/voltage 采样点，线性插值（端点外钳位）──\n"
        "ocp:\n"
        "  negative_electrode:\n"
        f"    # {neg_label}\n"
        "    stoich:\n"
        f"{fmt_flow_list(stoich)}\n"
        "    voltage:\n"
        f"{fmt_flow_list([neg_fn(s) for s in stoich])}\n"
        "  positive_electrode:\n"
        f"    # {pos_label}\n"
        "    stoich:\n"
        f"{fmt_flow_list(stoich)}\n"
        "    voltage:\n"
        f"{fmt_flow_list([pos_fn(s) for s in stoich])}\n"
    )


# 文件 → (负极 OCP, 正极 OCP, 负极标注, 正极标注)
TARGETS = {
    "chen2020.yml": (graphite_chen2020_ocp, nmc_chen2020_ocp,
                     "graphite_LGM50 (Chen2020)", "NMC_LGM50 (Chen2020)"),
    "dfn_chen2020.yml": (graphite_chen2020_ocp, nmc_chen2020_ocp,
                         "graphite_LGM50 (Chen2020)", "NMC_LGM50 (Chen2020)"),
    "spme_chen2020.yml": (graphite_chen2020_ocp, nmc_chen2020_ocp,
                          "graphite_LGM50 (Chen2020)", "NMC_LGM50 (Chen2020)"),
    "graphite_nmc111.yml": (graphite_ocp, nmc111_ocp,
                            "Graphite (Safari 2009)", "NMC111 (Lu 2013)"),
}


def main() -> None:
    params_dir = Path(__file__).resolve().parents[1] / "params"
    for name, (neg_fn, pos_fn, neg_label, pos_label) in TARGETS.items():
        path = params_dir / name
        text = path.read_text(encoding="utf-8")
        # 幂等：截掉旧的 ocp 段（从注释标记到文件尾）
        idx = text.find(OCV_MARKER)
        if idx >= 0:
            text = text[:idx].rstrip() + "\n"
        if not text.endswith("\n"):
            text += "\n"
        text += "\n" + ocp_block(neg_fn, pos_fn, neg_label, pos_label)
        path.write_text(text, encoding="utf-8")
        print(f"updated {path}")


if __name__ == "__main__":
    main()
