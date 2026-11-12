# FluxCell Python 验证环境（PyBaMM 基准）

使用 [uv](https://docs.astral.sh/uv/) 管理的 Python 环境，用于生成 PyBaMM DFN 基准数据并与 FluxCell 结果对比。

## 初始化

```powershell
# 在 python/ 目录下
uv sync          # 创建虚拟环境并安装依赖（首次运行）
```

## 目录结构

```
python/
├── pyproject.toml              # uv 项目配置（PyBaMM 等依赖）
├── .python-version              # Python 版本锁定
└── validation/
    ├── generate_reference.py   # 生成 PyBaMM DFN 基准数据 → reference_data/
    ├── compare_results.py       # 对比 FluxCell CSV 结果与基准数据
    └── reference_data/          # 基准数据输出（CSV，git 可提交）
```

## 使用方法

```powershell
# 1. 生成 PyBaMM 基准数据（多种 C-rate 放电曲线）
uv run python validation/generate_reference.py

# 2. FluxCell 仿真完成后，对比结果
uv run python validation/compare_results.py --fluxcell <path/to/fluxcell.csv> --crate 1C
```

## 验证标准

- 电压曲线最大绝对误差 < 10 mV
- 电压曲线平均绝对误差 < 2 mV
- 容量误差 < 1%
