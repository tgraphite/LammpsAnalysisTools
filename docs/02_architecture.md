# 02 · 系统架构

## 1. 分层

```
┌─────────────────────────────────────────────────────────────┐
│ CLI 层        lat <module> [-s data] [-f traj] [选项...]      │
│               统一 -s/-f 别名、分发、-h 自解释、错误分级        │
├─────────────────────────────────────────────────────────────┤
│ 模块层        modules/*: 纯函数 module(ctx, params) -> [Table]│
│               只依赖 core 抽象；不直接碰 IO 细节               │
├─────────────────────────────────────────────────────────────┤
│ 核心层 core/                                                  │
│   io/          lammps-data / lammpstrj 读；MDAnalysis 读 gmx(内)│
│   topology/    类型/电荷/质量/键/分子ID（用户提供或 data 自带）  │
│   selection/   静态 MDAnalysis 语义 + ndx；LAMMPS type 映射     │
│   pbc.py       正交/三斜 minimum image                          │
│   neighbor.py  cell list / KDTree                              │
│   grid.py      自动网格规则、体素累积、smooth                   │
│   correlate.py 自相关、块平均、线性拟合、误差                   │
│   output.py    唯一写文件处：Table 契约 + 序列化 + 原子写        │
│   errors.py / logger.py / context.py                           │
├─────────────────────────────────────────────────────────────┤
│ 绘图层        plotting/: plot1d / plot2d，只读 CSV              │
└─────────────────────────────────────────────────────────────┘
```

## 2. 数据流

```
输入文件 ──(io)──> Context(拓扑 + 逐帧坐标访问器 + 单位表 + 帧数)
                      │
用户选择/ndx ─(selection)─> 原子组
                      │
module(ctx, params) ──(neighbor/grid/correlate)──> [Table]
                      │
                   (output) ──> 原子写 CSV
                      │
                   (plotting) ──> 图
```

- **Context 不持有整条轨迹**：只提供 `for frame in ctx.frames(range, stride)` 的迭代器（流式，M9）。
- **模块是纯函数**：相同 Context+params → 相同 Table（除显式 seed）。

## 3. 模块接口契约

```python
@register("rdf")
class RDF(Module):
    # 从 config/lat.yaml 注入默认值；代码内 default 只能来自 config
    params_schema = {...}

    def run(self, ctx, *, sel1, sel2, nbins, cutoff, frame_range, ...) -> list[Table]:
        ...
```

- `Table`：`title / x(field,unit) / y(field,unit)[…] / scalars / warnings / meta`。
- 模块**不得**自己 open 文件；统一由 `output.py` 落盘。

## 4. 三条硬边界

1. `modules/*` 只依赖 `core` 抽象——将来替换 IO 后端不动算法。
2. `output.py` 是唯一写文件的地方。
3. `config/lat.yaml` 是唯一默认值来源；模块代码出现任何 default 数值即视为 bug。

## 5. 配置加载与优先级

```
CLI 参数  >  config/lat.yaml  >  FATAL
```

- YAML 含：`meta.units`、`grid` 默认、`output` 契约、`logger`、`modules.<name>.defaults`、`modules.<name>.help`（指向 `docs/MANUAL_ALGORITHMS.md`）。
- 找不到键 / 类型不符 / 单位缺失 → 立即 FATAL（无 backup）。

## 6. 错误分级

| 级别 | 行为 |
|---|---|
| FATAL | 打印原因，退出码非 0，不产生输出文件 |
| WARN | 继续，写入 header `warning:` 与屏显 |
| INFO | 屏显（`-v` 才出） |
| DEBUG | 审计用（seed/sha256/内存），仅开发者可见，成品不暴露 |

## 7. 目录职责

- `config/`：唯一配置。
- `docs/`：来龙去脉、元要求、契约、算法、模块文档（`-h` 同源）。
- `src/lat/`：实现。
- `tests/regression/`：两级回归；`tests/data/`：小型自造体系。

## 8. 路线图

| 阶段 | 交付 | 验收 |
|---|---|---|
| P0 | core：io/topology/selection/pbc/grid/correlate/output + `select` + `plot1d/2d` | 读 LAMMPS data/lammpstrj + gmx xtc；`select` 正确；输出契约成型 |
| P1 | `rdf` `density` `msd` `rms` | 对 40/11-systems 的 gmx xvg 数值回归 |
| P2 | `hbond` `agg` `cluster` `pistack` | hbond 自造体系点数；agg 用 Li-FSI；pistack 用现有测试 |
| P3 | `gyrate` `contact` `order` `probe` `residence` | 各自小体系 |
| P4 | `transport` `occupancy` | 速度列/区域事件 |

## 9. 测试契约

- **golden**：固定输入 + 输出 CSV 数值容差比对（回归见 05_regression.md）。
- **unit**：PBC/三斜距离、网格自动规则、归一化自检、自相关。
- **CLI**：`-s/-f` 别名、`-h`、FATAL 退出码、原子写幂等。
- **perf**：183k 原子体系峰值内存与耗时基准。
