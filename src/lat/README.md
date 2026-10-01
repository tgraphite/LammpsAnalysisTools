# src/lat — 实现（P0 进行中）

> 规则：模块只依赖 `core` 抽象；唯一写文件处是 `core/output.py`；默认值只来自 `config/lat.yaml`。
>
> 已实现：`core`（errors/logger/config/table/pbc/topology/frame/neighbor/grid/correlate/
> selection/context/output + io lammps-data/lammps-dump/gmx）、`modules`（`select`/`describe`/`rdf`）、
> `plotting`（`plot1d`/`plot2d`）、`cli`。/ 待实现：其余 P0 模块。

```
lat/
├── cli.py            # 统一 -s/-f、子命令分发、-h、错误分级
├── config.py         # 加载 config/lat.yaml（无 backup，缺键 FATAL）
├── core/
│   ├── io/           # lammps data/lammpstrj 读；gmx 内部读；CSV 写
│   ├── topology.py   # type/charge/mass/bond/mol（用户提供或 data 自带）
│   ├── selection.py  # 静态 MDAnalysis 语义 + ndx；LAMMPS type 映射
│   ├── pbc.py        # 正交/三斜 minimum image
│   ├── neighbor.py   # cell list / KDTree
│   ├── grid.py       # 自动网格 + 体素累积 + smooth
│   ├── correlate.py  # 自相关/块平均/拟合
│   ├── output.py     # Table 契约 + 序列化 + 原子写（唯一写盘处）
│   ├── errors.py     # FATAL/WARN/INFO
│   └── context.py    # Context + 流式帧迭代器
├── modules/          # 每个子命令一个纯函数 module(ctx, params) -> [Table]
│   ├── rdf.py density.py hbond.py rms.py msd.py agg.py cluster.py
│   ├── pistack.py probe.py residence.py transport.py occupancy.py
│   └── gyrate.py contact.py order.py select.py
└── plotting/
    └── plot1d.py plot2d.py
```
