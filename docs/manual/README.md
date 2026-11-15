# lat 模块手册

lat（LAMMPS Analysis Tool）把分子动力学后处理中高频、可复用的分析收敛为边界清晰的子命令。本目录是每个模块的 help manual；先读 [common.md](common.md) 了解公共参数、输入格式、选择语法与输出契约。

## 快速开始

    pip install lammps-analysis-tools
    lat --help
    lat rdf -s system.data -f run.lammpstrj -sel1 "type 2" -sel2 "type 5" -cutoff 10 -o out

## 模块索引

| 模块 | 用途 | 手册 |
|---|---|---|
| rdf | 径向分布函数与配位数 | [rdf.md](rdf.md) |
| density | 一维/二维密度剖面 | [density.md](density.md) |
| msd | 均方位移与扩散系数 | [msd.md](msd.md) |
| gyrate | 回转半径 Rg | [gyrate.md](gyrate.md) |
| rms | RMSD 与 RMSF | [rms.md](rms.md) |
| contact | 接触数与距离 | [contact.md](contact.md) |
| hbond | 氢键数量、寿命与几何分布 | [hbond.md](hbond.md) |
| agg | SSIP/CIP/AGG 离子聚合形态 | [agg.md](agg.md) |
| cluster | 配位数与团簇尺寸分布 | [cluster.md](cluster.md) |
| pistack | π–π 堆积 | [pistack.md](pistack.md) |
| freevol | 自由体积、孔径与连通性 | [freevol.md](freevol.md) |
| sasa | 溶剂可及表面积 | [sasa.md](sasa.md) |
| residence | 配体滞留时间 | [residence.md](residence.md) |
| transport | 扩散、Nernst–Einstein 与 Green–Kubo | [transport.md](transport.md) |
| occupancy | 区域占据、穿越与首次通过 | [occupancy.md](occupancy.md) |
| select | 选择检查与 .ndx 导出 | [select.md](select.md) |
| describe | 模块清单与契约 | [describe.md](describe.md) |
| plot1d | 一维 CSV 绘图 | [plot1d.md](plot1d.md) |
| plot2d | 二维 CSV 绘图 | [plot2d.md](plot2d.md) |

## 约定

- 所有模块共享 [公共参数](common.md#公共参数) 与 [输出契约](common.md#输出约定)。
- 默认值来自配置文件；每个模块手册的参数表由代码元数据生成，并由测试静态校验。
- 结果写入 CSV，便于脚本批处理与复现；绘图单独用 plot1d/plot2d。
- 命令示例中的 BASE 指 -s <topology> -f <trajectory>。

## 版本

当前版本 1.0.0。变更记录见 [../CHANGELOG.md](../CHANGELOG.md)。
