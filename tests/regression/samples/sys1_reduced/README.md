# 回归样例 · sys1_reduced

自包含的锂电池电解液样例：800 个分子、6600 个原子、5 帧，同时提供 LAMMPS 与 GROMACS 两种表示。
用于单元/格式回归和 `lat select`/`rdf`/`msd` 等模块的冒烟。

## 文件

| 文件 | 说明 |
|---|---|
| `sys1_reduced.data` | LAMMPS data，`atom_style full`（id mol type q x y z）+ Masses + Bonds |
| `sys1_reduced.lammpstrj` | LAMMPS dump，列 `id mol type q x y z` |
| `sys1_reduced.gro` / `.xtc` | GROMACS 拓扑/轨迹（坐标级） |
| `sys1_reduced.mapping.csv` | 逐原子映射：id/mol/type/element/mass/charge/resname/atomname |
| `sys1_reduced.selections.yaml` | 元素级/分子级选择示例 |

## 体系

- 每个 resname（LI/FSI/FEC/FEMC）取前 200 个分子：`nmol=800`，`natoms=6600`。
- 轨迹帧：GROMACS 原始帧号 `[0, 20, 40, 60, 80]`（dump 的 TIMESTEP = 帧号）。
- 体系保持电中性（200 Li⁺ + 200 FSI⁻）；FEC/FEMC 中性。

## 类型映射

| type | element | count | LAMMPS 选择 | GROMACS 选择 |
|---|---|---|---|---|
| 1 | H | 1200 | `type 1` | `name H*` |
| 2 | Li | 200 | `type 2` | `resname LI` |
| 3 | C | 1400 | `type 3` | `name C*` |
| 4 | N | 200 | `type 4` | FSI 的 N |
| 5 | O | 2000 | `type 5` | |
| 6 | F | 1200 | `type 6` | |
| 7 | S | 400 | `type 7` | |

分子：LI 1–200，FSI 201–400，FEC 401–600，FEMC 601–800。

## 自检结果

- 坐标逐帧最大偏差 `9.5e-7 Å`（5 帧）；盒一致；`natoms=6600`。
- 电荷存在，净电荷 `2e-4 e`；键表 `6000`。
- `type→element` 与 GROMACS 原子名推断元素逐原子一致率 `1.0000`。
- Li 选择：LAMMPS `type 2` = 200，GROMACS `resname LI` = 200。

## 已知限制

- `.gro` 不携带质量/电荷/键；质量/电荷相关量的权威来源是 `.mapping.csv` 与 LAMMPS data。
- MDAnalysis 不按 `.lammpstrj` 扩展名自动选 reader，必须显式 `format="LAMMPSDUMP"`。
- 该样例无速度列，`lat transport` 不可用（需要 `vx vy vz`）。
