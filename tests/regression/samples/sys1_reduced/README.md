# 回归样例 · sys1_reduced

由 `40_multi_electrolytes/11-systems/sys1_npt.{tpr,xtc}` 缩减并转换而来，用于 **两级输入回归**：
同一体系的 GROMACS 与 LAMMPS 表示，应给出相同的 lat 分析结果。

## 来源与生成

```
python ../../convert/convert_gmx_to_lammps.py \
  -s /mnt/f/CRO/40_multi_electrolytes/11-systems/sys1_npt.tpr \
  -f /mnt/f/CRO/40_multi_electrolytes/11-systems/sys1_npt.xtc \
  -o sys1_reduced --per-resname 200 --frames 5
python ../../convert/make_selections.py -p sys1_reduced
python ../../convert/verify_lammps_roundtrip.py -p sys1_reduced
```

- 每个 resname（LI/FSI/FEC/FEMC）取前 200 个分子：`nmol=800`，`natoms=6600`。
- 轨迹帧：GROMACS 原始帧号 `[0, 20, 40, 60, 80]`（dump 的 TIMESTEP = 帧号）。
- 体系保持电中性（200 Li⁺ + 200 FSI⁻）；FEC/FEMC 中性。

## 文件

| 文件 | 说明 |
|---|---|
| `sys1_reduced.data` | LAMMPS data，`atom_style full`（id mol type q x y z）+ Masses + Bonds |
| `sys1_reduced.lammpstrj` | LAMMPS dump，列 `id mol type q x y z` |
| `sys1_reduced.gro` / `.xtc` | 缩减版 GROMACS 拓扑/轨迹（坐标级） |
| `sys1_reduced.mapping.csv` | **权威逐原子映射**：id/mol/type/element/mass/charge/resname/resname/atomname |
| `sys1_reduced.selections.yaml` | 元素级/分子级选择 + 各 P0 模块候选选择 |
| `sys1_reduced.verify.json` | 回归自检结果（见下） |

## 类型映射（元素顺序按原子序数）

| type | element | count | LAMMPS 选择 | GROMACS 选择 |
|---|---|---|---|---|
| 1 | H | 1200 | `type 1` | `name H*` |
| 2 | Li | 200 | `type 2` | `resname LI` |
| 3 | C | 1400 | `type 3` | `name C*` |
| 4 | N | 200 | `type 4` | （FSI 的 N） |
| 5 | O | 2000 | `type 5` | |
| 6 | F | 1200 | `type 6` | |
| 7 | S | 400 | `type 7` | |

分子（LAMMPS 用 mol 区间）：`LI` 1–200，`FSI` 201–400，`FEC` 401–600，`FEMC` 601–800。

## 自检结果（全部通过）

- 坐标逐帧最大偏差 `9.5e-7 Å`（5 帧）；盒一致；`natoms=6600`。
- 电荷存在，净电荷 `2e-4 e`；键表 `6000`。
- `type→element` 与 GROMACS 原子名推断元素逐原子一致率 `1.0000`；元素计数完全一致。
- `Li` 选择：LAMMPS `type 2` = 200，GROMACS `resname LI` = 200。

## 已知限制

- `.gro` 不携带质量/电荷/键；MDAnalysis 读 `.gro` 时会把 `LI` 的质量猜成 `0`。
  因此 GROMACS 样例仅用于**坐标级**回归；质量/电荷相关量的权威来源是 `.mapping.csv`
  与 LAMMPS data。要跑含质量的 GROMACS 级对比，请用原始 `.tpr`（见 `docs/05_regression.md`）。
- MDAnalysis 不按 `.lammpstrj` 扩展名自动选 reader，必须显式 `format="LAMMPSDUMP"`。
- 该样例无速度列，`lat transport` 不可用（需要 `vx vy vz`）。
