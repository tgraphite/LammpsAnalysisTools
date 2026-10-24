# level1_algo — lat vs GROMACS 逐值回归

`compare_gmx.py` 在运行时用 `gmx` 生成 Ground Truth，再用 lat 的 Python API（拓扑只加载一次）
逐对比较 RDF 与配位数。

```
python compare_gmx.py --system sys1 \
  --base /mnt/f/CRO/40_multi_electrolytes/11-systems \
  --rmax 2.0 --bin 0.002 --pairs LI:LI,LI:FSI,LI:FEC,LI:FEMC --frames 90:101:1
```

- 真实体系：sys1–sys6（178k–188k 原子，sys4/sys6 含 PMBA 聚合物）。
- 参数：0.002 nm 网格，rmax 2 nm，帧 9000–10000 ps（11 帧）。
- 指标：`max|Δg|`、`rms`、`mean`、峰位/峰高、`cn_max|d|`、`cn_rel`。

## 结果（2026-09-18）

| 体系 | max\|Δg\| | rms | cn_rel | 备注 |
|---|---|---|---|---|
| sys1 | 0.134 | 0.021 | 0.16% | 纯电解液 |
| sys2 | 0.137 | 0.022 | 0.21% | +TTE |
| sys3 | 0.155 | 0.020 | 0.16% | +MBA |
| sys4 | 0.151 | 0.022 | 0.16% | **+PMBA 聚合物** |
| sys5 | 0.151 | 0.022 | 0.20% | MBA+TTE |
| sys6 | 0.184 | 0.024 | 0.20% | **PMBA+TTE** |

首峰高度可达 g≈12.7，故 max\|Δg\|≤0.184 ≈ **首峰 <1.5%**；mean 偏置 ~1e-5。结论：RDF/CN 与 `gmx rdf`
在分子混合物 + 大分子体系上事实上一致。密度剖面 <1%、MSD（mol_com）曲线 ~3–7%。

详见 `docs/09_test_campaign.md`。

---

## LAMMPS 输入：`gmx_ref_lammps.py`（2026-10）

上面这条链路只覆盖 **GROMACS 原生输入**。LAMMPS `data`+`lammpstrj` 这一半 IO
（原子顺序、`Bonds` 段、`make_whole`）此前**没有任何原生 gmx 参考**，这是本节新增的部分：

```
python gmx_ref_lammps.py --data system.in.data --traj run.langevin.out.lammpstrj \
    --out gmx_ref --rdf-pairs C:O,C:C,O:O --json gmx_ref/metrics.json
```

它先用 `../convert/lammps_to_gmx.py` 把 LAMMPS 体系转成 gmx 能读的
`sys.top/sys.gro/sys.ndx/sys.tpr/whole.xtc`（按 **id 顺序**建拓扑，每个 LAMMPS 分子一个
`moleculetype`，σ 按 `r=0.5·σ·2^(1/6)` 反解成 lat 的 vdW 表），再对同一批帧跑：

| 模块 | gmx 命令 | 关键口径 |
|---|---|---|
| gyrate | `gmx gyrate`（System / 每分子） | gmx 在给了 `.tpr` 时会自己 `rmpbc`，等价 lat 的 make-whole |
| rms | `gmx rms -nomw -fit rot+trans` | gmx 默认 `-mw yes`。不加 `-nomw` 会差 ~1.9% |
| rmsf | `gmx rmsf -fit yes` | `-fit` 是布尔，不是 `rot+trans` |
| sasa | `gmx sasa -probe -ndots` | 半径按 `0.5·σ·2^(1/6)` 对齐 lat 的 vdW 表；lat 侧现为 `lat sasa`（最小镜像）|
| rdf | `gmx rdf -bin -rmax`（元素 index 组） | lat 用 `-intramol` 才与 gmx 默认同口径 |
| hbond | `gmx hbond -num` | lat 给体已限制为电负性元素；**受体集合不同**（lat 含卤素/S，gmx 只认 O），含 Cl⁻ 的体系会高 +13.8%（见 `docs/MANUAL_ALGORITHMS.md` §12.5）|
| freevol（无 gmx 对应）| — | 只报 `∫dV/dr = 1`、孔心数、connectivity 汇总，供人工核对 |

### 结果：4 链 PET 柱状体系（6200 原子、501 帧、盒 80³、`whole.xtc`）

| 量 | lat | gmx 2022.5 | mean\|Δ\| | max\|Δ\| | 相对 |
|---|---|---|---|---|---|
| Rg(all) | 26.8319 | 26.8316 | **0.0003 Å** | 0.0004 | 0.00% |
| Rg(每链 ×4) | 25.2776 / 25.3872 / 24.0370 / 25.4916 | 25.2778 / 25.3874 / 24.0371 / 25.4918 | **0.0002 Å** | 0.0005 | 0.00% |
| RMSD（fit all, ref frame 0） | 8.2829 | 8.2830 | **0.0001 Å** | 0.0002 | 0.00% |
| RMSF（逐原子） | 2.3167 | 2.3167 | **0.0007 Å** | 0.0038 | 0.03% |
| SASA（0.1.0）| 18 316.9 | 16 033.4 | 2 283.6 Å² | 2 917.3 | **14.24%**（已修，见下） |
| SASA（**0.2.0**）| **15 946.5** | 16 033.4 | **94.1 Å²** | 485.6 | **+0.59%** |
| rdf C–O | — | — | 0.090 | 6.92 | CN 0.26% |
| rdf C–C | — | — | 0.098 | 7.73 | CN 0.27% |
| rdf O–O | — | — | 0.044 | 2.30 | CN 0.27% |
| hbond（每帧计数） | 6.6028 | 6.6168 | **0.0140** | 1.0 | 0.21% |

判读（沿用 `compare_gmx.py` 的约定：**积分量优先**）：

- **gyrate / rms / rmsf 事实上逐位相同**（差在 `.gro` 的 0.001 nm 量化级别）。
  修正键序前，同一体系 Rg 高估 3.60×、RMSD 3.50×、RMSF 4.53×（`docs/07_decision_log.md` §6.1）。
- **hbond 501 帧里只有 1 帧差 1 条键**，均值差 0.014/帧。修正给体定义前是 41.7 vs 4.3/帧。
- **rdf 的 CN 差 0.26–0.27%**；逐点 `max|Δg|` 落在 r≈1.2 Å 的**成键峰**上（该处 g≈78，
  0.02 Å 网格上单 bin 计数很少），r>2 Å 的 mean\|Δg\| 只有 0.05，r>4 Å 只有 0.009。
  把 gmx 网格平移半个 bin 会**变差**（max 20.9），说明 lat 的 bin 中心约定是对的。
- **SASA 的 14.24% 差是真实缺陷，0.2.0 已修**：`_sasa` 的 `cKDTree` 没带 `boxsize`，
  完全没用 PBC（误差 99% 落在贴面 3 Å 内的原子上）。改成最小镜像的
  `lat.core.periodic.sasa_shrake_rupley` 后，501 帧均值 15 946.5 vs 16 033.4 Å²（**+0.59%**，
  逐帧 −0.17%…−0.85%）；残余来自球面采样点分布不同。详见 `docs/MANUAL_ALGORITHMS.md` §12.5
  与 `docs/11_examples_recheck.md`。
