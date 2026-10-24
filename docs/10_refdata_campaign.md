# 10 · RefData 第三方订单兼容性与 gmx 回归战役

> 数据：`/mnt/e/RefData/{gromacs_real_projects,lammps_real_projects}`（复现小白用户订单场景）。
> GROMACS 侧 20 个订单、40 个 (topology,trajectory) 目录；LAMMPS 侧 22 个订单。
> 工具：`tests/compat/gmx_regress.py`（需 gmx）、`tests/compat/gmx_compat.py`（仅 lat）、
> `tests/compat/compat_run.py`（LAMMPS）。并行子代理分片执行。

## 1. 发现的缺陷与修复

| # | 现象 | 根因 | 处置 |
|---|---|---|---|
| 1 | `select` 崩溃：`unsupported format string passed to NoneType.__format__` | 粗粒化体系 element 为 None，`f"{el:<3}"` 崩 | 改 `str(el or '?')`；`select` 现显示 `element=?` |
| 2 | `hbond` 每帧 250 s、sys 时间 91 min（115k 原子） | 每个 D–H **暴力遍历全部受体** O(N_D·N_A) | 改用 `pairs_within` 周期性邻居搜索（DA/HA 两模式）；**250 s → ~5 s/帧** |
| 3 | 大盒系统 `probe` 直接 ERROR（网格 6.5e7 > 上限） | 网格点数硬上限 | 自动按 `(ntotal/cap)^{1/3}` 粗化 `grid_dr` 并 WARN，不再失败 |
| 4 | `-frame_range 0:3:1` 超出实际帧数时在 msd 里才报「fit needs ≥2 points」 | 已知帧数时未校验显式范围 | `resolve_frames` 在 `n_frames_known` 时校验上界，直接报清晰越界错误 |
| 5 | 截断的 xtc 报 MDAnalysis `OSError`/`AttributeError` 堆栈 | gmx reader 未兜底 | 包装为 `LatError`，提示文件截断/版本不匹配 |

## 2. 编译与调用约定（harness 修正，非 lat 缺陷）

- **`gmx rdf/msd` 的选择语法**：必须 `resname SOL`，裸 `SOL` 在部分 tpr 上触发
  gmx 内部错误 `listRanges does not have a first element with value 0`。
- **gmx 与 lat 的帧必须一致**：`gmx -dt` 与 lat `-frame_range:stride` 不等价；回归改用
  **同窗口 + stride=1**（`-b/-e`），density/gyrate/MSD 立即对齐。
- **tpr 版本高于 gmx 二进制**：RefData 中存在 VERSION **2025.1** 的 tpr，本机 gmx 2023.2
  读拓扑失败（选择为空 → 上述 gmx 内部错误），而 **lat/MDAnalysis 读取正常**。这类系统
  gmx 侧 ground truth 不可得，harness 记 `gmx_unavailable` 并跳过对比。
- **原子数配对**：绝不能把 tpr 与原子数不同的轨迹配对（如 `ion.tpr` vs `md.trr`）。

## 3. GROMACS 原生回归结果（窗口 = 末尾 20 帧，stride 1）

代表性系统（其余同量级）：

| 系统 | 原子 | rdf max rel | cn max rel | density max\|d\| (kg/m³) | Rg max\|d\| (nm) | MSD max\|d\| (nm²) |
|---|---|---|---|---|---|---|
| spdf（聚合物电解质） | 65021 | 0.04–0.05 | 0.0015 | ~16 | 3e-4 | 5e-5 |
| 25102824446/11 equ | 17009 | 0.007 | 0.0015 | ~17 | 2e-5 | 5e-3 |
| 25120554354/equ | 115728 | 0.012 | 0.0015 | ~46 | 9e-4 | 0.25 |
| 25112735551（小分子） | 6369 | 0.07 | 0.0016 | 尖峰见下 | 5e-6 | 0.85 |

- **RDF/CN**：CN 相对误差普遍 **≤0.2%**；g(r) 相对偏差大多 <5%。个别达到 0.2–0.4 的系统是
  **单帧最小化快照**上的尖锐首峰，或粗粒化/小体系，配合 CN 仍一致 → 属半 bin 网格/采样效应。
- **density**：均值偏置 ~0.1% 以内；max|d| 出现在**极尖锐界面/峰**（如 25112735551 剖面
  64→4037 kg/m³，11 帧），属 bin 边缘效应。
- **gyrate**：修好「跨 PBC 分子需 make-whole」（见下）后，最大偏差 **≤1e-3 nm**。
- **MSD**：正常扩散体系一致到 **<0.5%**；长步长/极短窗口（11 帧、步长 10 ns、位移 ~5 nm）
  有 3–5% 的有限采样差；gmx 在个别体系数值下溢（~1e-41）为其自身问题。

## 4. 由此修复的 PBC 物理缺陷（重要）

`gyrate` 原本用 wrapped 坐标，对跨盒的聚合物 Rg **系统性偏低**（nafion 3.7 vs gmx 6.8 nm）。
新增 `core/unwrap.py`：按**键连接生成树**把分子复原（gmx `rmpbc` 语义），而非「相对首原子 MIC」
（对长于半盒的链会失败）。`gyrate` 与 `rms` 均改用之。修复后 nafion Rg 偏差 0.00015 nm。

## 5. LAMMPS 侧

`tests/compat/compat_run.py` 覆盖 8 个有效 (data,lammpstrj) 对：**0 失败**。此前
`probe`（0.2.0 起为 `freevol`）网格过大错误已由自动粗化修复。其余非通过项均为 SKIP
（无速度列 → transport；无 Bonds → hbond/agg/cluster/residence；帧数不足 → rms/msd）。

## 6. 已知限制（未解决，记录在案）

- `hbond` 的 `_acf`：0.2.0 已改为增量式并丢弃已结束的氢键对（120k 对 × 601 帧从分钟级到
  0.04 s），`present` 矩阵仍为 (键数 × T)，超长轨迹建议按窗口/子采样运行。
- 粗粒化体系（无元素、无化学键）不支持 `freevol`/`sasa`/`hbond`（报错清晰）。
- 反应轨迹逐帧原子数变化 → lat 报错并建议限定窗口（与 09 一致）。
- gmx 版本低于 tpr 版本时无法取得 gmx GT；lat 仍可分析。
