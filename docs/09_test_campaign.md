# 09 · 测试战役报告（端到端 / 静态分析 / 幂等 / 文档）

> 日期：2026-09-18。范围：`/opt/lat` P0 全模块；真实体系来自
> `/mnt/f/CRO/40_multi_electrolytes/11-systems`（sys1–sys6，178k–188k 原子，含聚合物 PMBA）。

## 1. 端到端回归：lat vs GROMACS（真实电解液–聚合物体系）

方法：`tests/regression/level1_algo/compare_gmx.py` 在运行时用 `gmx` 生成 Ground Truth，
再用 lat 的 Python API（拓扑只加载一次）逐对比较。网格 0.002 nm，rmax 2 nm，帧 90–100（11 帧）。

### 1.1 RDF 与配位数（ref = LI）

| 体系 | 原子 | 最差 max\|Δg\| | rms | mean | 最大 cn 相对误差 | 备注 |
|---|---|---|---|---|---|---|
| sys1 | 183567 | 0.134 | 0.021 | ~1e-5 | 0.16% | 纯电解液 |
| sys2 | 178950 | 0.137 | 0.022 | ~1e-5 | 0.21% | +TTE |
| sys3 | 188607 | 0.155 | 0.020 | ~1e-5 | 0.16% | +MBA |
| sys4 | 188607 | 0.151 | 0.022 | ~1e-5 | 0.16% | **+PMBA 聚合物** |
| sys5 | 183990 | 0.151 | 0.022 | ~1e-5 | 0.20% | MBA+TTE |
| sys6 | 183990 | 0.184 | 0.024 | ~1e-5 | 0.20% | **PMBA+TTE** |

- 各对含 LI–LI / LI–FSI / LI–FEC / LI–FEMC / LI–TTE / LI–MBA / LI–PMBA（按体系存在）。
- 峰值高度最大到 g≈12.7，故 max\|Δg\|≤0.184 相当于 **首峰 < 1.5%**；mean 偏置 ~1e-5（无系统偏差）。
- **最大配位数相对误差 ≤ 0.2%**（在 rmax=2 nm 处）。
- 结论：RDF/CN 在**分子混合物 + 大分子（PMBA）体系**上与 `gmx rdf` 事实上一致。

### 1.2 密度剖面（sys1，100 slabs，z）

| 量 | 结果 |
|---|---|
| max\|Δρ\| | 16.2 kg/m³（≈0.95%） |
| rms | 4.1 kg/m³（0.24%） |
| mean 偏置 | +0.008 kg/m³ |
| 平均密度 | gmx 1717.9 vs lat 1717.9 kg/m³ |

### 1.3 MSD / 扩散（sys1，`-seltype mol_com`，stride 10）

| 对象 | max\|ΔMSD\| | rms | MSD@9000ps gmx | lat |
|---|---|---|---|---|
| LI | 0.0132 nm² | 0.0043 | 0.0368 | 0.0379 nm² |
| FSI | 0.0136 nm² | 0.0044 | 0.0298 | 0.0306 nm² |

- 曲线一致到 **~3–7%**。为此实现了 GROMACS 轨迹的**去 PBC 跳跃连续化**（xtc 无 image flags，
  gmx 内部默认做 nojump；lat 现在逐帧最小镜像重建 `xu/yu/zu`）。

### 1.4 未与 gmx 对拍的模块

`agg/cluster`（gmx 无对应）、`probe`、`pistack`、`occupancy`、`hbond`（有 gmx hbond，但判据已按源码修正，
未做逐值对拍；数据集无合适给体）。这些模块以合成体系数值断言 + 真实样例冒烟为准。

## 2. 静态分析（子代理，对照 gmx 源码 / 文献）

### 2.1 RDF/PBC（子代理 A）

| 条目 | 判定 |
|---|---|
| 归一化 `g = factor·ΣcV/(F·n_a·n_b·shell)`，`cn=Σc/(F·n_a)` | 与 gmx 等价（常数体积/小 bin 极限） |
| 体积加权：lat 算术 `<cV>` vs gmx 调和 `<c>N/⟨1/V⟩` | minor（NPT 0.5% 波动→~2.5e-5） |
| shell：`4πr²Δ` vs 精确 `(4/3)πΔ(r³)` | minor；**已改为精确 shell** |
| abscissa：lat 用 bin 中心，gmx xvg 用左边缘 | minor，是实测 max\|Δg\| 的主因 |
| 同选择自对 | 等价（lat factor=2 对应 gmx 有序对） |
| `-excl`/`-intramol` 不改 N·M 归一 | correct |

**已修**：精确 shell；同集合改用集合相等判定；不同但重叠的选择剔除自身对。

### 2.2 agg / hbond / transport（子代理 B）

| # | 位置 | 严重度 | 问题 | 处置 |
|---|---|---|---|---|
| 1 | `hbond.py` | major | 用了 H 处夹角（D–H···A>150°），gmx 用 D 处 H–D···A≤30° | **已改为 gmx 判据** |
| 2 | `hbond.py` | major | 无 PBC 最小镜像 | **已修（传 box，MIC）** |
| 3 | `hbond.py` | major | ACF 用常数 total 归一→伪线性衰减 | **已修（按有效 origin 无偏归一）** |
| 4 | `transport.py` | major | `D_si` 用 1e-4（cm²/s）却代入 m²/s → σ_NE ×1e4 | **已修为 1e-8** |
| 5 | `transport.py` | major | `frame_interval` 未随 stride/frame_range 缩放 | **已修（`Context.effective_dt`）** |
| 6 | `residence.py` | major | 同上 ACF 归一 | **已修** |
| 7 | `_graph.py` | minor | `-mode center` 的 COM 未做 PBC 展开 | **已修** |
| 8 | `agg.py` | minor | 文档所述边数守恒告警未实现 | 改为「全连通→cutoff 过大」告警 + 去重 |
| 9 | `transport.py` | minor | cumsum 缺 acf[0] 半权 | **已补梯形项** |
| 10 | `_graph.py` | 定义 | 组分法 SSIP/CIP/AGG 与文献「阳离子配位数」口径在 A2B1 等桥接处不同 | **保留**（`agg.md` 明确定义为组分法；见 §2.3） |
| 11 | `pistack.py` | minor | `cos_tol` 死代码；`-level` 未用 | 记录，未改 |

### 2.3 关于 SSIP/CIP/AGG 口径（重要声明）

`agg` 采用**组分法**：连通分量 `A1B0→SSIP`、`A1B1→CIP`、其余→AGG。文献（Efaw 2023）用**阳离子配位数**：
CN=0/1/≥2。二者在无桥接时一致；出现 A2B1 桥接时，组分法把整个 AGG 组分计入，而配位数法会给每个 Li CN=1。
两者物理含义不同，均已公开定义。后续可加 `-label cn|cluster` 开关（记为 P1）。

## 3. 幂等性测试

`tests/unit/test_idempotency.py`：对 `rdf density(1d/2d) gyrate contact rms msd hbond agg cluster
pistack probe(ffv/sasa) residence occupancy` 各运行两次，比较除 `# generated:` 时间戳外的**全部内容**
（数据 + 稳定 header），并验证第二次触发 `.bak.1` 轮转；`select/describe` 比较 stdout 逐字节一致。
**16 passed。**

## 4. 文档 / 用户体验（无上下文子代理，禁读源码）

子代理仅用 `README/docs/config/-h/describe` 操作样例，结论：

| 级别 | 问题 | 处置 |
|---|---|---|
| H1 | `-h` 不自解释（无默认值/取值/示例），与「文档即教材」承诺不符 | **已改**：`-h` 现带 `[default: …]`、枚举 `choices`、epilog 指向 `docs/MANUAL_ALGORITHMS.md` 与 `lat describe` |
| H2 | 文档示例用 `resname LI`，但 LAMMPS 拓扑只有 `type/element/mol` | 待补：在 `03_io_selection.md` 明确 LAMMPS 原生选择子 |
| H3 | WARNING 只进 header，不屏显（违反 M7） | **已改**：CLI 将表警告去重后经 logger 输出到 stderr |
| H4 | 输出路径打到 stdout（违反契约） | **已改**：改走 stderr `notice` |
| H5 | `r_first_min/cn_first_shell` 在稀释样例上误导 | **已改**：始终给 `r_first_peak/g_first_peak`；无有效首极小则省略并 WARN |
| H6 | `msd` 在 wrapped 坐标上仍给 D、R²=1 | **已改**：GROMACS 轨迹自动去跳跃；缺失时 WARN 屏显 |
| H7 | MSD 文档 COM 默认自相矛盾 | 待补：文档统一为 `remove_com=true` |
| H8 | `-fit_range auto` 公式与实测不符 | 待补：文档改为「`[max(1,0.1T), max(+2,0.5T))`」 |
| M3 | `06` 误写 vdW 表单位 nm | 待补：应为 Å |
| M5 | 命名示例含 selection 字面值，与规范冲突 | 待补 |

> 子代理同时确认：`density`、`agg` 结果可被正确物理解释；`rdf` 曲线可解释但工具自带标量需谨慎；
> 稀释样例的 `msd` 本不可用——现已能给出屏显警告。

## 5. 本轮额外发现并修复的缺陷

- **PBC 边界崩溃**：`Box.wrap` 浮点落在 `frac=1.0` 时，`cKDTree(boxsize)` 报「greater than box」。
  真实 183k 体系一运行即崩。已修（fractional clip `1-1e-9`）。
- **hbond 空帧崩溃**：部分帧无氢键时 `np.vstack` 维度错误（E2E 冒烟抓到）。已修。
- **hbond ACF 病态拟合**：不衰减时 `tau_fit~1.7e15`。已加 `0<tau<=T` 否则 `nan`+WARN。
- **性能**：hbond 接受体循环向量化，全样例从分钟级降到 1.2 s。

## 5b. 真实项目兼容性战役（/mnt/d/CurrentProjects + /mnt/f/CRO）

用 `tests/compat/compat_run.py` 对真实项目目录浅扫描（`--maxdepth 3`，自动按**原子数匹配**
配对 data+dump），对每个体系跑 15–16 个模块。D 盘 21 个体系、F 盘 43_reaxff_cu_diag / 45_nbn /
48_og3 / 51_FeCHONP_shear / 53_dehydrate 共 75+ 体系。**修完后全部通过**（剩余失败均为真实数据问题，
见下）。

### 本战役发现并修复的缺陷

| # | 现象 | 根因 | 处置 |
|---|---|---|---|
| 1 | `ValueError: Some input data are greater than the size of the periodic box` | `cKDTree(boxsize=L)` 要求坐标 ∈ [0,L)，而 `wrap()` 返回 [origin, origin+L)（LAMMPS 原点非零即崩）；负坐标同样报错 | 邻居搜索前减 `box.origin`（`neighbor.py`） |
| 2 | 1 帧读取耗时 19 s（39408 原子） | `DumpFrameSource.frames` 过滤后仍读到 EOF；`_count_frames` 逐行扫全文件 | 到最大请求帧即 break；帧计数改 1 MiB 分块 `bytes.count` |
| 3 | 显式 `frame_range` 仍触发全文件帧计数 | `Context` 构造时强制 `n_frames` | 帧计数**惰性化**；`resolve_frames` 对显式 start:stop 不再要求总数 |
| 4 | `unknown atom_style 'full/kk'` | 未剥离 KOKKOS/GPU/OpenMP 后缀 | `st.split("/")[0]` 归一化（`lammps_data.py`） |
| 5 | `could not convert string to float: 'Ni'` | dump 含 `element` 等字符串列 | `_columns()` 容错解析字符串列（`lammps_dump.py`） |
| 6 | probe 在 4096 原子体系 118 s、39408 原子超时 | `_occupancy` Python 三重循环逐体素 | 每原子局部块 numpy 向量化（probe 120 s→2–3 s，快 ~40–60×） |
| 7 | `no vdW radius for ['Mo'/'Te'/'Ti']` | vdW 表只有 26 种元素 | 补全 reax_tools 全周期表（~118 元素） |
| 8 | OVITO data 无 `Masses` → probe/density 报错 | 无 type→element / 无质量 | dump 有 `element` 列时回退推断（`Context.element_map`）；密度错误提示 `-quantity number` |
| 9 | 跨帧原子数变化时错误信息无帧号 | — | 报「frame N has X atoms but topology has Y」并建议 `-frame_range` |

### 被判定为「数据本身问题」的失败（非 lat 缺陷）

- **反应/脱水轨迹原子数逐帧变化**（51_FeCHONP_shear reaxff、53_dehydrate continuous、
  48_og3 Fe3O4）：dump 内 `NUMBER OF ATOMS` 从 3911→3931 波动。lat 明确报错并建议限定窗口。
- **data 与 dump 不属同一阶段**（GCMC 初始 MOF vs 含吸附质快照、ETH.data vs 乙醇/水混合 dump）：
  harness 按原子数过滤后跳过或报清晰错误。

### 进度提示（用户要求）

`Context.frames` 现对每次读取输出（stderr，默认可见）：
`reading N frame(s) from <file>` → `processed frame k/N (t s)`（首帧、末帧、且每 ≥2 s 一次）→
`done: N/N frames in t s`；CLI 另打印 `reading topology <file>` / `parsed X atoms in t s`。
用于区分「卡死」与「仅慢」。

## 6. 建议的后续测试（尚未做）

1. `level1_algo` 常态化：把 `compare_gmx.py` 接入 CI（抽样 2 体系 × LI–FSI 即可）。
2. gmx `hbond`/`gyrate`/`sasa`/`current` 逐值对拍（需合适体系与判据对齐）。
3. `transport` 的 GK 平台与 NE 与 `gmx current` 对拍（需速度单位约定）。
4. 三斜盒真实体系回归；`neighbor` 非 same 分支在 cutoff>半高时的去重。
5. 183k 体系的内存/耗时基准（当前 msd 全 101 帧 >120 s，需流式化）。
