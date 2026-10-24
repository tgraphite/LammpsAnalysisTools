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
pistack freevol(ffv/connectivity) sasa residence occupancy` 各运行两次，比较除 `# generated:`
时间戳外的**全部内容**（数据 + 稳定 header），并验证第二次触发 `.bak.1` 轮转；
`select/describe` 比较 stdout 逐字节一致。**16 passed**（0.2.0 把 `probe` 拆成 `freevol`+`sasa` 后，用例数不变）。

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
   进度（2026-10）：`hbond`/`gyrate`/`rms`/`rmsf`/`rdf` 已完成，`sasa` 已定位缺陷，剩 `current`。
3. `transport` 的 GK 平台与 NE 与 `gmx current` 对拍（需速度单位约定）。
4. 三斜盒真实体系回归；`neighbor` 非 same 分支在 cutoff>半高时的去重。
5. 183k 体系的内存/耗时基准（当前 msd 全 101 帧 >120 s，需流式化）。

---

## 7. LAMMPS→GROMACS 原生对拍战役（2026-10）

### 7.1 动机

第 5b 节的真实项目战役与 `gmx_regress.py` 都只覆盖 **GROMACS 原生输入**。LAMMPS `data`+`dump`
这一半 IO（Atoms 段排序、`Bonds` 段、`make_whole`）**从来没有原生 gmx 参考**——正是这个盲区
让 §7.2 的两个 bug 长期存活（两个 harness 分别在硬编码的 `/usr/local/gromacs/bin/gmx`
上跑不起来，本机实际是 `/usr/bin/gmx` 2022.5；现已改为 `$GMX`→`PATH` 解析）。

### 7.2 新增

- `tests/regression/convert/lammps_to_gmx.py` —— 尚无任何工具能做这件事。
  按 **id 顺序**建 gmx 拓扑（`Bonds` 段引用 atom id，id 恰为 1..N 时可直接映射），
  每个 LAMMPS 分子一个 `moleculetype`，σ 按 `r = 0.5·σ·2^(1/6)` 反解成 lat 的 vdW 表，
  轨迹写多帧 `.gro` 再让 `gmx trjconv` 转 `.xtc`（`.gro` 定义即 nm，绕开单位歧义）。
- `tests/regression/level1_algo/gmx_ref_lammps.py` —— 对拍台：gyrate（整体 + 每分子）、rms、rmsf、
  sasa、rdf（元素对）、hbond，全部在同一批帧上。

### 7.3 结果（4 链 PET，6200 原子，501 帧，盒 80³）

| 量 | lat | gmx 2022.5 | mean\|Δ\| | max\|Δ\| | 相对 |
|---|---|---|---|---|---|
| Rg(all) | 26.8319 | 26.8316 | 0.0003 Å | 0.0004 | 0.00% |
| Rg(每链 ×4) | 25.2776 / 25.3872 / 24.0370 / 25.4916 | 25.2778 / 25.3874 / 24.0371 / 25.4918 | 0.0002 Å | 0.0005 | 0.00% |
| RMSD | 8.2829 | 8.2830 | 0.0001 Å | 0.0002 | 0.00% |
| RMSF | 2.3167 | 2.3167 | 0.0007 Å | 0.0038 | 0.03% |
| hbond（每帧计数） | 6.6028 | 6.6168 | 0.0140 | 1.0 | 0.21% |
| rdf C–O / C–C / O–O（CN） | 48.44 / 120.22 / 46.59 | 48.31 / 119.89 / 46.47 | — | — | 0.26 / 0.27 / 0.27% |
| SASA | 18 316.9 | 16 033.4 | 2 283.6 Å² | 2 917.3 | **14.24%（缺陷）** |

### 7.4 由此发现并修复

1. **键表索引未随 id 排序重排**（`lammps_data.parse`）。Atoms 段被 OVITO 写成非 id 序，
   6476 条“键”只有 880 条是真的，平均最小镜像“键长” 8.14 Å。后果：`make_whole` 连错原子，
   **Rg 高估 3.60×、RMSD 3.50×、RMSF 4.53×**；87 帧 \|ΔRg\|>3 Å 中 87/87 都伴随一次整盒镜像翻转
   （最大键矢量帧间跳变 80.0 Å，而正确键图只有 3.1 Å 的 H 热运动）。
   修复后 `tests/unit/test_bond_atom_order.py` 同时覆盖解析层与 `gyrate` 端到端。
2. **`hbond` 把 C–H 当给体**。聚合物里每个 CH₂/CH₃ 碳都成了给体：41.7 vs gmx 4.3 每帧。
   修复（给体限制为 N/O/F/S/Cl/Br/I）后 501 帧只有 1 帧差 1 条键。
   `sys1_reduced` 只有 C–H，因此样本级 hbond 测试改为「无给体 → 报错」，与 gmx 的
   `No Donors found` 一致；几何分布的空帧路径改由合成 O–H 体系覆盖。

### 7.5 由此发现、并在 0.2.0 修掉（原记"本次不改"）

- `probe -mode sasa` 不用 PBC → 高估 5–17%（`gmx sasa` 与独立周期实现互证，误差 99% 在贴面原子）。
- PSD：EDT 不周期（最大偏差 17.4 Å）、bin 宽随数据变（0.944 Å vs 头部写的 0.5）、
  积分只有 1.5e-4、局部极大按平台计数、只用第 0 帧。
- connectivity：不闭合环面，free 侧默认 6 邻接；FFV 83% 时除“贯通”外无分辨力
  （非最大分量总和 84 Å³ = 盒的 0.0165%）。
- `gyrate -sel all` 在多分子周期体系里的固有歧义（实测 17%）——**保留**，是量的定义问题，不是 bug。

以上四条在 0.2.0 全部落地（SASA/PSD/connectivity 重做 + 周期几何统一 + 模块拆分），
另加了 `contact`/`residence -mode center` 的 wrapped 质心与 `pistack` 的裸环心两处同类缺陷。
详见 `07_decision_log.md` §7 与 `12_0.2.0_release.md`。

### 7.6 顺带发现的轨迹问题（非 lat 缺陷，属数据）

该 PET 的初始结构（OVITO 写）有 **433 对原子间距 < 1 Å，最近 0.112 Å**：
LAMMPS step 0 的 PE = 1.996e16 kcal/mol、T = 1472 K，前 ~10 帧是从 ~4000 K 淬火回 300 K 的瞬态
（gmx 也看到同一现象：SASA 第 0→1 帧从 24 774 涨到 52 115 Å²）。
**任何覆盖全 501 帧的均值都混入了这段瞬态，分析应从 step ≥ 10000（frame ≥ 10）起算。**

---

## 8. LAT 0.2.0 战役（2026-10）：周期几何 + 模块拆分 + examples 复跑

### 8.1 起因

§7 的对拍暴露出：周期几何在 0.1.0 里有**四套并行实现**，三处把 wrapped 坐标当物理坐标。
用户追问"其它模块有没有同类问题"，逐模块审计后确认 6 处（见 §8.2），并拍板：
不留旧版本、接受数值变化、examples 重新生成、没用到 LAT 的项目不管。

### 8.2 审计结论（0.1.0 代码）

| 模块 | 判定 |
|---|---|
| `rdf` | ✅ 正确（`-intramol` 统计成对 MIC 距离，wrapped 坐标下本就正确——**更正** §7 的"可疑"笔记）|
| `density` `occupancy` `select` `describe` | ✅ 只需 wrapped |
| `hbond` `rms` `gyrate` | ✅ 已对拍 gmx（0.21% / 1e-4 Å / 3e-4 Å）|
| `msd` `transport` | ⚠️ 位移接近半盒时有歧义（文档已声明，本轮不改）|
| `agg` `cluster` | ✗ `-mode center` 用"首原子 MIC"求质心（长于半盒的分子会折叠）|
| `contact` `residence` | ✗ `-mode center` 直接对 wrapped 坐标取平均（完全没有 PBC）|
| `pistack` | ✗ 环心/法向用裸坐标（跨面时整盒偏移）|
| `probe` | ✗ SASA 裸 `cKDTree`；PSD/connectivity 网格不闭合环面 |

### 8.3 判据：整体平移不变性

新增 `tests/regression/level2_format/test_pbc_invariance.py`：把体系整体平移一个非格点矢量（+7.3 Å），
逐列比较输出，物理量必须逐位相同。0.1.0 代码下 contact（min/center 的 com_dist：
4.0 vs 6.0，真值 6.0）、residence(center)、pistack 全部失败；修好后全绿。
这条判据比"看数值合不合理"可执行得多。

### 8.4 与 gmx 的二次对拍（回答"gmx 怎么做"）

- **SASA**：0.2.0 周期化后 PET 体系 frame 0/1/2 相对 `gmx sasa` 为 −0.17% / −0.83% / −0.85%
  （0.1.0 为 +14.2%）。
- **hbond 给体**：`gmx hbond` 对 `El_C`（给体组）/`El_O` 输出
  `Found 0 donors and 1128 acceptors` + `No Donors found`；对 `System`/`System` 输出
  `Found 8 donors and 1128 acceptors`（= 4 条 PET 链的 8 个端羟基）。与 lat 的给体规则完全一致。
- **hbond 大体系**：32800 原子的水/表面活性剂体系，lat 13908.8 vs gmx 13808.5（同 40 帧），**+0.7%**；
  而 0.1.0 的存档值 4340.3 与**两个版本都差 3.2×**——该文件是更早的实现产出的，属于必须重算的产物。
- **质心（`contact -mode center`）**：`gmx distance -select 'com of group Mol1 plus com of group Mol2'`
  在 **wrapped 与 whole 轨迹上给出完全相同的值**（1.408564870 nm），即 gmx 在给 `.tpr` 时
  自己 rmpbc；其 `com of group` 是**质量加权**质心。lat 的几何质心与之差 0.0097 Å（mean），
  换成质量加权质心后差 0.0025 Å（max 0.0052 Å，已达 `.xtc` 的 0.001 nm 量化极限）。

### 8.5 交付物

- `src/lat/core/periodic.py`（四层入口）+ `tests/unit/{test_periodic_core,test_no_ad_hoc_pbc,test_config_sync}.py`；
- `lat freevol` + `lat sasa`（`probe` 删除）；
- `docs/12_0.2.0_release.md`（发布说明）、`docs/11_examples_recheck.md`（examples 复跑报告）、
  `examples/_lat_v020_recheck/{recheck.py,MAPPING.tsv,donor_inventory.py}`（可复跑的映射 + 工具）。
