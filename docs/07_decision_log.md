# 07 · 决策日志与优先级

> 本文件记录已确认决策（开发前冻结）与优先级。所有「待定」项已在此关闭。

## 1. 优先级定义（用户口径）

| 级别 | 含义 |
|---|---|
| **P0** | 一会就要开发（所有已确认项） |
| **P1** | 待验收 P0 后开发（必要性 / 是否重复功能尚需判断） |
| **P2** | 以后验证了关键点再开发（算法或定义不确定） |

**全部已确认功能 = P0。** 唯一调整：
- `gyrate`(rg)、`contact` = **P0**；
- `order` = **P2**（参考轴定义不确定）；
- 以分子体积自动生成 `occupancy` zone = **P2**（多孔材料控制区指认问题）；
- 其余「必要性/重复性」未决项 = P1（当前为空，待 P0 验收后审计）。

## 2. P0 模块清单

`rdf` `density` `hbond` `rms` `msd` `agg` `cluster` `pistack` `probe` `residence`
`transport` `occupancy` `gyrate` `contact` `select` `plot1d` `plot2d` `describe`
+ `core`（io/topology/selection/pbc/neighbor/grid/correlate/output/errors/context）。

## 3. 已确认决策

### 单位 / 时间
- 默认长度 **Å**（real 与 metal 都是 Å；lj 太小众不默认）。未显式给 `-unit` → 报 **INFO**，按默认 real/Å。
- 时间 **全部帧号**，忽略 LAMMPS 时间。寿命/滞留以帧输出，标「k × 帧间距」。

### 输入 / 格式（含源码依据）
- `atom_style` 数据列（LAMMPS 2023.2 `fields_data_atom`）：
  `atomic`=id type x y z(5)；`charge`=id type q x y z(6)；
  `molecular/bond/angle`=id molecule type x y z(6)；`full`=id molecule type q x y z(7)。
  **5→atomic、7→full 可确定；6 列歧义**（charge vs molecular）→ 启发式判断，仍不确定 ERROR 要求 `-style`。
- **不支持多文件/glob 拼接**：一轨迹一模拟，用户自行保证。
- `frame_range`：`start:stop:stride`，默认全帧。
- MSD/transport 解包：`xu yu zu` > `ix iy iz` > 最小镜像增量+WARN。
- **dump 原子顺序按 `id` 对齐拓扑**：dump 顺序是 LAMMPS 内部顺序，与 data 文件（按 `id` 升序）不一致；
  逐帧按 `id` 重排，否则 `type` 与坐标错配、g(r) 失真（见 13_ethanol_water 案例）。
  缺 `id` 列 → 假定一致并 WARN；`id` 集合不匹配 → ERROR。
- **hbond 几何中点**：`mid = A + 0.5·MIC(D−A)` 再 wrap；曾误写为 `0.5·(MIC(D−A)+A)`，
  导致二维分布被压到 1/4 盒（见 14_surfactant 案例）。
- PBC 支持正交/三斜；`-pbc` 隐藏默认开。
- **GROMACS 可内部直读，不对外暴露**（供回归）。

### 拓扑 / 元素
- **不推断分子**；分子内对需用户提供 `mol` 列或 GROMACS 残基编号。
- **rdf 默认排除同分子内对**（不计算分子内 RDF）；`-intramol` 打开（包含分子内）。
  有键但无分子 ID 时无法排除 → WARN，须显式 `-intramol`。
- `type→element`：默认按质量 **±0.2%** 匹配元素表；多解/无解 → **ERROR 中断**。
  行尾注释不可靠，不采用；可用 `-types` 文件显式覆盖。
- vdW 半径表采用 `reax_tools` 的 `ELEMENT_ATOMIC_RADII`（Bondi 型），整表进 `lat.yaml:vdw`；
  probe 用纯 vdW，键感知才用 `0.6*(ri+rj)`。

### 选择器
- 静态 MDAnalysis 语义 + `ndx`；无交互输入。
- 读入 `-index`/`-n` 时，`-sel1 -sel2` 等**允许是 ndx 中的组名**。
- **selection 字面值只存在于输出文件头**（因其可能为 `file.ndx -> water` 或 `type 1 2 3`）。
- `lat select` 打印命中统计，可导出 ndx。

### 算法口径
- **绝对轴分布（density 1D/2D、hbond 几何）**：照 **gmx density 口径**——逐帧按**分数坐标**分箱、
  轴按**平均盒长**标定；相关轴盒长偏离均值 **>±2% WARN**。非正交/畸变盒的截面近似由用户负责。
- **平滑只在绘图**：任何计算命令都输出**原始结果**（含 density；ACF 算不出寿命也照常输出）；
  `-smooth {none,gaussian,savgol}` 只在 `plot1d`/`plot2d` 生效。
- **rdf**：`cutoff` 默认 **10 Å**（不随盒子膨胀），n bins 200；`sel1≡sel2` 默认不算自对。
- **hbond**：对齐 gmx；**D···A 3.5 Å**（`-da` 默认）、H–D–A **>150°**；两组须完全相同或完全不重叠。
- **agg/cluster**：`cutoff` **必填**（不设默认）；真·残基级二部图；边数守恒 WARN。
- **probe**：默认 **H2 探针，半径 1.445 Å**（H2 动力学直径 2.89 Å 的一半）；PSD 用探针-网格近似并声明。
- **transport**：强制 `lammpstrj` 含速度列，否则 ERROR；默认去质心漂移；NE + GK 都给。
- **occupancy**：区域语法 `-zone x/y/z:min-max[,min-max][ ...]`，支持多维度约束；边界半开 `[a,b)`。
- **msd**：默认 3D、去质心漂移。

### 输出
- 命名：`<basename>.<tag>.<mode>.out.csv`，命名即意图。
  例：`mysystem.rdf.cn.out.csv`、`mysystem.rg.out.csv`、`mysystem.rdf_Li_FSI.out.csv`。
- **CSV 头冻结**：`program/module/command` + `selections/definitions/frame_range/scalars/warning`
  + **绘图文本契约** `title/xlabel/ylabel/columns`（`columns` 的 y 名即 legend）。
- **`plot1d` 输入用 `-f`**，不保留位置参数；多图/分组模式已放弃（不记入文档）。
- **`plot1d` 列选择/柱状**：`-columns` 按名只画部分 y 列（x 恒为第一列，永远绘制，不可选为 y）；
  `-bar` 分组柱、`-stack` 归一化堆积（占据比、y 轴 0–1）；`-bar`/`-stack` 互斥、仅限 `plot1d`、不平滑。
- **同名冲突**：旧文件改名 `.bak.1`、`.bak.2`…（不覆盖、不报错）。
- `-debug`：额外落完整 JSON 审计；**对外文档不暴露 debug 模式**。
- 不用退出码；分级 **INFO/WARNING/ERROR**（ERROR 中断），全部到 stderr。

### 配置 / 工程
- 有且仅有程序根目录下的**一个 `lat.yaml`**（开发人员用）；用户可配置的物理参数全部有命令行开关。
- 计算/绘图解耦（`plot1d`/`plot2d` 只读 CSV）。
- 回归：两级输入（算法 GT + 格式兼容）；**开发期重装 InterMol 或无头 VMD** 做 gmx→lammps 转换。

## 4. 遗留技术前置

- **gmx→lammps 转换链：已完成**，改为自研 MDAnalysis writer，无需 InterMol/VMD。
  脚本 `tests/regression/convert/`，样例 `tests/regression/samples/sys1_reduced/`（自检全通过）。
- 转换后原子顺序/type 编号的选择器映射：已随样例落盘（`*.mapping.csv`、`*.selections.yaml`）。
- CI 规模（全量 vs 抽帧子集）——开发时实测后定，只要前后行为意图一致。

## 5. 明确不做

能量/自由能、PMF、MMPBSA、对接、反应网络、固体/金属力学、DXA/位错、SDF、谱(VDOS/SAXS)、
phonon、相场、DPD/CG、AIMD/QM、PKA。

## 6. LAMMPS→GROMACS 原生对拍（2026-10）与其后的修复

起因：此前**没有任何 LAMMPS→GROMACS 的转换路径**（`tests/compat/gmx_regress.py` 只吃 GROMACS
原生 tpr/xtc），所以 LAMMPS data 解析器、原子顺序、键表这一半 IO 从来没有原生 gmx 参考。
新增 `tests/regression/convert/lammps_to_gmx.py`（转换器）与
`tests/regression/level1_algo/gmx_ref_lammps.py`（对拍台：gyrate/rms/rmsf/sasa/rdf/hbond）。

### 6.1 已修（本次固化）

| # | 缺陷 | 影响 | 修法 | 回归 |
|---|---|---|---|---|
| 1 | `lammps_data.parse` 里 `bonds` 建在 **Atoms 段行号**空间，随后其余逐原子数组按 **id 升序**重排，键表没有跟着重排 | Atoms 段非 id 排序时（OVITO 写的文件）键图变成随机重连：6476 条“键”只有 880 条是真的，平均最小镜像“键长” 8.14 Å。`make_whole` 据此摆正跨盒分子 → **gyrate 高估 3.6×、RMSD 高估 3.5×、RMSF 高估 4.5×**，并产生“Rg 突变后又回落”的假象（87 帧 \|ΔRg\|>3 Å，其中 87/87 都伴随一次整盒长度的镜像翻转） | `lammps_data.py` 中把同一置换施加到 `bonds` 上 | `tests/unit/test_bond_atom_order.py`（旧代码下 Rg=9.0 而非 1.0，必失败） |
| 2 | `hbond._donors` 把**任何带 H 的重原子**当给体（未限制元素） | 聚合物里每个 CH₂/CH₃ 碳都成为给体，氢键计数比 gmx 高约一个量级（PET：41.7 vs 4.3/帧） | 给体限制为电负性集合 `{N,O,F,S,Cl,Br,I}`（与受体同集）；无给体时保持 ERROR，对齐 gmx 的 `No Donors found` | `tests/unit/test_bond_atom_order.py` 之外的 3 个样本测试改写为“样本无给体 → 报错”+ 合成 O-H 空帧几何测试；gmx 对拍逐帧完全相同 |

判据：`tests/regression/level1_algo/gmx_ref_lammps.py` 在 4 链 PET（6200 原子、501 帧）上，
修正后与原生 gmx 2022.5 的逐帧差为 Rg(all) 3e-4 Å、每链 Rg 2e-4 Å、RMSD 1e-4 Å、RMSF 7e-4 Å、
hbond **0**（完全一致）。

### 6.2 已确认但本次不改行为（**0.2.0 已全部落地，见 §7**）

- `probe -mode sasa` 不用 PBC → 绝对值高估 5–17%（原生 `gmx sasa` 与独立周期实现互证）。
- PSD：距离变换不周期、bin 宽随数据变（0.944 Å vs 头部写的 0.5）、不作体积归一（积分 1.5e-4）、
  局部极大按平台计数、只用第 0 帧。
- connectivity：不闭合环面、free 侧默认 6 邻接、在 FFV 83% 的体系里除“贯通”外无分辨力。
- `gyrate -sel all` 在多分子周期体系里的固有歧义（可达 17%）。

改这些会改变已发布的数值，属于行为变更，不跟随本次“键序 + 给体”两个 bugfix 一起做。

### 6.3 对拍时踩到的工具口径（写进 `gmx_ref_lammps.py` 注释，避免重复踩）

- `gmx gyrate/rms/rmsf/sasa` 在给了 `.tpr` 时**自己会 rmpbc**（等价 make-whole）；给 `.gro` 当 `-s`
  就不会。所以 gmx↔lat 必须在“都给 tpr”的口径下比。
- `gmx rms` 默认 `-mw yes`（质量加权叠合），lat 的 Kabsch 无权重 → 必须 `-nomw`，
  否则 C/H/O 聚合物上系统性差 ~1.9%（加 `-nomw` 后降到 0.001%）。
- `gmx rmsf -fit` 是布尔（yes/no），不是 `rot+trans`。
- `gmx sasa` 的 vdW 半径是 `r = 0.5·σ·2^(1/6)`，转换器据此把 σ 反解成 lat 的 vdW 表。
- 转换器按 **id 顺序**建 GROMACS 拓扑（Bonds 段引用 id，id 恰为 1..N 时可直接映射），
  每个 LAMMPS 分子一个 `moleculetype`；`[ molecules ]` 不能写成“一个 6200 原子的分子 ×4”。
- 轨迹写成多帧 `.gro` 再让 `gmx trjconv` 转 `.xtc`：`.gro` 定义就是 nm，绕开 MDAnalysis
  XTC writer 的基准单位换算（本地实测会把坐标再除以 10）。

---

## 7. LAT 0.2.0：周期几何统一 + 模块拆分（2026-10）

用户在 §6.2 的四个"待决策"项上拍板：**不保留旧版本、接受 SASA 数值变化并重新生成产物、
跨盒非整分子的质心先 WARN（本轮结束后汇报 gmx 的做法）、examples 先写到新目录并保留映射、
没用到 LAT 的项目不管**。据此完成下列改动。

### 7.1 架构决策：周期几何只有一个入口

0.1.0 有**四套并行**的周期几何实现，其中三处把 wrapped 坐标当物理坐标用。0.2.0 新增
`src/lat/core/periodic.py`，把分子级（`whole_positions`/`residue_com(s)`）、原子对级
（`periodic_neighbors`/`periodic_nearest`/`periodic_vector`/`periodic_distance`）、
体素网格级（`grid_axes`/`voxel_occupancy`/`periodic_distance_transform`/`label_periodic`）、
曲面级（`sasa_shrake_rupley`）四层收成一个入口。

守卫测试 `tests/unit/test_no_ad_hoc_pbc.py` 禁止任何模块出现
`cKDTree` / `distance_transform_edt` / `ndimage.label` / `generate_binary_structure` /
`np.pad` / 裸 `minimum_image`，并强制 PBC 敏感模块按下表引用符号：

| 模块 | 必须使用 | 原因 |
|---|---|---|
| `_graph.py`（agg/cluster 的 center 模式）| `residue_coms` | 原实现 `ref + MIC(pos-ref)` 对长于半盒的链会折叠 |
| `contact.py` | `residue_com` `periodic_distance` | 原实现直接对 wrapped 坐标求平均 |
| `residence.py` | `residue_com(s)` | 同上 |
| `pistack.py` | `whole_positions` | 环心/法向原用裸坐标，跨面即整盒偏移 |
| `hbond.py` | `periodic_vector` | 统一入口（几何原本已用 MIC，无行为变化）|
| `freevol.py` / `sasa.py` | 网格层 / 曲面层 | 见 7.2 |

`rdf.py` **故意不在表里**：`-intramol` 统计的是**成对 MIC 距离**，在 wrapped 坐标上本来就正确。
这一条是本次审计对 §6 笔记的**更正**（原先记为"可疑"）；
`tests/regression/level2_format/test_pbc_invariance.py` 用"整体非格点平移后逐列不变"固定了它，
同时该测试在旧代码下会在 contact（两种模式）、residence(center)、pistack 上失败。
判据是平移不变性，而不是"看起来对不对"。

### 7.2 模块拆分：`probe` → `freevol` + `sasa`

`probe` 把两类算法绑在一起：体素网格类（ffv/psd/porosity/connectivity，需要网格、O(N_voxel)）
与曲面类（sasa，Shrake–Rupley，不需要网格），只共享"vdW + 探针半径"这一行输入。
`-mode all` 因此必须两边都跑，而 PSD/connectivity 极少用（PSD 每帧还要再算一次距离变换，
PET 实测 6.6 s/帧 vs FFV 0.5 s/帧）。

- 新增 `lat freevol -mode {ffv,psd,connectivity}`（默认 ffv）与 `lat sasa`；
- **删除 `lat probe`，不留兼容壳**：`porosity` 与 `ffv` 数值完全相同（0.1.0 里两列同值），
  故不再作为独立模式；传入旧名会得到 `unknown module`；
- 文件命名：`run.probe.ffv.out.csv` → `run.freevol.out.csv`（主表用模块名），
  `run.probe.psd.out.csv` → `run.freevol.psd.out.csv`，`run.probe.sasa.out.csv` → `run.sasa.out.csv`；
- 网格模式不再接受 `-sel`（0.1.0 收下 `-sel` 但在网格模式里**静默忽略**；现在直接报错）。

### 7.3 PSD 重定义（数值会变）

按 Gelb–Gubbins 口径重写：距离变换对 **vdW 占据**做（$r=d$，孔半径 = 能放下的最大球半径，
0.1.0 用 $d-r_{probe}$ 多减了一个探针半径）；可及性用**加探针**的占据；局部极大用
26 连通**平台去重**；每个自由格点按**周期最近孔心 Voronoi** 归属并归一化，使 $\int dV/dr\,dr=1$；
bin 宽固定为 `grid_dr`；逐帧累加（文件头 `frame_range` 与实际一致）。

### 7.4 connectivity 重做

`label_periodic` 闭合环面（wrap 补齐 + "同一物理格点"标签合并），默认 **26 邻接**
（6 连通固体的补集本该用 26），输出改为逐帧汇总（`n_components / n_percolating /
v_largest / frac_largest / v_percolating`），`-components` 才写逐分量表。

### 7.5 行为变更与硬错误策略

- 网格类分析要求**正交且固定盒**：三斜直接报错；盒漂移 >2% 报错，0.1%–2% WARN 并把网格钉在首帧。
  0.1.0 用首帧的网格配每帧的盒子，等于悄悄改体积。
- `contact`/`residence -mode center`：选区不是整分子时，用"首原子 MIC"质心并报出**参考原子歧义**
  （>0.1 Å 才 WARN）。要精确就选整分子；不静默给近似值。
- SASA 的 `-sel` 语义明确为"单独报告哪部分原子的贡献"，表面始终用全部原子构建，
  输出 `sasa_total` + `sasa_sel` 两列。

### 7.6 工程收尾

- 版本 0.2.0（`__init__.py` / 两份 `lat.yaml` / `pyproject.toml`）；`load()` 在
  `meta.version` 与包版本不一致时 WARN；`tests/unit/test_config_sync.py` 断言
  `config/lat.yaml` 与 `src/lat/data/lat.yaml` **逐字节相同**（0.1.0 时二者已漂移）。
- 删除陈旧副本 `build/`（其 `lammps_data.py` 与 `src/` 不一致，是个陷阱）、所有 `__pycache__`、
  `.pytest_cache`；删除死配置键（`meta.units.energy`、`io.data_style`、`io.dump_fields`、
  `grid.nbins_1d`、`grid.dr_3d`、`grid.max_bins_axis`）与死代码 `grid.axis_bins`。
- 测试里的硬编码版本号改为读 `lat.__version__`。
- 转换器 `tests/regression/convert/lammps_to_gmx.py` 修复三处真实缺陷（它是对拍工具，
  自身错了会给出错误的"真值"）：Atoms 段列布局按 `Atoms # <style>` 注释解析（原先只认 5 列 atomic /
  7 列 full）、忽略 `Pair Coeffs` 等额外段（原先会被当成质量行）、`.gro` 原子名只用元素符号
  （原名 `C12345` 超过 5 字符会破坏定长格式）、盒与坐标统一取自**轨迹首帧**（原先用 data 文件
  的盒 + 轨迹的坐标，NPT 体系会静默缩放），`-pbc whole` 失败时回退到 wrapped 轨迹并 WARN。

### 7.7 用户决策（2026-10-17 收口）

| 编号 | 决策 | 结论 |
|---|---|---|
| D1 | `probe` 是否留兼容壳 | **删除**，不留壳；`porosity`（与 ffv 同值）不再单列 |
| D2 | SASA 数值变化 | 接受（对 gmx 从 +14.24% 变 +0.59%），产物全部重算 |
| D3 | 跨盒非整分子的质心 | **先 WARN**；本轮已实测 gmx 的做法（`rmpbc` + 质量加权质心；二维周期网络直接报错），见 `docs/11_examples_recheck.md` §2 |
| D4 | examples 产物如何落地 | 先写独立测试目录 + 保留映射，确认后覆盖回原目录（已执行，见 §7.8）|
| D5 | 未使用 LAT 的项目 | 不动 |
| **D6** | `hbond` 受体集合是否对齐 gmx（gmx 只认 O）| **保持现状**：lat 用 `{N,O,F,S,Cl,Br,I}`，在 `MANUAL_ALGORITHMS.md` §12.5/§17 与 `12_0.2.0_release.md` §6 明确声明"含卤素/S 的体系不可与 `gmx hbond` 直接比数值"（LiCl 体系 +13.8%）|
| D7 | examples 里的开发/测试痕迹 | **全部清除**：复跑目录、gmx 参考目录、调查期脚本与图、`LAT_0.2_PLAN.md` 一律删除；工具与映射归档到 `tests/regression/examples_recheck/` 与 `docs/data/` |

### 7.8 examples 覆盖执行记录（2026-10-17）

- 覆盖 **73** 个产物（25 affected + 48 verify），删除 **5** 个无效产物（`hbond_surfactant_water.*`，
  sel1 无电负性 H）与 1 个残留目录（`31/analysis/lat_stage/`）。
- 4 个模块改名产物随之改名：`run.probe.{ffv,psd,connectivity,sasa}.out.csv` →
  `run.freevol{,.psd,.connectivity}.out.csv` / `run.sasa.out.csv`（同名 `.png` 一并改）。
- 16 张受影响的图用 `lat plot1d`/`plot2d` 依据新数据重画；73 个 CSV 头里的 `-o` 路径改回正常形式
  （复跑时指向测试目录）。
- 修复前后对比图 16 张存 `docs/figures/0.2.0/`（`tests/regression/examples_recheck/make_before_after.py` 可重生成）。
- `32_molecular_crystal/plot_commands.sh` 里失效的 `lat probe ... -sel` 行改为 `lat freevol -grid_dr 0.25`。
- 项目树 `lammps_academia/examples` 不再包含任何 LAT 开发/测试痕迹（复跑目录、gmx 参考、调查脚本/图、
  计划文件全部删除；`*.out.csv.bak.N` 中仅保留用户 2026-09/10 自建的那批）。
