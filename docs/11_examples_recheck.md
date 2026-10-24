# examples 产物复跑报告（LAT 0.1.0 → 0.2.0）

本文件回答一个问题：**`/mnt/d/CurrentProjects/lammps_academia/examples` 里已有的 LAT 产物，
哪些必须重算、重算后变了多少、为什么变、外部工具是否支持新值。**

- 逐文件映射（归档）：`docs/data/0.2.0_examples_mapping.tsv`（84 个 LAT 产物 100% 覆盖）
- 复跑/覆盖/绘图工具（归档）：`tests/regression/examples_recheck/`
  （`recheck.py` 49 条命令：17 条 affected / 32 条 verify，共 73 个产物文件；
  `apply.py` 覆盖回原目录，默认 dry-run，只有 `--apply --delete` 才写/删；
  `donor_inventory.py` 按元素清点"带 H 的重原子"即给体候选；
  `make_before_after.py` 生成修复前后对比图）
- 修复前后对比图：`docs/figures/0.2.0/`（16 张，见 §6）
- 执行状态：**已于 2026-10-17 覆盖回原目录**（73 个文件；删除 5 个 invalid + 1 个残留目录），
  复跑目录本身已删除，项目树里不再保留任何 LAT 开发/测试痕迹。

统计：**affected 25 个产物**（0.2.0 的代码改动会改数值）、**verify 48 个**（命令重跑以证明代码路径未变；
其中 5 个的旧文件本身有缺陷，见 §3）、
**invalid 5 个**（命令化学上不成立，建议删除）、**scaffolding 27 个**（0.1.0 时代的探索残留，建议删除）。

---

## 1. 结论表（affected）

| 产物 | 旧均值 | 新均值 | 变化 | 原因 | 外部证据 |
|---|---|---|---|---|---|
| `31/analysis/run.gyrate` | 96.536 Å | **26.8319 Å** | −3.60× | `Bonds` 按 atom id 解释后未随 id 排序重排 | `gmx gyrate` 26.8316 Å（差 3e-4 Å）|
| `31/analysis/run.rmsd` | 28.974 Å | **8.2829 Å** | −3.50× | 同上 | `gmx rms -nomw` 8.2830 Å（1e-4 Å）|
| `31/analysis/run.rmsf` | 10.484 Å（max 36.4）| **2.3167 Å**（max 8.90）| −4.53× | 同上 | `gmx rmsf -fit yes` 2.3167 Å（7e-4 Å）|
| `31/analysis/run.probe.ffv` → `run.freevol` | 0.828822 | **0.828822** | **逐位不变** | 占用判据本来就按轴 MIC 回绕，只是换了文件名 | — |
| `31/analysis/run.probe.psd` → `run.freevol.psd` | 40 bins、宽 0.944 Å、∫=1.5e-4、主峰 0.53 Å | **69 bins、宽 0.5 Å、∫=1.000、主峰 1.75 Å** | 重定义 | 见 `MANUAL_ALGORITHMS.md` §12.2 | 距离变换对 9×8×7 非均匀网距的 O(N²) 暴力校验（单测）|
| `31/analysis/run.probe.connectivity` → `run.freevol.connectivity` | 30812 行、每帧最多 193 分量、`percolating` 列无分辨力 | **501 行汇总**，11–42 分量、1–3 贯通、最大分量占比 0.99953–1.00000 | 重做 | 环面闭合 + 26 邻接 + 逐帧汇总 | 环面连通分量单测（跨面合并/角落 26 邻接/贯通）|
| `31/analysis/run.probe.sasa` → `run.sasa` | 18316.9 Å² | **15946.5 Å²** | **−12.9%（对 gmx 从 +14.24% 降到 +0.59%）** | 邻域搜索与遮挡判定未用 MIC | `gmx sasa` 16033.4 Å²；逐帧 −0.17%…−0.85% |
| `22/8x/hbond_any.*`（5 个）| n_hbond 4340.3 | **13908.8** | +3.20× | 存储文件来自更早实现；**代码层 0.1.0→0.2.0 只差 +0.36%** | `gmx hbond`（同 40 帧）13808.5 → **+0.7%** |
| `22/8x/hbond_any.hbond.acf` | 11 行（lag 0–10）| **40 行**（lag 0–39）| 表被截断 | 旧文件与 40 帧窗口不符（头部还写 `frame_range: all`）| — |
| `24/licl/any.hbond.*`（5 个）| n_hbond 2369.5 | **6773.1** | +2.86× | 存储文件来自更早实现；**代码层 0.1.0→0.2.0 逐位相同**（该体系 8080 个 H 全部连在 O 上，无 C–H 给体）| `gmx hbond` 6087.2；差 **+13.8%，来源是受体集合（见 §4）** |
| `24/licl*/{li,cl}_go.contact`（4 个）| com_dist 见 §2 | **+0.007 Å** | 质心改为整分子重建 | 只有 `com_dist` 变；`n_contacts`/`min_dist` 逐位不变 | `gmx distance`（`com of group`）PET 体系差 0.0097 Å（几何质心）/0.0025 Å（质量加权）|
| `24/licl*/{li,cl}_gr.contact`（4 个）| 同上 | **逐位不变 + 新增 WARN** | 石墨烯层是二维周期共价网络，跨盒时无自洽像 | 见 §2 | `gmx distance` 在该体系**直接报错**（`inconsistent shifts over periodic boundaries`）|
| `22/8x/oil.density.2d` | 0.635 kg/m³ | **44.472 kg/m³** | **+70×** | 旧文件少了**一个盒长**的归一因子（44.328/0.635 = 69.8 ≈ L）| 解析值 mass/V = 44.328 kg/m³ |
| `22/8x/surfactant_density_2d.density.2d` | 2.772 kg/m³ | **194.187 kg/m³** | +70× | 同上 | 解析值 193.561 kg/m³ |

---

## 2. 质心类分析（contact / residence `-mode center`）：gmx 到底怎么做

用户要求"本轮结束后汇报 gmx 的实现"。用 PET 体系（4 条链跨越盒面，`gmx` 可用）与 24 号体系
（`gmx` 不可用）分别实测：

**gmx 侧**

1. `gmx distance -select 'com of group "Mol1" plus com of group "Mol2"'`
   在 **wrapped 轨迹与 `-pbc whole` 轨迹上给出完全相同的值**（1.408564870 nm，平均绝对差 0）
   → gmx 在给了 `.tpr` 时会自己 `rmpbc`（删掉 tpr 就不会），等价于 lat 的 `make_whole`。
2. 该值是**质量加权质心**：lat 的几何质心差 0.0097 Å（mean），
   把 lat 换成质量加权后差 **0.0025 Å（max 0.0052 Å）**——已达 `.xtc` 的 0.001 nm 量化极限。
3. gmx 没有"min/center"两套模式：`gmx pairdist` 用 `-refgrouping atom|res|mol` 决定
   逐原子 MIC 距离还是成组质心距离；两者都走 MIC + rmpbc。
4. **二维周期共价网络（电极/凝胶）gmx 直接失败**：24 号体系上
   `gmx distance` 与 `gmx trjconv -pbc whole` 都报
   `inconsistent shifts over periodic boundaries in a molecule type consisting of 14456 atoms`。
   即：这种分子图无法自洽地"摆正"，gmx 拒绝给数。

**lat 侧（0.2.0 的策略）**

- 选区是**整分子**且键图可自洽重建 → 沿键图复原后求质心（精确）。
- 键图**跨盒不自洽**（环绕盒一圈）→ 退回 wrapped 质心 + **WARN**，明确写出
  "沿被回绕轴的位置无意义，gmx 会拒绝该情形"。
  24 号体系的 `type 9`（石墨烯层，960 原子一个分子）正是这种情形，所以 `_gr` 的 4 个产物
  **数值逐位不变**（旧值本来就是 wrapped 均值），只是多了 WARN；`_go`（GO 片层，分子不跨盒不自洽）
  则改用精确质心，com_dist 变 **+0.007 Å**。
- 选区**不是整分子**（如 `type 1-8` 里的部分原子）→ 用"首原子 + MIC"质心，
  并按参考原子依赖度给 WARN（阈值 0.1 Å）；实测 examples 里没有触发。
- 一致性检测函数 `lat.core.periodic.bond_graph_consistent()` 用**与 gmx 相同的判据**
  （重建后每条键的长度是否等于其 MIC 长度），单测覆盖"跨盒环 → False / 去掉闭合键 → True"。

---

## 3. verify：未改动代码路径（多数逐位不变，例外见末两行）

48 个产物的**命令**在 `_lat_v020_recheck/` 下重跑，与原文件逐值比较（按位置比较，忽略列名）。
这一组里 lat 的代码路径与 0.1.0 相同，所以"数值不变"是预期；**变了的 5 个都说明旧文件本身有问题**：

| 结果 | 数量 | 说明 |
|---|---|---|
| **逐位完全相同** | 32 | 含全部 `rdf`**数值**、`agg`（3 个，含 3320 行的 speciation）、`cluster.cn`、`msd`（3 个）、`density`（24 号 8 个 + 32 号 2 个）|
| 列名不同、数值逐位相同 | 11 | 旧文件用 `Water (O) - Water (O)` / `Li - Solv` 之类的人类可读列名，现为 `g_O_O` / `g_Li_type_6-18`；**值 max_rel = 0.000e+00** |
| 行数不同（旧文件不完整）| 3 | `23/li_solv.rdf.cn` 100 行 vs 200（旧文件与自己的 `nbins=200` 头不符）、`23/li_fsi.rdf.cn` 99 行 vs 200、`23/li_solv.cluster.size` 9 行 vs 10 |
| **数值不同（旧文件归一化错）** | 2 | `22/8x` 的两个 2D 密度（见 §1 末两行）：**旧值是解析值的 1/70**，新值 = 解析值 |

其中 `li_solv.cluster.size` 的 9 vs 10 行就是"旧文件缺 `size=1` 一档"：它把 34502 个
单分子/单体素分量算进分母却不出行 → `p` 之和只有 **0.2337**；新文件 `p` 之和 = **1.0000**。

> 列名变化与那 3 个"行数不同"的旧文件都来自更早的 lat 版本（`git log` 只有一个 commit，
> 无法逐行复现当时的代码）。**新文件彼此自洽**（`nbins` 与行数一致、分布积分 = 1），
> 因此建议按新文件覆盖。

---

## 4. hbond：两个必须分开讲的问题

1. **给体定义**（0.1.0 的 `elem[i] == "H"` 使每个 C–H 都成为给体）。
   实测影响：
   - PET 体系（链状，端羟基）：0.1.0 41.7/帧 vs gmx 4.3/帧 → 修好后 6.6028 vs 6.6168（**0.21%**）；
   - 22 号体系（水 + 表面活性剂）：旧代码 13885.7 vs 新代码 13840.3（**−0.36%**，
     C–H 只带来 45 条额外氢键），gmx 13808.5（新代码 **+0.7%**）；
   - 24 号体系：8080 个 H **全部**连在 O 上（`donor_inventory.py` 实测），
     两种定义下给体数都是 4080 → **逐位相同**。
2. **受体集合**（0.2.0 与 0.1.0 相同，`_ACCEPTORS = {N,O,F,S,Cl,Br,I}`）。
   `gmx hbond` 在这三个体系里都把受体算成 **O 一个元素**（PET 1128、22 号 8480、24 号 4136），
   而 lat 在 24 号体系还接受 **144 个 Cl⁻** → 6440 个受体对多出约 **+13.8%** 氢键。
   即 24 号的 `any.hbond` 与 gmx 的差**全部来自受体集合**，与给体修复无关。

**建议（留给用户决策）**：lat 的受体集合（含卤素/S）在化学上更正确（Cl⁻ 是经典氢键受体），
gmx 的 N/O 集合是蛋白取向的约定。可选：
(a) 保持 lat 口径，在文档和产物头部注明"含卤素体系不可与 `gmx hbond` 直接比数值"（**推荐**，已写入 §5）；
(b) 加 `-acceptors gmx|polar` 开关，默认 `polar`，需要对齐 gmx 时切 `gmx`。

---

## 5. 结论与建议动作

1. **25 个 affected 产物按新值覆盖**（`apply.py --apply`）。其中 6 个 `_gr` 类只有 WARN 变化，
   数值逐位不变——覆盖与否都不影响结果，覆盖后能把"此质心沿回绕轴无意义"的告警固化下来。
2. **5 个 invalid 产物建议删除**：`22/8x/hbond_surfactant_water.*`（实测 `mean_n_hbond=19.9`，
   即 1920 个 C–H 给体只凑出 ~20 条"氢键"，全部是 C–H…O 伪迹）。
   `-sel1 type 2 3 4 5 6` 里有 1920 个 C–H、0 个电负性 H，lat 0.2.0 直接报错，
   `gmx hbond` 对同样的组输出 `Found 0 donors ... No Donors found`。旧文件是 C–H…O 伪迹。
3. **27 个 scaffolding 文件**（`31/analysis/latfix.*`、`latcheck.*`、`run.{ffv,psd,sasa}.out.csv`、
   `latcheck_bond_order.patch`、`gmx_convert3.py`、`bonds_fixed.tmp.npy`）：0.1.0 调查期的中间产物，
   会话中已随复跑清理删除（它们不在 2026-10-16 的 examples 快照里，不可恢复也不需要恢复）。
   仍在原地的同类残留只有一个目录 `31/analysis/lat_stage/`（`/opt/lat` 的临时镜像），
   以及 3 个失去 CSV 的孤儿图 `run.{ffv,psd,sasa}.out.png`——都列在 `MAPPING.tsv` 里，
   `apply.py --delete` 会一并清掉。
4. **没有用到 LAT 的目录不动**：`11_cu_minimal`、`12_cu_tensile`、`13_sio2`、`14_pure_water`、
   `toassign_*` 里的 `.out.csv` 都是 LAMMPS `fix print` 或第三方工具的输出（首行不是 `# program: LAT`），
   不在本次范围。
5. `32_molecular_crystal/plot_commands.sh` 里的 `lat probe ... -sel "type 1-7"` 已失效
   （`probe` 删除；且网格模式本来就不接受 `-sel`），需改成 `lat freevol -grid_dr 0.25`
   （该命令从未产出过文件，不影响已有产物）。

**已执行（2026-10-17）**：

```bash
python3 apply.py --apply --delete     # 73 copy, 6 delete（5 个 invalid + lat_stage/）
# 之后：重命名 4 个模块改名产物（run.probe.* → run.freevol* / run.sasa），
#       用 `lat plot1d`/`plot2d` 重画**全部 50 张** LAT 一维/二维图，
#       改写 73 个 CSV 头里的 -o 路径，删除 examples 下的复跑目录与全部开发痕迹。
```

重命名产物（0.1.0 → 0.2.0）：`run.probe.ffv.out.csv` → `run.freevol.out.csv`、
`run.probe.psd.out.csv` → `run.freevol.psd.out.csv`、
`run.probe.connectivity.out.csv` → `run.freevol.connectivity.out.csv`、
`run.probe.sasa.out.csv` → `run.sasa.out.csv`（同名 `.png` 一并改）。

## 6. 修复前后对比图（`docs/figures/0.2.0/`）

每张图把 0.1.0（红）与 0.2.0（蓝）画在同一坐标系里，标题给出变化原因：

| 图 | 看点 |
|---|---|
| `31_..._run.gyrate.png` | Rg 从 96.54 Å 掉到 26.83 Å（键表 id 重排）|
| `31_..._run.rmsd.png` | RMSD 28.97 → 8.28 Å |
| `31_..._run.rmsf.png` | RMSF 逐原子锯齿（36 Å 峰值）→ 与结构对应的 2.32 Å |
| `31_..._run.probe.sasa.png` | SASA 18317 → 15946 Å²（贴面原子丢邻居）|
| `31_..._run.probe.psd.png` | PSD 从 0.944 Å bin、∫=1.5e-4 变成 0.5 Å bin、∫=1 |
| `31_..._run.probe.connectivity.png` | 分量数/帧 193 → 42，并给出 frac_largest≈0.99996 |
| `31_..._run.probe.ffv.png` | 对照组：逐位不变 |
| `22_..._hbond_any.hbond.png` / `.acf.png` / `.dist_xy.png` | 氢键数 4340 → 13909（gmx 13808）；ACF 从 11 行补到 40 行 |
| `22_..._oil.density.2d.png` / `surfactant_density_2d...png` | 共用色标下旧图几乎全空（少一个盒长因子，×70）|
| `23_..._li_solv.cluster.size.png` | 旧图缺 size=1 一档（p 之和 0.2337）|
| `23_..._li_solv.rdf.cn.png` / `li_fsi.rdf.cn.png` | 旧表只有 100/99 行，新表 200 行 |
| `24_..._any.hbond.png` | 氢键数 2369 → 6773（与 gmx 的差全部来自受体集合 Cl⁻）|
