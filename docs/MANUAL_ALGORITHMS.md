# lat 算法原理手册（v0.1）

> 目标：读这一份文档即可理解 lat 每个模块**在算什么、怎么算、为什么这么算**。
> 本文档是教材而非罗列：先讲共同地基，再逐模块讲原理、公式、实现口径与陷阱。

---

## 0. 共同地基

### 0.1 数据模型
一次计算 = **拓扑**（类型/质量/电荷/键/分子）× **轨迹**（逐帧坐标）。
轨迹流式逐帧读取（内存里只保留当前帧），跨帧统计在遍历中累加。时间轴用**帧号**，
不解析文件里的物理时间；寿命等量的单位是“帧”，在文件头注明“1 帧 = 配置的 frame_interval ps”。

### 0.2 周期性边界与最小镜像（PBC/MIC）
盒子用 3×3 矩阵 `H` 表示（列向量为 a,b,c），支持正交与三斜：
- 分数坐标 `s = (r-lo)·H⁻¹`；最小镜像把每个分量 `s -= round(s)` 后再回到笛卡尔。
- 所有距离/矢量（RDF、邻居、氢键、质心、Rg、MSD）都使用 MIC，避免跨界伪距离。
- **整分子复原（make-whole）**：当分子横跨盒面时，wrapped 坐标会把一个分子拆成两半。
  按键连接做生成树遍历，把每个原子相对其成键邻居用 MIC 摆正（等价 gmx `rmpbc`）。
  仅按“相对首原子 MIC”对长于半盒的链会失败，因此必须用键连接。`gyrate`/`rms` 依赖此法。
- **键表按 atom id 解释，且必须与拓扑一起重排**：`Bonds` 段写的是 **atom id**（不是 Atoms 段的
  行号），而 `Topology` 的所有逐原子数组都按 **id 升序**存放。data 文件的 Atoms 段**不保证**按 id
  排序（OVITO 写的文件就不是），因此解析时必须把键索引施加与坐标相同的置换。
  漏掉这一步会让键图变成随机重连，`make_whole` 去连接相距 8–45 Å 的非键合原子对，
  gyrate / rms / hbond **静默出错**（实测 OVITO 写的 PET：6476 条“键”里只有 880 条是真的，
  Rg 高估 3.6×、RMSF 高估 4.5×）。回归：`tests/unit/test_bond_atom_order.py`。

### 0.3 周期几何的唯一入口（`lat.core.periodic`，0.2.0）
0.1.0 里"周期几何"有四套并行实现，其中三处把 wrapped 坐标当成物理坐标用（SASA 丢邻居、
`contact/residence -mode center` 直接对 wrapped 求平均、`pistack` 用裸坐标算环心），
只有 `Box.minimum_image` / `pairs_within` / `make_whole` 三处是原则性的。

0.2.0 起，所有周期判断必须走一个入口，共四层：

| 层 | 函数 | 用途 |
|---|---|---|
| 分子 | `whole_positions` / `residue_com(s)` | 整分子复原、质心（可限定 `seeds` 只遍历选区所在分子）|
| 原子对 | `periodic_neighbors` / `periodic_nearest` / `periodic_vector` / `periodic_distance` | 周期邻域、周期最近点、MIC 位移/距离 |
| 体素网格 | `grid_axes` / `voxel_occupancy` / `periodic_distance_transform` / `label_periodic` | 周期占据、**环面**距离变换、**环面**连通分量（仅正交盒）|
| 曲面 | `sasa_shrake_rupley` | 周期 Shrake–Rupley |

架构守卫测试 `tests/unit/test_no_ad_hoc_pbc.py` 会失败于任何模块自建 `cKDTree` /
`distance_transform_edt` / `ndimage.label` / `np.pad` / 裸 `minimum_image`，
并强制 PBC 敏感模块按下表引用对应符号（`rdf -intramol` 除外：它统计的是**成对** MIC 距离，
在 wrapped 坐标上本来就正确，`tests/regression/level2_format/test_pbc_invariance.py` 固定了这一点）。

判据：把整个体系做一次非格点平移，所有**物理量必须不变**。
`test_pbc_invariance.py` 对 contact/min·center、residence/min·center、pistack、agg、cluster、rdf 逐列对拍两个平移位置。

### 0.4 邻居搜索
在 cutoff 内找原子对：
- 正交盒：`scipy.spatial.cKDTree(boxsize=L)`，坐标先平移到 `[0,L)`（注意盒子原点非零）。
- 三斜盒：把 B 复制到相邻镜像格点后再建树（镜像数由面高 `V/area` 决定）。
- 自身配对（`same`）只返回无序对一次；不同集合返回所有对。

### 0.5 归一化自检
凡“分布 + 积分”成对的量（RDF/配位数、密度/总质量、电流自相关/电导），
实现内部保证两者数学自洽，并在偏差超阈值时 `WARNING`（不静默修正）。

### 0.6 输出契约
每个模块返回若干 `Table`，写为 CSV。文件名 `<基名>.<标签>[.<模式>].out.csv`；
文件头含 units、selection、frame_range、scalars、warnings、command、version。
屏幕只出进度/警告/路径。物理参数**全部**来自配置或命令行，代码不含隐藏默认值。

### 0.7 符号与缩写（全文通用）

| 记号 | 含义 |
|---|---|
| `r_i` | 第 i 个原子的笛卡尔坐标（Å） |
| `lo` | 盒子的下角坐标；`L` 为盒子三边长 |
| `cutoff` | 距离截断：只统计距离小于它的原子对 |
| MIC | 最小镜像约定（minimum image），处理周期性 |
| COM | 质心（center of mass） |
| `N_A` / `N_B` | 选择 A / B 的**原子数**；**阿伏伽德罗常数记作 `N_Avo`**（避免混淆） |
| `V_f`、`V` | 第 f 帧的盒子体积 / 平均体积 |
| `shell` | 球壳体积 `(4/3)π(r_hi³-r_lo³)` |
| `N_frames` | 参与统计的帧数 |
| `step`（步长） | `frame_range` 第三段或 `-stride`，即每隔多少帧取一帧 |
| `dt_eff` | 相邻**被选帧**的时间间隔 = `frame_interval × step`（ps） |
| vdW | 范德华半径；探针半径 `r_probe` 另计 |
| FFV | 自由体积分数（free volume fraction） |
| PSD | 孔径分布（pore size distribution） |
| SASA | 溶剂可及表面积（solvent accessible surface area） |
| SSIP / CIP / AGG | 溶剂分离离子对 / 接触离子对 / 聚集体 |
| NE / GK | Nernst–Einstein / Green–Kubo |
| R² | 线性拟合的决定系数（1=完全线性） |
| `image flags`、`xu` | LAMMPS 的镜像计数与“已解包坐标” |
| 6/26 邻接 | 三维连通分量用面邻接(6)或面+棱+角邻接(26) |

---

## 1. RDF 与配位数（rdf）

**原理**：$g_{AB}(r)=\dfrac{1}{N_A}\sum_{i\in A}\sum_{j\in B}\dfrac{\langle\delta(r-|r_i-r_j|)\rangle}{4\pi r^2\Delta r\,\rho_B}$，
配位数 $n_{AB}(r)=4\pi\rho_B\int_0^r g\,r'^2dr'$（$\rho_B=N_B/V$）。

**实现**：逐帧近邻统计距离直方图；精确球壳体积 `shell=(4/3)π(r_hi³-r_lo³)`；
`g = f·Σ(counts_f·V_f) / (N_frames·N_A·N_B·shell)`，`cn = f·cumsum(Σ counts_f)/(N_frames·N_A)`；
`f=2` 当 sel1≡sel2（无序对还原有序），否则 1。
`V_f` 是第 f 帧的盒子体积（逐帧密度平均）。在**恒定体积**下二者严格满足
`cn=∫4πρ_B g r²dr`（`ρ_B=N_B/V`）；恒压（盒子涨落）时 g 按逐帧 V 加权、cn 不含体积因子，
二者是同一口径下的近似自洽，归一化自检只用于发现实现错误。

**首壳标量**：首峰=第一个 `g>1` 的局部极大；首壳极小=其后第一个 `g<0.5` 的局部极小；
找不到就省略标量并 WARN（避免噪声伪极小）。

**陷阱**：$r>$ 半盒长时 g 失真（cutoff 应 < 半盒，越界会 WARN）；多组分用 B 的自身数密度；
**同分子对默认排除**：默认**不计算**同一分子内的原子对（需要输入的 `mol` 列或 GROMACS 残基编号）；
`-intramol` 打开后**包含**分子内对。若拓扑有键但无分子 ID，无法排除会 WARN，此时须显式加 `-intramol`。
**判读用配位数而非 max|Δg|**（见第 16 章）。

---

## 2. 密度剖面（density）

**原理（gmx density 口径）**：**逐帧按分数坐标分 bin**，固定 `nbins`，最后把轴按**平均盒长**标定。
即第 f 帧内 `s = (x − lo)/L`（wrap 到 [0,1)），落进固定的 `[0,1]` 分箱；该分数薄层的体积恰为
`V_f/nbins`，所以 `rho_bin = Σ_f (hist_bin,f · nbins / V_f) / N_frames`。这样 NPT 下 bin 随盒缩放，
原子不会掉出固定边界；轴坐标 = `mean(lo) + 分数中心 · mean(L)`。2D 同理，单元格体积 `V_f/(nbins²)`。

- **盒变 WARN**：任一相关轴的帧盒长相对其均值偏离 > **±2%** 时 WARN，提示已用分数分箱、轴按平均盒长标定。
- 质量密度 `kg/m³ = (g/mol/Å³)·1e27/N_Avo`。**非正交/畸变盒的截面近似由用户自行负责**（照抄 gmx 口径）。
- **平滑不在这里做**：计算命令永远输出原始结果；`-smooth` 只在 `plot1d`/`plot2d` 生效。

**陷阱**：界面/尖峰处逐点差大但均值才反映系统偏差；`-quantity number` 用于无法得到质量的体系。

---

## 3. MSD 与扩散（msd）

**原理**：$\mathrm{MSD}(\tau)=\langle|\mathbf r(t+\tau)-\mathbf r(t)|^2\rangle$，
Einstein 关系 $D=\text{slope}/6$（3D）或 `/4`（2D）。

**实现**：解包坐标优先 `xu/yu/zu`，其次 `x + i·L`，否则 wrapped（WARN）。
GROMACS 轨迹无 image flags，读取时逐帧 MIC 重建连续坐标（与 gmx nojump 一致）。
全时间原点平均：`MSD(τ)=(1/(T-τ))Σ_t (1/M)Σ_i |r_i(t+τ)-r_i(t)|²`，并对 τ>0 全部 lag 取点。
`-seltype mol_com` 用分子质心；默认去除所选组 COM 漂移（`-noremove_com` 关闭）。
`D` 由 `fit_range` 内 lag（默认 `[max(1,0.1T), max(start+2,0.5T))`）对 MSD 做**最小二乘直线**拟合，
斜率 `/6`（3D）或 `/4`（2D）；给出 R²，<0.8 时 WARN。`D` 默认单位 `Å²/frame`；配置 `frame_interval`
后另给 `D_cm2_per_s = D_Å²/frame ÷ dt_eff[ps] × 1e-4`。**≥3 帧**。

**陷阱**：帧间隔极大时单步位移接近半盒，PBC 解包本身有歧义，任何实现都会有几 % 偏差；
无 image flags 且非 GROMACS 时结果可能饱和。

---

## 4. 回转半径（gyrate）

$R_g=\sqrt{\sum m_i|\mathbf r_i-\mathbf r_{cm}|^2/\sum m_i}$。
**关键**：先按键连接整分子复原，再算 $R_g$；否则跨盒聚合物会被系统性低估
（实测 nafion 未复原 3.7 nm vs gmx 6.8 nm，修复后偏差 <1e-3 nm）。无质量时用等权并 WARN。
键表本身必须正确解释（见 §0.2 的 id 重排）；键图错了这一步会给出比“未复原”更离谱的**高估**。
与 `gmx gyrate` 的逐值对拍见 `tests/regression/level1_algo/README.md`（Rg 逐帧一致到 3e-4 Å）。

**多分子选择（`-sel all`）的固有歧义**：Rg 不是周期体系的平移不变量。每个分子可以被放在任意
周期像里，`make_whole` 只能把“某个分子内部”摆正，分子之间的相对位置仍依赖锚点约定。
实测 4 条链、盒长 80 Å：LAT 的锚点给 27.4 Å，改成“按各链质心回绕”给 32.0 Å（差 17%）。
**结论：周期盒里报单链（或按分子分别算）的 Rg，不要用 `-sel all` 的均值下结论。**

---

## 5. RMSD / RMSF（rms）

**原理**：Kabsch 最优旋转对齐：$R=U\,\mathrm{diag}(1,1,\mathrm{sign}(\det(UV^\top)))\,V^\top$；
对齐后逐帧求 RMSD，逐原子求涨落 RMSF。

**实现**：`-fit` 指定对齐原子（默认 `-sel`）；`-ref first|average|<帧号>`；先整分子复原再做平移/旋转。
参考 `average` 是“逐帧坐标的算术平均”（`-ref average` 时 RMSD 把每帧对齐到该平均结构）。
RMSD 在主选组上算；RMSF 是**对齐后各帧坐标相对其时间平均位置的涨落**（逐原子
`sqrt(⟨|x_i-⟨x_i⟩_t|²⟩_t)`）。**≥2 帧**。

---

## 6. 接触（contact）

`-mode min`：sel1–sel2 中距离 < cutoff 的原子对数/帧；
`-mode center`：两组分别取每分子 COM 后的最近距离/帧。用 MIC。

---

## 7. 氢键（hbond）

**判据（对齐 gmx）**：D···A ≤ 3.5 Å（`-distance_mode DA`，或 `HA` 用 H···A），
且**在给体 D 处**的 H–D···A 夹角 ≤ 30°（即 $\cos\angle(D\to H,D\to A)\ge\cos30°$）。
所有矢量 MIC。需要键信息识别 D–H；`sel1`/`sel2` 必须相同或完全不重叠。

**给体必须是电负性重原子**（N/O/F/S/Cl/Br/I，与受体同一集合）。**C–H 不是氢键给体**：
把“任何带 H 的重原子”当给体，会让聚合物里每个 CH₂/CH₃ 碳都成为给体，计数比 gmx 高一个量级
（PET 实测 41.7 vs 4.3 每帧；修正后逐帧完全相同）。无电负性给体时**直接 ERROR**，
与 gmx 的 `No Donors found` 行为一致（`sys1_reduced` 这类只有 C–H 的体系即属此列）。

**统计**：逐帧计数；存在函数自相关
$C(\tau)=\langle h(t)h(t+\tau)\rangle/\langle h(t)\rangle$，分间断（可离开再回）与连续（不允许中断）；
寿命由两路给出（单位=帧）：**积分法** `τ_int=Σ_{τ≥1} C(τ)`；**指数拟合法** 对 `τ≥1` 的
`ln C(τ)` 做直线拟合，`τ_fit=-1/斜率`，仅当 `0<τ_fit≤T` 且 `R²` 尚可时报告，否则记 `nan` 并 WARN
（拟合窗默认从 τ=1 到曲线末端）。**自相关必须对有效起点无偏**（除以 $T-\tau$），
否则永久存在的键会出现 $1-\tau/T$ 的伪衰减。

**实现要点**：对每个 D–H，先用周期性邻居搜索找 cutoff 内的候选受体（O(N·邻居)），
再算角度，避免“每个 D–H 扫全部受体”的 O(N²)。几何分布取 D–A 中点投影到 xz/yz/xy；
中点为 **A 与“离 A 最近的 D 的周期像”的中点**，即 `mid = A + 0.5·MIC(D−A)`，并 wrap 回盒内
（跨边界不丢点）。**注意不要写成 `0.5·(MIC(D−A)+A)`**——那会把中点压到半个盒子里。

---

## 8. 聚合物种 SSIP/CIP/AGG 与团簇（agg / cluster）

**本质**：邻接 → 图 → 连通分量 → 按端点计数分类。
1. 对 sel1（中心，如 Li⁺）与 sel2（反离子，如 FSI⁻）按 `-mode min|center`、`-cutoff` 建邻接；
2. 建**残基级二部图**（一个多齿离子计一条边）；
3. 求连通分量，记录每个分量 $(i,j)$ = A 数、B 数；
4. 分类（对 A）：$A_1B_0\to$SSIP，$A_1B_1\to$CIP，其余 $\to$AGG。

**为什么不用贪心最近邻**：一对一“掐尖”会漏掉桥接离子，系统性低估 AGG；必须建真图。
`cluster` 共用引擎，`-focus` 时退化为“中心离子周围 B 的配位数分布”。

**陷阱**：1:1 电解质两侧平均配位数相等但 SSIP/CIP/AGG 分数不对称，输出必须声明中心与被计数方；
cutoff 过大导致全连通时会 WARN。

---

## 9. 滞留时间（residence）

参考组 `-ref` 与配体组 `-lig`：配体到参考距离 < cutoff 记 $h=1$。
间断/连续自相关同氢键；寿命同样给积分（`τ_int=Σ_{τ≥1}C`）与指数拟合（对 `τ≥1` 的 `ln C` 拟合，
`τ_fit=-1/斜率`）两路，单位=帧；另给每个配体连续 $h=1$ 的区间长度分布。
配体按分子分组；判据半开 `d<cutoff`。

---

## 10. 输运系数（transport）

**前置**：轨迹必须含 `vx/vy/vz`（否则报错）。
- **Einstein/Nernst–Einstein**：各物种 $D_i$ 由 MSD 得，
  $\sigma_{NE}=\frac{F^2}{RT}\sum z_i^2 c_i D_i$，迁移数 $t_i=\frac{z_i^2c_iD_i}{\sum_j z_j^2c_jD_j}$。
- **Green–Kubo**：电流 $\mathbf J(t)=\sum q_i\mathbf v_i$，
  $\sigma_{GK}=\frac{1}{3Vk_BT}\int_0^\infty\langle\mathbf J(0)\cdot\mathbf J(\tau)\rangle d\tau$，输出运行积分曲线。

**单位链（可手算复核）**：`step` 见 §0.6；`dt_eff = frame_interval × step`（ps）。
`D[m²/s] = D[Å²/frame] ÷ dt_eff[ps] × 1e-10`（因 1 Å²=1e-20 m²、1 ps=1e-12 s）。
`D[cm²/s] = D[m²/s]×1e4`。速度用 `-velocity_unit`（`A/ps`/`A/fs`）换算到 m/s 才能算电流。
电荷 `q_i` 单位为元电荷 e，电流 `J=Σ q_i v_i`。
**worked example**：`frame_interval=1 ps`、`step=10` → `dt_eff=10 ps`；若 `D=0.5 Å²/frame`，则
`D=0.5/10×1e-8 m²/s=5×10⁻¹⁰ m²/s=5×10⁻⁶ cm²/s`。
缺 `frame_interval` 只出 D（`Å²/frame`），NE/GK 跳过并 WARN。NE 忽略离子关联，通常高估电导，需声明。

---

## 11. 区域占据与穿越（occupancy）

区域语法 `轴:min-max[,min-max]`，半开区间 `[min,max)`，多轴笛卡尔积。
每帧给每个原子/分子打区域号（0=区域外），派生：区域占据数、穿越计数、首次通过时间、富集
（区域物种密度/体系平均）。区域体积按区间与盒子交集计。

---

## 12. 自由体积 / 孔径 / 连通性 / SASA（`freevol` + `sasa`）

0.2.0 把原来的单模块 `probe` 拆成两个，因为它们是**两类不同算法**，只共享"vdW 半径 + 探针半径"：

- **`lat freevol`**（体素网格类）：FFV / PSD / 连通性。占用 = 格点中心落在任一原子的
  $(r_{vdW}+r_{probe})$ 球内，各轴按 MIC 回绕。
- **`lat sasa`**（曲面类）：Shrake–Rupley，不需要网格。

拆的理由：`probe -mode all` 把 3 个体素网格模式和 1 个曲面模式绑在一起跑；PSD/连通性很少用，
而 PSD 每帧还要额外算一次距离变换（实测 PET 体系 6.6 s/帧 vs FFV 0.5 s/帧）。
`probe` 已删除，**不留兼容壳**：传入 `probe` 会直接报"unknown module"。

### 12.1 freevol

**共同前提**：网格分析要求**正交且固定**的盒子。笛卡尔体素网格只在正交盒上周期自洽；
三斜盒直接报错；所选帧内盒子漂移 >2% 报错，0.1%–2% 给 WARN 且网格锁定在首帧。
每轴 `n=round(L/dr)`、格点取 bin 中心、`-grid_dr` 默认 0.5 Å；总格点超过上限
（ffv 5e7、psd/connectivity 2.5e7）时按 $(n_{total}/cap)^{1/3}$ 放大 `dr` 并 WARN（因此实际 `dr` 可能大于请求值）。

- **ffv**：$FFV=$ 空闲格点/总格点，逐帧输出 `ffv_mean`/`ffv_std`。占用判据与 0.1.0 相同，数值可逐位复现。
- **connectivity**：自由格点的连通分量，**环面闭合**，默认 **26 邻接**（6 连通固体的补集本该用 26）。
  输出逐帧汇总 `n_components / n_percolating / v_largest / frac_largest / v_percolating`，
  `-components` 才另写逐分量表（0.1.0 无条件逐帧写全部分量）。
- **psd**：见 12.2。

### 12.2 PSD 的定义与归一化（0.2.0 重定义）

1. 距离变换对**未加探针的 vdW 占据**做（Gelb–Gubbins 口径），得每个格点到最近固体表面的距离 $d$；
2. 可及性用**加了探针**的占据集：$free=\lnot occ(r_{vdW}+r_{probe})$；
3. 在 free 区取 $d$ 的局部极大为孔心：用 26 连通**平台去重**（立方格点的 EDT 到处是山脊/平台，
   直接用 `==` 会把整片平台都当孔心——0.1.0 因此报出 0.055 Å 的伪主峰）；
4. 每个自由格点按**周期最近孔心**归属（Voronoi 分配），权重 = 该孔心分到的自由体积；
5. $dV/dr$ 用固定 bin 宽（`= grid_dr`）与范围 `[0, -psd_rmax]`（默认由数据定），
   归一化到 $\int (dV/dr)\,dr = 1$；文件头 `scalars: integral` 给出实测积分。

与 0.1.0 的差别（每条都在文件头/本节声明）：0.1.0 用 $r=d-r_{probe}$（多减一个探针半径）、
bin 宽随数据范围变化（PET 实测 0.944 Å，而文件头写 `grid_dr=0.5`）、只统计孔心格点
（$\int = 1.5\times10^{-4}$，即只有 0.015% 的自由体积被计入）、且只用第 0 帧却在文件头写全部帧。
0.2.0 逐帧累加（PET 回归用 `-stride 10` 覆盖 0–500 帧），文件头的 `frame_range` 与实际数据一致。

### 12.3 connectivity 的可读性

0.1.0 在 FFV≈83% 的体系里逐帧输出 ~165 个分量（501 帧 = 8 万行），其中 112 个是单体素栅格粉尘，
最大分量恒为盒体积的 83%、恒标 `percolating=1`——表里看不出"结构是否贯通"。
0.2.0 逐帧汇总后 `frac_largest` 直接回答这个问题：PET 体系 frame 0 为 0.99989、frame 1 为 0.99998，
即自由空间基本是**一个贯通网络 + 少量孤立小孔**（分量数 35 → 17 是这些小孔的消长，不是数值噪声）。

### 12.4 sasa

Shrake–Rupley：每个原子以 $r_{vdW}+r_{probe}$ 为半径撒 M 个 Fibonacci 球面点（默认 92），
落在任何邻居球内即被遮挡，$SASA = 4\pi r^2 \times$ 暴露点比例。

- **邻域搜索与遮挡判定都用 MIC**。0.1.0 用裸 `cKDTree(pos)`，贴面原子丢周期性邻居 → 系统性高估：
  PET 体系对 `gmx sasa` 偏高 **+14.2%**，修正后 **−0.2%…−0.8%**（§12.5）。
- 表面**始终用全部原子**构建，`-sel` 只决定"单独报哪部分原子的贡献"；输出两列
  `sasa_total` 与 `sasa_sel`。0.1.0 只用被选中的原子做遮挡，会高估选区面积。
- 迁移：`lat probe -mode sasa -sel X` → `lat sasa -sel X`。

### 12.5 与原生 gmx 对拍（2026-10，已改行为）

`tests/regression/level1_algo/gmx_ref_lammps.py` 把 LAMMPS 输入转成 GROMACS 后与原生 gmx 逐帧对拍。
PET 柱状体系（6200 原子 × 501 帧，`gmx 2022.5`）：

| 量 | lat | gmx | 平均\|Δ\| | max\|Δ\| | 相对 |
|---|---|---|---|---|---|
| gyrate（全体）| 26.8319 Å | 26.8316 Å | 0.0003 Å | 0.0004 | 0.001% |
| gyrate（单链 1–4）| 25.2776 / 25.3872 / 24.0370 / 25.4916 Å | 25.2778 / 25.3874 / 24.0371 / 25.4918 | 0.0002 Å | 0.0005 | 0.001% |
| RMSD（`gmx rms -nomw`）| 8.2829 Å | 8.2830 Å | 0.0001 Å | 0.0002 | 0.001% |
| RMSF（`gmx rmsf -fit yes`，6200 原子）| 2.3167 Å | 2.3167 Å | 0.0007 Å | 0.0038 | 0.032% |
| **SASA（0.1.0）** | **18316.9 Å²** | **16033.4 Å²** | **2283.6 Å²** | **2917.3** | **+14.24%** |
| hbond（给体数）| 6.6028 | 6.6168 | 0.0140 | 1.0 | 0.21% |
| RDF 配位数（C-O / C-C / O-O）| 6.2959 / 6.3286 / 4.9969 | 6.2818 / 6.3284 / 4.9968 | 0.09 / 0.10 / 0.04 | — | 0.26 / 0.27 / 0.27% |

> SASA 一行是 **0.1.0** 的实测值（说明缺陷量级）。0.2.0 的周期化实现复测（同一 501 帧窗口）：
> **全 501 帧均值 15946.5 Å² vs `gmx sasa` 16033.4 Å²，平均绝对偏差 94.1 Å² = +0.59%**
> （逐帧 −0.17%…−0.85%）；残余偏差来自球面采样点分布（gmx 与 lat 各自用 92 点但分布不同），
> 详见 `docs/11_examples_recheck.md`。
> 大体系复测：32800 原子的水/表面活性剂体系，`lat hbond -sel1 all -sel2 all` 均值
> 13908.8 vs `gmx hbond`（同 40 帧窗口）13808.5，差 **+0.7%**。

**gmx 的给体判据（与 lat 0.2.0 一致）**：只有电负性重原子（N/O/F/S/Cl/Br/I）上连的 H 才算给体。
实测（PET 体系，`gmx hbond` 输出原文）：

- `El_C`(给体组) / `El_O` → `Found 0 donors and 1128 acceptors` + `No Donors found`（**直接报错**）；
- `El_O` / `El_C` → 同样 0 donors（因为 H 不在两组之并集里）；
- `System` / `System` → `Found 8 donors and 1128 acceptors`（正是 4 条 PET 链的 8 个端羟基）。

**gmx 的选择语义（与 lat 不同，见 §17）**：`gmx hbond` 在**两组之并集**内识别给体/受体，
统计的是两组"之间"的氢键（两个方向都算）；lat 的 `-sel1` 是**给体组**、`-sel2` 是**受体组**，
方向明确。`-sel1 all -sel2 all` 时两者等价。

**受体集合的差异（2026-10 新增，已确认未改）**：lat 的受体是 `{N,O,F,S,Cl,Br,I}`，
gmx 在这三个体系里只把 **O** 算作受体（PET 1128、22 号 8480、24 号 4136 个受体原子）：

| 体系 | lat − gmx | 原因 |
|---|---|---|
| PET（只有 C/H/O）| +0.21% | 集合相同 |
| 水/表面活性剂（有 S）| +0.7% | 硫酸根 S 很少满足几何判据 |
| LiCl/电极（144 个 Cl⁻）| **+13.8%** | Cl⁻ 是强受体，gmx 不认 |

即：**含卤素（或硫）的体系里，lat 的氢键计数会高于 gmx，且差值可解释**。
lat 的口径在化学上更完整（Cl⁻/S 确实是受体）；要逐值对齐 gmx 需要一个 `-acceptors gmx` 开关，
本轮未加（见 `docs/12_0.2.0_release.md` §6）。

---

## 13. π–π 堆积（pistack）

每个残基的选原子定义一“环”：质心 + 由 SVD 得到的法向。
两环“堆积”判据：质心距离 < `-distance`（默认 5.5 Å）且法向夹角 < `-angle`（默认 30°）。
输出每帧计数与逐对几何（距离、夹角）。

---

## 14. 区域/选择工具与绘图

- `select`：显示选择命中数、类型/元素分布、电荷和、分子数、残基名；可导出 `.ndx`。
- `describe`：列模块与参数（参数、默认来源、说明）。
- `plot1d`/`plot2d`：只读 CSV 画图，与计算解耦。

---

## 15. 与 GROMACS 的一致性（实现口径对照）

- RDF：gmx 逐帧 `N_sel/V` 平均、精确球壳；lat 同样精确球壳、逐帧密度。
  唯一系统性差异是 **abscissa 约定**：gmx 报告 bin 左边缘 `i·Δr`，lat 报告 bin 中心 `(i+½)Δr`。
- MSD：gmx 默认不去 COM；lat 默认去所选组 COM（`-noremove_com` 对齐 gmx）。
- 氢键：已对齐 gmx 的给体处角度判据与 3.5 Å。
- 密度：两者均按帧内盒的垂直截面积归一。

---

## 16. 精度、误差来源与判读（重要）

“照抄公式”仍会有可观差异，来源按重要性排序：

1. **有限采样（主导）**。g(r) 是系综平均，轨迹长度/粒子数有限时逐点有涨落。
   例：spdf 体系某原子对最大 |Δg|=0.045，而该 r 处**逐帧标准误**为 0.025（≈1.8σ）；
   同时配位数（积分）相对误差仅 **0.02%**。因此逐点最大值不是好指标。
2. **网格/栅格约定**。bin 中心 vs 左边缘、bin 宽度、精确球壳 vs 近似，在**尖锐峰/界面**处
   把小偏差放大成大的 max|Δg| 或 max|Δρ|。实测：spdf 的 MO1–OH1 在最大 |Δg|=0.045 处，把 gmx
   横坐标平移 ±1 pm（半 bin）后 max|Δg| 仍在 0.045–0.064，说明该点**不是**半 bin 平移主导，而是
   噪声与网格叠加；只有界面/峰足够陡时网格约定才成为主因（密度剖面 64→4037 kg/m³ 即属此类）。
3. **帧集/网格对齐**。不同工具用不同窗口或抽帧方式时差异巨大；
   对拍必须使用**同一组帧**（本仓库回归脚本已强制 stride=1 同窗口）。
4. **PBC 解包**。极短轨迹 + 极大步长时位移接近半盒，解包有歧义（MSD 可差 3–5%）；
   跨盒分子必须整分子复原（曾在 `gyrate` 造成 >40% 的真实错误，已修）。
5. **归一化实现细节**。逐帧密度平均 vs 平均密度（NPT 下 ~1e-4）；积分端点/半权。
6. **对拍工具自身问题**。如 gmx 在极小 MSD 上数值下溢（~1e-41）；本机 gmx 2023.2 无法读取
   RefData 中 VERSION 2025.1 的 tpr（选择为空并报内部错误），而 lat/MDAnalysis 可读。
   对拍前须确认参考工具能真正处理同一输入。

> 说明：上面 1–6 是按**在本仓库实测样本上的出现频率与量级**的经验排序，其中 1、4、6 有具体数字
> 支撑，2、3、5 为定性判断（在不同体系上可能互换主次）。

**结论/判读规则**：
- 用**积分量**判断正确性：配位数、平均密度、Rg 均值、D（斜率）。
- “积分量差异 >1–2% 视为缺陷”是**判读约定**（依可见系统偏差经验选定），不是物理常数；
  逐点最大差异单独列出但不单独定性。
- 报告差异时必须同时给出：帧数、粒子数、bin 宽、窗口、以及参考工具版本。

---

## 17. 适用范围与未实现

- **适用**：正交/三斜周期盒；LAMMPS `data`+`dump`；GROMACS `tpr`+`xtc/trr`（内部回归用）。
- **不适用/限制**：
  - 粗粒化、无元素体系不支持 freevol/sasa/hbond（化学身份缺失，报错清晰）。
  - 反应轨迹原子数逐帧变化：报错并建议限定一致窗口。
  - `hbond` 存在的自相关在数千帧长轨迹上为 O(T²)，建议窗口/子采样。
  - **`freevol` 的体素网格要求正交且固定盒**（三斜/盒子漂移 >2% 直接报错，0.1%–2% WARN）；
    PSD/connectivity 都是网格近似，不是 Zeo++ 的 Voronoi/几何算法。
  - **`contact` / `residence -mode center` 的选择若"跨盒且非整分子"**，质心依赖参考原子，
    lat 会给 WARN 并报出歧义大小（阈值 0.1 Å）；要精确请选整分子。gmx 在等价场景下
    靠 `rmpbc` + 质心（质量加权）实现，见 `docs/11_examples_recheck.md` §4。
  - **`hbond` 的选择语义与 gmx 不同**：lat 的 `-sel1` 是给体组、`-sel2` 是受体组（有方向）；
    `gmx hbond` 在两组的**并集**里识别给体/受体，统计两组之间的氢键（方向无关）。
    两者只在 `-sel1 all -sel2 all` 时等价。当 `-sel1` 里没有电负性 H 时 lat 报错
    （gmx 在并集也没有 H 时同样报 `No Donors found`）。
  - `-mode center` 的质心是**几何中心**（未按质量加权），gmx 的 `com of group` 是质量加权，
    实测差异 0.01 Å（PET 体系，见 `docs/11_examples_recheck.md` §4）。
  - **`gyrate -sel all` 在多分子周期体系里不是良定义量**（周期像歧义可达 17%，§4）。
  - `occupancy` 的“按分子体积自动生成 zone”未实现。
  - 能量/自由能、PMF、对接、phonon、AIMD 等不在范围内（见 `docs/01_meta_requirements.md`）。
