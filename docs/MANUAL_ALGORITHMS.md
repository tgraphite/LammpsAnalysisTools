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

### 0.3 邻居搜索
在 cutoff 内找原子对：
- 正交盒：`scipy.spatial.cKDTree(boxsize=L)`，坐标先平移到 `[0,L)`（注意盒子原点非零）。
- 三斜盒：把 B 复制到相邻镜像格点后再建树（镜像数由面高 `V/area` 决定）。
- 自身配对（`same`）只返回无序对一次；不同集合返回所有对。

### 0.4 归一化自检
凡“分布 + 积分”成对的量（RDF/配位数、密度/总质量、电流自相关/电导），
实现内部保证两者数学自洽，并在偏差超阈值时 `WARNING`（不静默修正）。

### 0.5 输出契约
每个模块返回若干 `Table`，写为 CSV。文件名 `<基名>.<标签>[.<模式>].out.csv`；
文件头含 units、selection、frame_range、scalars、warnings、command、version。
屏幕只出进度/警告/路径。物理参数**全部**来自配置或命令行，代码不含隐藏默认值。

### 0.6 符号与缩写（全文通用）

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

## 12. 自由体积 / 孔隙 / SASA（probe）

共用“探针 + vdW 球 + 空间占据”内核，`-mode` 切换：
- **ffv**：盒内布网格，格点若落在任一原子 $(r_{vdW}+r_{probe})$ 球内则算占据，
  FFV=空闲/总。标记时各轴周期回绕（贴边原子不漏）。
- **psd**：对空闲空间做距离变换 $d(x)$，孔半径 $r=d-r_{probe}$，局部极大为孔心，按体素体积加权成 PSD
  （探针-网格近似，与 Zeo++ 的 Voronoi 不同，不做外接）。
- **porosity**：空闲体积/总体积。
- **connectivity**：空闲格点连通分量（6/26 邻接），输出体积与是否贯通。
- **sasa**：Shrake–Rupley，按 $(r_{vdW}+r_{probe})$ 球面采样 M 点，暴露点比例×球面积。

**可复现参数**：网格每轴 `n=round(L/dr)`、格点取 bin 中心，`-grid_dr` 默认 0.5 Å；当总格点 >5e7 时
自动按 `(ntotal/5e7)^(1/3)` 放大 `dr` 并 WARN（因此实际 `dr` 可能大于请求值）。PSD 用
`distance_transform_edt`，以 `r=d-r_probe`，取 3×3 局部极大为孔心，按体素体积 `dr³` 加权成 `dv/dr`。
SASA 用 Shrake–Rupley，球面采样点数 `M`=`config.modules.probe.sasa_points`（默认 92）。
vdW 半径来自配置元素表（全周期表）；探针默认 1.445 Å（H₂）。缺元素会报错并提示补
`Masses`/`element`/`-types`。connectivity 当前**不闭合 PBC 环面**（按首帧报告）。

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
  - 粗粒化、无元素体系不支持 probe/hbond（化学身份缺失，报错清晰）。
  - 反应轨迹原子数逐帧变化：报错并建议限定一致窗口。
  - `hbond` 存在的自相关在数千帧长轨迹上为 O(T²)，建议窗口/子采样。
  - `probe` 的 PSD/connectivity 为网格近似；connectivity 暂不闭合 PBC 环面。
  - `occupancy` 的“按分子体积自动生成 zone”未实现。
  - 能量/自由能、PMF、对接、phonon、AIMD 等不在范围内（见 `docs/01_meta_requirements.md`）。
