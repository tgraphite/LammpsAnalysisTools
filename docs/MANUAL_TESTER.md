# lat 测试人员手册（v0.1）

本手册面向**测试人员**：如何在新设备上安装、如何运行每个功能、如何判断结果是否正常、
如何选择测试输入、以及发现问题时如何反馈。读完本手册即可独立完成一轮验收测试。

> 术语：**拓扑** = 原子类型/质量/键（LAMMPS `*.data` 或 GROMACS `*.tpr`）；
> **轨迹** = 逐帧坐标（LAMMPS `*.lammpstrj` 或 GROMACS `*.xtc/.trr`）；
> **选择** = 用一句话描述“算哪些原子”（见第 4 节）。

---

## 1. 安装（一次即可，无需每次 source）

在任意新设备上：

```bash
# 拿到交付包里的 lat/ 目录后（或已有的 /opt/lat）
cd lat
bash install.sh
export PATH="$HOME/.local/bin:$PATH"   # 若安装脚本提示不在 PATH
lat --version                          # 应输出 lat 0.1.0
```

若机器自带 Python < 3.9（如 CentOS 7/8 的 3.6），`install.sh` 会自动用 `uv`
下载一个私有 Python 3.11，无需 root、无需手动装 Python。

`install.sh` 会：建立独立虚拟环境 `.venv`、安装全部依赖、在 `~/.local/bin/lat`
放置启动器。安装后**直接**用 `lat ...`，不需要 `source activate`，也不需要 `cd` 到源码目录。

可配置的环境变量：

| 变量 | 作用 | 默认 |
|---|---|---|
| `LAT_CONFIG` | 指定自定义配置文件（覆盖内置 `config/lat.yaml`） | 内置配置 |
| `LAT_VENV` | 安装时虚拟环境路径 | `<repo>/.venv` |
| `LAT_BIN` | 安装时启动器输出目录 | `~/.local/bin` |
| `PYTHON` | 安装时使用的解释器 | `python3` |

依赖：Python ≥ 3.9、`numpy scipy pyyaml matplotlib MDAnalysis ase`（`install.sh` 自动安装）。

### 通过 SSH 在新设备试装（示例别名 `comp-portal`）

```bash
scp -r /opt/lat comp-portal:/tmp/lat
ssh comp-portal 'cd /tmp/lat && bash install.sh && ~/.local/bin/lat --version'
# 冒烟：用自带样例跑 rdf
ssh comp-portal '~/.local/bin/lat rdf -s /tmp/lat/tests/regression/samples/sys1_reduced/sys1_reduced.data \
  -f /tmp/lat/tests/regression/samples/sys1_reduced/sys1_reduced.lammpstrj \
  -sel1 "type 2" -sel2 "type 5" -cutoff 8 -frame_range 0:3:1 -o /tmp/smoke'
```

---

## 2. 30 秒上手

仓库自带一个 6600 原子的样例，可直接跑：

```bash
S=/opt/lat/tests/regression/samples/sys1_reduced/sys1_reduced
lat rdf -s $S.data -f $S.lammpstrj -sel1 "type 2" -sel2 "type 5" -cutoff 8 -o /tmp/demo
lat --help                      # 顶层帮助与模块列表
lat describe -module rdf        # 模块输入/输出/参数；具体默认值看 `lat rdf -h`
```

- 通用参数：`-s`（拓扑，必填）、`-f`（轨迹）、`-o`（输出基名，缺省用拓扑文件名）、
  `-unit`、`-style`、`-types`、`-index`、`-frame_range start:stop[:stride]`、`-stride`、`-v`（显示 INFO 进度）。
- **结果全部写入 CSV 文件**；`stdout` 基本为空，进度/警告/错误路径全部走 **stderr**
  （`[lat INFO/WARNING/ERROR] ...`）。错误时返回码为 1。
- 输出文件名规则：`<基名>.<标签>.<模式>.out.csv`，例如 `out.rdf.out.csv`、`out.rdf.cn.out.csv`。
- **标量结果（配位数、D、R²、turnover 数等）在 CSV 的 `# scalars:` 头注释里**，例如
  `# scalars: r_first_peak/A=2.06; g_first_peak=11.38; cn_first_shell=0.13`。数值列在
  `# --- data ---` 之后。
- 重复运行同一 `-o` 会先把旧文件备份为 `.bak.1`、`.bak.2`…（不会覆盖丢失）。

---

## 3. 支持的输入格式

| | 拓扑 | 轨迹 |
|---|---|---|
| LAMMPS（首选） | `*.data` / `*.lmp`（含 `Masses`、`Atoms`；含 `Bonds` 可做氢键/聚合） | `*.lammpstrj`（`dump`） |
| GROMACS | `*.tpr` / `*.gro` / `*.pdb` | `*.xtc` / `*.trr` / `*.gro` |

- LAMMPS 的 `atom_style` 自动识别（`atomic`/`charge`/`molecular`/`full`，支持 `full/kk` 等加速器后缀）；
  六列歧义时才需要 `-style`。
- 类型→元素采用质量匹配（±0.2%，IUPAC 原子量）；不确定时用 `-types 文件` 显式指定，
  或让 dump 带 `element` 列。
- 轨迹坐标支持 `x/y/z`、`xu/yu/zu`、`xs/ys/zs`、`xsu/ysu/zsu`；带 `ix/iy/iz` 可解包。
- 若某帧与拓扑原子数不一致（反应轨迹中途增删原子），会报明确错误并建议用
  `-frame_range`/`-stride` 限定到一致窗口。

---

## 4. 选择（selection）速查

写法：MDAnalysis 风格表达式，可用 `and`/`or`/`not` 与括号。

| 输入 | 可用关键字 | 例子 |
|---|---|---|
| LAMMPS | `type`, `element`, `mol`, `id`, `index` | `-sel "type 2"`，`-sel1 "type 1 2" -sel2 "type 3 4"` |
| GROMACS | `resname`, `name`, `resid`, `type`, `element`, `mol` | `-sel1 "resname LI" -sel2 "resname FSI"` |

> 注意：`resname` **只适用于 GROMACS**；LAMMPS 拓扑没有残基名。
> `all` 表示全部原子。**单选择模块**（gyrate/rms/msd/residence 的 `-lig`、probe 的 `-sel` 等）
> 用 `-sel`；**双选择模块**（rdf/contact/hbond/agg/cluster）用 `-sel1`/`-sel2`，此时用 `-sel` 会
> 因参数歧义报错（`-sel` 可匹配 `-sel1/-sel2`）。

---

## 5. 模块操作手册（逐个）

下表中的 `BASE` 指 `-s` / `-f` 两个通用参数。除非注明，均支持 `-frame_range`。

### 5.1 rdf — 径向分布函数与配位数
```bash
lat rdf BASE -sel1 "type 2" -sel2 "type 5" -cutoff 10 [-nbins 200] [-intramol] -o out
```
输出 `out.rdf.out.csv`（`r, g_…`）与 `out.rdf.cn.out.csv`（`r, cn_…`）。
**怎么判断正常**：`g(r)` 在近距为 0、首峰后回到 1 附近；`cn` 单调上升到平台；
首峰位置应等于该原子对的典型键长/接触距离。`sel1` 与 `sel2` 是同一集合时会自动按有序对处理。
**同分子对默认不计算**（默认排除）；需要分子内 RDF 时加 `-intramol`。若拓扑有键但无分子 ID，会 WARN 提示无法排除。

### 5.2 density — 密度剖面
```bash
lat density BASE -axis z [-nbins 200] [-quantity mass|number] -o out
lat density BASE -dim 2 -axis z ...        # 2D 网格（长表 x,y,rho）
```
输出 `out.density.out.csv`（1D）或 `out.density.2d.out.csv`（2D）。质量密度单位 `kg/m³`。
NPT 盒变时按分数分箱、轴按平均盒长标定，盒长偏离均值 >±2% 会 WARN。
**怎么判断正常**：均匀液体应为水平线；出现明显分层/界面峰应与你对体系的预期一致。
**平滑只在绘图做**：`lat plot1d -f out.density.out.csv -smooth gaussian -o out.png`（CSV 始终是原始值）。

### 5.3 msd — 均方位移与扩散系数
```bash
lat msd BASE -sel "type 2" [-dim 3|2] [-seltype atom|mol_com] [-noremove_com] [-fit_range a:b] -o out
```
输出逐 lag 的 MSD 及标量 `D, R²`（在 `# scalars:` 头）。**需要 ≥3 帧**。
**怎么判断正常**：长轨迹上 MSD 对 lag 近似线性、`R²` 接近 1。**短轨迹（几帧）MSD 可能非单调、
R² 无意义，属正常**，不要据此判失败。`D` 默认单位 `Å²/frame`；配置了 `frame_interval` 时
额外给 `D_cm2_per_s`。`-noremove_com` 可关闭整体质心漂移扣除。
若轨迹无 `xu`/`ix`（LAMMPS）且非 GROMACS，会屏显 WARNING（wrapped 坐标下 MSD 可能饱和）。

### 5.4 gyrate — 回转半径 Rg
```bash
lat gyrate BASE -sel "resname POLYMER" -o out
```
输出 `frame, rg`。**怎么判断正常**：链状分子 Rg 在 nm 量级；跨盒大分子不会异常缩小
（内部已做“整分子复原”）。

### 5.5 rms — RMSD / RMSF
```bash
lat rms BASE -sel "resname PROTEIN" [-fit "name CA"] [-ref first|average|<帧号>] [-rmsf] -o out
```
输出 `out.rmsd.out.csv`；`-rmsf` 时另出 `out.rmsf.out.csv`。**需要 ≥2 帧**。

### 5.6 contact — 接触数/距离
```bash
lat contact BASE -sel1 "type 1" -sel2 "type 3" -cutoff 5 [-mode min|center] -o out
```
输出每帧原子对数（`min`）或质心距离（`center`）。

### 5.7 hbond — 氢键
```bash
lat hbond BASE [-sel1 "resname SOL"] [-sel2 "resname SOL"] [-distance 3.5] [-distance_mode DA|HA] [-angle 30] -o out
```
输出 `out.hbond.out.csv`（计数/帧）、`out.hbond.acf.out.csv`（间断/连续自相关与寿命）、
`out.hbond.dist_{xz,yz,xy}.out.csv`（几何分布）。**需要键信息与 ≥2 帧**。
判据：D···A ≤ 3.5 Å 且 H–D···A 夹角 ≤ 30°（在给体处，与 gmx 一致）。

### 5.8 agg / cluster — 离子聚合（SSIP/CIP/AGG）与团簇
```bash
lat agg BASE -sel1 "type 2" -sel2 "type 5" -cutoff 3.5 [-mode min|center] -o out
lat cluster BASE -sel1 "type 2" -sel2 "type 5" -cutoff 3.5 -o out
```
`agg` 输出 `out.agg.out.csv`（每帧各类团簇数）与 `out.agg.labels.out.csv`（SSIP/CIP/AGG 分数）。
`-cutoff` **必填**。需要键/分子信息。

**列的含义**：`A_iB_j` = 含 i 个中心离子 A、j 个反离子 B 的连通团簇的数量。标签按中心 A 计：
`SSIP` = 孤立 A（A₁B₀，溶剂分离离子对），`CIP` = A₁B₁（接触离子对），`AGG` = 其余（聚集体）。
**怎么判断正常**：同一物种的 `SSIP+CIP+AGG` 分数应≈1；`-cutoff` 应与 rdf 首壳极小一致；
屏显的 `charge imbalance` 是**警告不是错误**（两侧电荷不必相等，取决于选择）。

### 5.9 residence — 滞留时间
```bash
lat residence BASE -lig "type 5" -ref "type 2" -cutoff 3.5 [-mode min|center] -o out
```
输出占据数、自相关、寿命（单位=帧）与连续区间分布。

### 5.10 transport — 输运系数
```bash
lat transport BASE -sel "type 2" -temperature 300 [-velocity_unit A/ps|A/fs] -o out
```
**轨迹必须含速度列 `vx vy vz`**（否则报错）。输出扩散系数、Nernst–Einstein 电导、Green–Kubo
电导积分与迁移数。需要 `-temperature` 与配置里的 `frame_interval` 才能给出电导。

### 5.11 occupancy — 区域占据与穿越
```bash
lat occupancy BASE -sel "type 2" -zone "z:0-10,10-20 x:0-5" -o out
```
区域语法 `轴:min-max[,min-max]`，半开区间 `[min,max)`。输出占据数、穿越计数、首次通过、富集。

### 5.12 probe — 自由体积 / 孔隙 / SASA
```bash
lat probe BASE -mode ffv|psd|porosity|connectivity|sasa|all [-radius 1.445] [-grid_dr 0.5] -o out
```
输出 `out.probe.ffv.out.csv`、`out.probe.psd.out.csv`、`out.probe.connectivity.out.csv`、
`out.probe.sasa.out.csv`。`sasa` 可用 `-sel` 限定子集。需要能推断元素（否则报错并提示用 `-types`）。

### 5.13 pistack — π–π 堆积
```bash
lat pistack BASE -sel1 "type 3" [-sel2 "type 3"] [-distance 5.5] [-angle 30] -o out
```
`sel1` 每个残基应给出一个芳环（≥3 个环原子）。输出计数与逐对几何。

### 5.14 select / describe / plot1d / plot2d
```bash
lat select BASE -sel "type 2"                 # 屏显选择结果，可导出 .ndx
lat describe [-module rdf]                     # 列出模块/参数
lat plot1d -f out.rdf.out.csv -o rdf.png        # 输入用 -f（无位置参数）
lat plot2d -f out.density.2d.out.csv -o d.png
lat plot1d -f out.agg.labels.out.csv -columns A_SSIP,A_CIP,A_AGG -o spec.png   # 只画指定 y 列
lat plot1d -f out.agg.labels.out.csv -bar   -o spec_bar.png     # 分组柱状图
lat plot1d -f out.agg.labels.out.csv -stack -o spec_stack.png   # 归一化堆积（占据比）
```
`-columns` 的 x 轴恒为第一列；`-bar`/`-stack` 互斥、仅限 `plot1d`、不做平滑。

---

## 6. 正式测试流程

### 6.1 三级证据（缺一不可）
1. **单元/合成**：`pytest tests/unit` —— PBC、密度均匀性、随机游走 MSD、刚体 RMSD≈0、
   手算氢键/SSIP-CIP-AGG/ππ、区域穿越、FFV 等。
2. **端到端冒烟**：用自带样例 `tests/regression/samples/sys1_reduced/` 逐个模块跑通。
3. **数值对拍**：`tests/regression/level1_algo/compare_gmx.py`（RDF/CN 对 gmx，CN 误差应 ≤0.2%）。

一键全量：`pytest tests`（当前 64 项，约 1 分钟）。

### 6.2 真实订单回归（scripts）
- GROMACS 对 gmx：`python tests/compat/gmx_regress.py --root <项目目录> --frames 20 --stride 1 --json out.json`
  （若本机 gmx 版本低于 tpr 版本会记 `gmx_unavailable`，属正常）。
- GROMACS 仅 lat（崩溃猎捕）：`python tests/compat/gmx_compat.py --root <项目目录> --json out.json`
- LAMMPS：`python tests/compat/compat_run.py --root <项目目录> --maxdepth 6 --json out.json`
- **F 盘是 USB 机械盘**：务必只指定具体目录并限制 `--maxdepth`，不要全盘扫描。

### 6.3 测试输入选择建议
- 至少覆盖：正交盒与三斜盒；`atomic`/`charge`/`molecular`/`full` 四种 atom_style；
  全原子与小分子/聚合物；含/不含速度列；含/不含 Bonds；单帧与多帧；大盒与“盒子接近 cutoff”的小盒。
- 重点边界：单帧轨迹、空选择、盒边小于 2×cutoff、类型缺质量/元素、粗粒化体系。

---

## 7. 判读与容差（重要）

**“照抄算法”也可能出现 >1% 的差异，绝大多数来自采样与网格约定，不是实现错误。** 判读原则：

1. **先看积分量（配位数、平均密度、Rg 均值），再看逐点最大值。**
   逐点最大值对噪声极敏感：某真实体系 MO1–OH1 的最大 |Δg| 为 0.045，但该处逐帧标准误是
   0.025（≈1.8σ），而配位数相对误差只有 0.02%。
2. **RDF 逐点差异**主要来自：参考/选择粒子数少（聚合物、离子）导致的统计噪声；bin 中心 vs
   左边缘的半 bin 约定；尖锐峰处 g 的梯度被放大。**用配位数（积分）判断，而非 max|Δg|**。
3. **密度剖面**的 max|Δ| 常出现在尖锐界面/峰（剖面跨度可达数十倍），均值偏置才是系统误差。
4. **MSD** 的 3–5% 差异来自有限时间原点采样；当帧间隔极大（位移接近半盒）时 PBC 解包本身有歧义。
5. **已知真实差异**：`gyrate`/`rms` 若不整分子复原会低估 Rg（v0.1 已修复，见算法手册第 4 节）。

如果**积分量**（CN、平均密度、Rg 均值、D）与参考软件差异 >1–2%，才应作为缺陷上报。

---

## 8. 故障排查

| 现象 | 原因/处理 |
|---|---|
| `[lat ERROR] selection '...' matched 0 atoms` | 选择写错，或用了 LAMMPS 不支持的 `resname`；先用 `lat select` 验证 |
| `unknown atom_style` | 用 `-style` 指定；加速器后缀已自动处理 |
| `no per-type masses available` | 拓扑无 `Masses`；密度改用 `-quantity number` |
| `probe: no vdW radius for [...]` | 元素无法推断；补 `Masses`、dump `element` 列或 `-types` |
| `transport: trajectory must contain velocity columns` | 轨迹无 `vx/vy/vz`（如 xtc）；用带速度的 dump/trr |
| `dump has N atoms, topology has M` | 拓扑与轨迹不是同一体系，或反应轨迹原子数变化；换正确拓扑或限定窗口 |
| `XDR read error` | 轨迹文件截断/损坏 |
| 运行“卡住” | 加 `-v` 看逐帧进度；大体系先缩小 `-frame_range` |
| 输出没更新 | 旧文件被备份为 `.bak.N`；检查 `-o` 基名 |

---

## 9. 反馈缺陷时请附上

1. 完整命令（含 `-s/-f` 的绝对路径）；
2. `lat --version` 与 `python --version`；
3. 完整 stderr（`[lat ERROR]/[lat WARNING]` 原文）；
4. `lat select -s ... -sel all` 的屏显（原子数/类型/键/帧数/盒）；
5. 期望值及来源（gmx 命令与输出、或解析解）；
6. 输入文件允许的最小复现样例（可裁剪帧数与原子数）。

---

## 10. 版本与范围

- **版本号 0.1**：功能全、接口稳定；未实现项明确标注（见算法手册“未实现/限制”）。
- 本轮已清理：无用参数（`pistack -level`）、死代码；文档合并为《算法原理手册》与本手册。
- 支持与不支持的边界见《算法原理手册》各模块“适用范围”。
