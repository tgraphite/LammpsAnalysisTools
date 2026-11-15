# 公共参数与约定

本文汇总所有分析模块共用的命令行参数、输入格式、选择语法和输出契约。每个模块的专属参数见对应手册。

## 命令形状

    lat <module> -s <topology> [-f <trajectory>] [module options] [-o <basename>]

- 单横线参数，风格对齐 GROMACS：-s 是拓扑，-f 是轨迹。
- 分析模块必须有 -s；静态分析可以不传 -f，但至少需要一帧坐标。
- 不提供交互式输入；同一命令重复运行结果一致。

## 公共参数

| 参数 | 说明 |
|---|---|
| -s, -topo | 拓扑文件：LAMMPS .data/.lmp，或 GROMACS .tpr/.gro/.pdb |
| -f, -traj | 轨迹文件：LAMMPS .lammpstrj，或 GROMACS .xtc/.trr/.gro |
| -o | 输出基名；缺省用拓扑文件名（不含扩展名） |
| -unit | 长度单位：real、metal、lj；缺省用配置 io.default_unit 并报 INFO |
| -style | 强制 LAMMPS atom_style（Atoms 段列数有歧义时使用） |
| -types | type→元素映射文件，每行 “type element” |
| -index, -n | GROMACS .ndx 文件，为选择表达式提供命名组 |
| -frame_range | 帧窗口 start:stop[:stride]，stop 不含；缺省全部帧 |
| -stride | 每隔多少帧取一帧，缺省 1 |
| -v | 在 stderr 显示 INFO 进度 |
| -debug | 额外写一份 JSON 审计文件 |
| -h | 显示该模块的完整帮助（含默认值） |
| -V, --version | 显示版本号 |

注意：PBC 默认开启且不可关闭；lat 的所有距离都是最小镜像距离。

## 输入格式

| | 拓扑 | 轨迹 |
|---|---|---|
| LAMMPS（首选） | .data / .lmp，含 Masses、Atoms；含 Bonds 才能做氢键/整分子分析 | .lammpstrj（dump） |
| GROMACS | .tpr / .gro / .pdb | .xtc / .trr / .gro |

- 自动识别 atom_style：atomic、charge、molecular、full，支持 full/kk 等加速器后缀；六列歧义时才需要 -style。
- 类型→元素默认按质量 ±0.2% 匹配；不确定时用 -types 显式指定，或让 dump 带 element 列。
- 轨迹坐标支持 x/y/z、xu/yu/zu、xs/ys/zs、xsu/ysu/zsu；带 ix/iy/iz 可解包。
- 原子数逐帧变化（反应轨迹）会报错并建议用 -frame_range/-stride 限定到一致窗口。

## 选择语法

MDAnalysis 风格子集，支持 and、or、not 与括号；相邻关键字隐含 and。

| 关键字 | 含义 | LAMMPS | GROMACS |
|---|---|---|---|
| all | 全部原子 | 是 | 是 |
| type | 原子类型 | 是 | 是 |
| element | 元素符号 | 是（需推断或 -types） | 是 |
| mol | 分子 ID | 是（需分子信息） | 是 |
| id | 原子 id | 是 | 是 |
| index | 0-based 原子序号 | 是 | 是 |
| name | 原子名 | 否 | 是 |
| resname | 残基名 | 否 | 是 |
| resid | 残基号 | 否 | 是 |

例子：

    -sel "type 2"
    -sel1 "type 1 2" -sel2 "type 3 4"
    -sel1 "resname LI" -sel2 "resname FSI"
    -sel "all and not type 1"

单选择模块（gyrate/rms/msd/sasa 等）用 -sel；双选择模块（rdf/contact/hbond/agg/cluster/pistack）用 -sel1/-sel2。先运行 lat select 验证表达式。

## 输出约定

- 结果全部写入 CSV；stdout 保持为空，进度/警告/错误走 stderr。
- 文件名：<basename>.<tag>[.<mode>].out.csv，例如 out.rdf.out.csv、out.rdf.cn.out.csv、out.density.2d.out.csv。
- 同名文件已存在时，旧文件先改名为 .bak.1、.bak.2……，不会静默覆盖。
- 文件头以 # 开头，包含 program、module、command、selections、definitions、frame_range、scalars、warning、title、xlabel、ylabel、columns，最后是 # --- data ---。
- 标量（配位数、D、R²、寿命、积分等）在 # scalars: 一行，用分号分隔。
- 退出码：0 成功；1 参数/输入/计算错误；130 用户中断。警告不会改变退出码。

## 单位

- 长度 Å；质量 g/mol；电荷 e；温度 K；密度 kg/m³（质量）或 1/Å³（数）。
- 时间轴是帧号，不解析文件中的物理时间。寿命/滞留时间单位是帧；配置 meta.frame_interval（ps）后额外给出换算标量。
- D 默认 Å²/frame；配置 frame_interval 后同时给 cm²/s。电导 S/m。

## 配置与环境变量

默认值来自 lat 配置：源码树是 config/lat.yaml，安装后是包内 data/lat.yaml。命令行覆盖配置。

| 变量 | 作用 |
|---|---|
| LAT_CONFIG | 使用自定义配置文件 |
| LAT_VENV | install.sh 创建的虚拟环境路径（默认 <repo>/.venv） |
| LAT_BIN | install.sh 放置 lat 启动器的目录（默认 ~/.local/bin） |
| PYTHON | install.sh 使用的解释器 |

## 常见错误

| 现象 | 处理 |
|---|---|
| selection matched 0 atoms | 选择写错，或用了 LAMMPS 不支持的 resname/name；先用 lat select 验证 |
| unknown atom_style | 用 -style 指定 |
| no per-type masses available | 拓扑无 Masses；密度改用 -quantity number |
| no vdW radius for [...] | 元素无法推断；补 Masses、dump 的 element 列或 -types |
| trajectory must contain velocity columns | 轨迹无 vx/vy/vz；换带速度的 dump/trr |
| dump has N atoms, topology has M | 拓扑与轨迹不是同一体系，或反应轨迹原子数变化 |
| 输出没更新 | 旧文件被备份为 .bak.N，检查 -o 基名 |

## 相关

- 模块总览与手册索引：[README.md](README.md)
- 版本变更：[../CHANGELOG.md](../CHANGELOG.md)
