# lat — LAMMPS Analysis Tool

[![PyPI](https://img.shields.io/pypi/v/lammps-analysis-tools.svg)](https://pypi.org/project/lammps-analysis-tools/)
[![Python](https://img.shields.io/pypi/pyversions/lammps-analysis-tools.svg)](https://pypi.org/project/lammps-analysis-tools/)

面向 LAMMPS `data` / `lammpstrj` 的轨迹分析工具。命令风格对齐 GROMACS（`-s` 拓扑、`-f` 轨迹、单横线参数），
输出统一为带审计头的 CSV，适合脚本批处理、回归与作图。

## 特性

- **直接读 LAMMPS**：`*.data` / `*.lmp` + `*.lammpstrj`，自动识别 `atom_style`，按质量推断元素。
- **GROMACS 兼容输入**：`.tpr/.gro/.pdb` + `.xtc/.trr/.gro` 可直接分析，便于交叉验证。
- **周期几何一致性**：所有距离、质心、整分子复原、网格与曲面都走同一套最小镜像实现。
- **可复现输出**：文件名即意图；CSV 头记录命令、选择、参数、警告与标量；同名输出自动备份。
- **自解释**：`lat --help` 列出模块，`lat <module> -h` 给出参数与默认值，`lat describe -module <name>` 给出完整契约。

## 安装

需要 Python ≥ 3.9。

    pip install lammps-analysis-tools

安装后命令是 `lat`：

    lat --version
    lat --help

升级：

    pip install --upgrade lammps-analysis-tools

从源码或离线安装：

    git clone https://github.com/tgraphite/LammpsAnalysisTools.git
    cd LammpsAnalysisTools
    bash install.sh
    export PATH="$HOME/.local/bin:$PATH"

`install.sh` 会建独立虚拟环境、安装依赖，并把 `lat` 放到 `~/.local/bin`；系统 Python < 3.9 时会用
`uv` 引导一个私有 Python 3.11。也可以直接 `pip install .`。

依赖：numpy、scipy、pyyaml、matplotlib、MDAnalysis、ase。

## 快速开始

    # 径向分布函数 + 配位数
    lat rdf -s system.data -f run.lammpstrj -sel1 "type 2" -sel2 "type 5" -cutoff 10 -o out

    # 沿 z 的密度剖面
    lat density -s system.data -f run.lammpstrj -axis z -nbins 200 -o out

    # 均方位移与扩散系数
    lat msd -s system.data -f run.lammpstrj -sel "type 2" -o out

    # 画图（输入是上一步的 CSV）
    lat plot1d -f out.rdf.out.csv -o rdf.png

结果写入当前目录或 `-o` 指定的基名：

    out.rdf.out.csv
    out.rdf.cn.out.csv

## 支持输入

| | 拓扑 | 轨迹 |
|---|---|---|
| LAMMPS（首选） | `.data` / `.lmp`：Masses、Atoms；含 Bonds 才能做氢键/整分子分析 | `.lammpstrj`（dump） |
| GROMACS | `.tpr` / `.gro` / `.pdb` | `.xtc` / `.trr` / `.gro` |

- 自动识别 `atomic`、`charge`、`molecular`、`full`，支持 `full/kk` 等加速器后缀；六列歧义时用 `-style`。
- 类型→元素默认按质量 ±0.2% 匹配；不确定时用 `-types` 显式指定，或让 dump 带 `element` 列。
- 坐标支持 `x/y/z`、`xu/yu/zu`、`xs/ys/zs`、`xsu/ysu/zsu`；带 `ix/iy/iz` 可解包。
- 反应轨迹原子数逐帧变化时会报错，并提示用 `-frame_range` / `-stride` 限定到一致窗口。

## 模块总览

| 模块 | 用途 | 手册 |
|---|---|---|
| `rdf` | 径向分布函数 g(r) 与配位数 | [rdf](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/rdf.md) |
| `density` | 一维/二维密度剖面 | [density](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/density.md) |
| `msd` | 均方位移与扩散系数 | [msd](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/msd.md) |
| `gyrate` | 回转半径 Rg | [gyrate](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/gyrate.md) |
| `rms` | RMSD / RMSF | [rms](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/rms.md) |
| `contact` | 两组间接触数与距离 | [contact](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/contact.md) |
| `hbond` | 氢键数量、寿命与空间分布 | [hbond](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/hbond.md) |
| `agg` | SSIP / CIP / AGG 离子聚合形态 | [agg](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/agg.md) |
| `cluster` | 配位数与团簇尺寸分布 | [cluster](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/cluster.md) |
| `pistack` | π–π 堆积 | [pistack](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/pistack.md) |
| `freevol` | 自由体积、孔径分布、连通性 | [freevol](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/freevol.md) |
| `sasa` | 溶剂可及表面积 | [sasa](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/sasa.md) |
| `residence` | 配体滞留时间 | [residence](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/residence.md) |
| `transport` | 扩散、Nernst–Einstein、Green–Kubo | [transport](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/transport.md) |
| `occupancy` | 区域占据、穿越、首次通过 | [occupancy](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/occupancy.md) |
| `select` | 检查选择、导出 .ndx | [select](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/select.md) |
| `describe` | 模块清单与契约 | [describe](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/describe.md) |
| `plot1d` | 一维 CSV 绘图 | [plot1d](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/plot1d.md) |
| `plot2d` | 二维 CSV 绘图 | [plot2d](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/plot2d.md) |

完整手册（公共参数、输入、选择语法、输出契约）：[docs/manual](https://github.com/tgraphite/LammpsAnalysisTools/tree/main/docs/manual)。

## 公共参数

| 参数 | 说明 |
|---|---|
| `-s, -topo` | 拓扑文件（分析模块必填） |
| `-f, -traj` | 轨迹文件 |
| `-o` | 输出基名；缺省用拓扑文件名 |
| `-unit` | 长度单位：`real`、`metal`、`lj`；缺省用配置 |
| `-style` | 强制 LAMMPS `atom_style` |
| `-types` | type→元素映射文件，每行 `type element` |
| `-index, -n` | GROMACS `.ndx`，为选择提供命名组 |
| `-frame_range` | 帧窗口 `start:stop[:stride]`，stop 不含 |
| `-stride` | 每隔多少帧取一帧 |
| `-v` | 显示 INFO 进度 |
| `-debug` | 额外写 JSON 审计 |
| `-h` | 模块帮助（含默认值） |
| `-V, --version` | 版本号 |

PBC 始终开启，所有距离使用最小镜像。时间轴是帧号，不解析文件中的物理时间。

## 选择语法

MDAnalysis 风格子集，支持 `and`、`or`、`not` 与括号；相邻关键字隐含 `and`。

| 关键字 | 含义 | LAMMPS | GROMACS |
|---|---|---|---|
| `all` | 全部原子 | 是 | 是 |
| `type` | 原子类型 | 是 | 是 |
| `element` | 元素符号 | 是（需推断或 `-types`） | 是 |
| `mol` | 分子 ID | 是（需分子信息） | 是 |
| `id` | 原子 id | 是 | 是 |
| `index` | 0-based 原子序号 | 是 | 是 |
| `name` | 原子名 | 否 | 是 |
| `resname` | 残基名 | 否 | 是 |
| `resid` | 残基号 | 否 | 是 |

    -sel "type 2"
    -sel "type 2 or type 5"
    -sel1 "resname LI" -sel2 "resname FSI"
    -sel "all and not type 1"

单选择模块用 `-sel`；双选择模块用 `-sel1`/`-sel2`。先运行 `lat select` 验证命中。

## 输出

- 结果全部写入 CSV；stdout 为空，进度/警告/错误走 stderr。
- 文件名：`<basename>.<tag>[.<mode>].out.csv`。
- 文件头以 `#` 开头，包含 program、module、command、selections、definitions、frame_range、
  scalars、warning、title、xlabel、ylabel、columns，最后是 `# --- data ---`。
- 标量（首壳配位数、D、R²、寿命、积分等）在 `# scalars:` 行，分号分隔。
- 同名文件先备份为 `.bak.1`、`.bak.2`……，不会静默覆盖。
- 退出码：0 成功，1 参数/输入/计算错误，130 中断；警告不改变退出码。

## 绘图

`plot1d` / `plot2d` 只读 CSV，计算与绘图解耦。

    lat plot1d -f out.rdf.out.csv -o rdf.png
    lat plot1d -f out.agg.labels.out.csv -columns A_SSIP,A_CIP,A_AGG -o spec.png
    lat plot1d -f out.agg.labels.out.csv -stack -o spec_stack.png
    lat plot2d -f out.density.2d.out.csv -o density2d.png

`-smooth`（`none`/`gaussian`/`savgol`）只在绘图阶段生效，CSV 始终是原始值；
`-bar` 与 `-stack` 互斥，且不做平滑。

## 配置

默认值来自 `lat.yaml`：源码树是 `config/lat.yaml`，安装后是包内 `data/lat.yaml`。
命令行覆盖配置。

| 环境变量 | 作用 |
|---|---|
| `LAT_CONFIG` | 使用自定义配置文件 |
| `LAT_VENV` | `install.sh` 创建的虚拟环境路径 |
| `LAT_BIN` | `install.sh` 放置 `lat` 的目录 |
| `PYTHON` | `install.sh` 使用的解释器 |
| `LAT_NO_UPDATE_CHECK` | 设为 `1` 关闭自动更新检查 |
| `LAT_UPDATE_URL` | 自定义 JSON 版本端点 |
| `LAT_UPDATE_INDEX` | 自定义 PEP 503 简单索引（如内网镜像） |
| `LAT_UPDATE_INTERVAL_HOURS` | 覆盖检查间隔（默认 24 小时） |
| `LAT_UPDATE_TIMEOUT_S` | 覆盖单源网络超时（默认 1.5 秒） |

## 更新检查

lat 在交互终端启动时做一次**只提示、不安装**的版本检查：

- 默认每 24 小时最多检查一次；结果缓存在 `$XDG_CACHE_HOME/lat/update.json`。
- 检查在后台线程执行，单源超时 1.5 秒，不阻塞命令，也不影响退出码。
- 同时查询 PyPI、`PIP_INDEX_URL`、`LAT_UPDATE_URL` / `LAT_UPDATE_INDEX`，取最高版本。
- 只在 stderr 打印一行提示，例如：
  `[lat NOTICE] new version 1.2.0 available (current 1.1.0); upgrade with: pip install -U lammps-analysis-tools`
- 非交互终端、`CI` 环境、`--version` / `-h` / `describe` 自动跳过。
- 关闭：`LAT_NO_UPDATE_CHECK=1`，或配置 `update.enabled: false`。
- 手动检查：`lat --check-update`；单次跳过：`lat --no-update-check <module> ...`。

lat 不会自动下载或执行任何更新；升级始终由你通过 pip/uv 完成。

## 文档

- 逐模块手册：[docs/manual](https://github.com/tgraphite/LammpsAnalysisTools/tree/main/docs/manual)
- 公共参数与输出契约：[common.md](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/manual/common.md)
- 版本变更：[CHANGELOG.md](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/docs/CHANGELOG.md)

## 常见问题

**选择命中 0 个原子？** 用 `lat select` 检查表达式；LAMMPS 拓扑通常没有 `resname`/`name`。

**元素推断失败？** 补 `Masses`，或让 dump 带 `element` 列，或用 `-types`。

**氢键报 “no donor”？** 给体只认 N/O/F/S/Cl/Br/I 上的 H；C–H 不是给体。

**MSD 提示没有 image flags？** 轨迹只有 wrapped 坐标，MSD 可能饱和；改用带 `xu/yu/zu` 或 `ix/iy/iz` 的 dump。

**freevol 报三斜盒？** 体素网格只支持正交且固定的盒子；改用分子/原子对分析或正交化盒子。

## 许可

Proprietary。详见 [LICENSE](https://github.com/tgraphite/LammpsAnalysisTools/blob/main/LICENSE)。
