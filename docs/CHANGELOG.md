# Changelog

本文件记录 lat 的版本变更。开发文档只保留本文件与 [MODULE_STATUS.md](MODULE_STATUS.md)；面向用户的使用说明见 [README.md](../README.md)，逐模块手册见 [manual/](manual/README.md)。

## [1.0.0] — 2026-11-15

首个公开发布版本。

### 新增

- 以 PyPI 包名 lammps-analysis-tools 发布，命令行入口仍为 lat。
- 每模块一份帮助手册（docs/manual/），参数表由代码元数据生成；新增静态校验测试，手册与模块参数/输出不一致时测试失败。
- 面向用户的 README：安装、模块总览、公共参数、输入格式、选择语法、输出契约与常见问题。
- contact 的 -mode center 现在真正按残基质心建接触图；此前该参数被接受但不影响结果。
- 打包元数据：项目主页、仓库、问题跟踪、分类器、关键词、MANIFEST.in 与 LICENSE。

### 变更

- 版本号 0.2.0 → 1.0.0；配置文件 meta.version 同步。
- 文档体系收敛为 CHANGELOG + MODULE_STATUS + 逐模块手册；删除内部过程文档、战役报告、决策日志与对比图。
- lat describe 在源码树和 pip 安装后都能定位到手册（源码 docs/manual，安装包 data/manual）。
- CLI 帮助补齐公共参数说明；所有模块参数都有 help 文本。
- 默认值配置文件删除未使用的键（io.pbc/unwrap/triclinic、rdf.normalization/compute_coordination、agg.edge_conservation_tol、transport.require_velocity/remove_com/methods、modules.order），保留的键全部有代码路径。

### 修复

- contact -mode center 与文档不一致：现在 center 模式按残基质心距离计数，min 模式按原子对计数；两个选择相同时按无序对处理，避免自配对。
- 参数帮助中带竖线的文本会破坏 Markdown 参数表；已改写为分号/逗号描述。
- 删除无效的 -pbc/-nopbc 开关（PBC 始终开启，开关从未生效）。
- 清理仓库中误提交的 __pycache__/*.pyc、开发用启动脚本与内部回归工具。

### 移除

- docs/00–12、BYPASS.md、module_contracts.md、MANUAL_ALGORITHMS.md、MANUAL_TESTER.md、docs/data、docs/figures：内容已合并进 CHANGELOG、MODULE_STATUS 与逐模块手册。
- 仓库根目录的 lat 开发包装脚本；安装后使用 pip 生成的 lat 入口，源码树使用 python -m lat 或 install.sh。
- tests/regression 中仅用于历史战役的 convert/examples_recheck/level1_algo 工具与 tests/compat 脚本。

## [0.2.0] — 2026-10

相对 0.1.0 的行为变更与缺陷修复。该版本没有公开发布，作为 1.0.0 的行为基线。

### 修复（会改变数值）

- Bonds 段按 atom id 解释，并在 Atoms 段不按 id 排序时同步重排；修复 gyrate/rms/hbond 的键图错乱（OVITO 写的 PET 上 Rg 曾高估 3.6×）。
- hbond 只把 N/O/F/S/Cl/Br/I 上的 H 当作给体；C–H 不再算给体（PET 上曾 41.7 vs gmx 4.3 每帧）。
- sasa 的邻域搜索与遮挡判定改为最小镜像；修复贴面原子丢邻居导致的系统性高估（PET 上曾 +14.2%）。
- PSD 改为周期距离变换、固定 bin 宽、按自由体积归一、26 连通平台去重；积分从 1.5e-4 修正到 1。
- connectivity 改为环面闭合、默认 26 邻接、逐帧汇总；PET 分量数从 165 降到 35。
- contact/residence -mode center 不再对 wrapped 坐标直接求平均；pistack 环心/法向、agg/cluster 质心统一走周期几何入口。
- hbond 连续自相关改为增量算法；601 帧、120k 对的 ACF 从分钟级降到 0.04 s，数学结果不变。
- 键图跨盒不一致时给出 WARN 并退回可用近似，而不是悄悄给出位移 23 Å 的质心。

### 新增

- lat.core.periodic：分子/原子对/体素网格/曲面四层周期几何唯一入口。
- lat freevol -mode ffv|psd|connectivity 与 lat sasa；删除 lat probe，不留兼容壳。
- 架构守卫 test_no_ad_hoc_pbc.py、周期核心单测、配置同步测试、平移不变性回归。
- 0.2.0 行为基线：`lat --version` 输出 0.2.0；config/lat.yaml 与 src/lat/data/lat.yaml 逐字节一致。

## [0.1.0] — 2026-09

内部首个可用版本。

### 新增

- LAMMPS data/lammpstrj 读取、类型/质量/电荷/键/分子拓扑、选择表达式、PBC 最小镜像。
- 模块：rdf、density、hbond、rms、msd、agg、cluster、pistack、freevol（原 probe）、sasa、residence、transport、occupancy、gyrate、contact、select、describe。
- 统一 CSV 输出契约、原子写、同名备份、stderr 日志、lat.yaml 单一默认值来源。
- plot1d/plot2d 绘图；install.sh 一键安装（含 uv 私有 Python 引导）。
