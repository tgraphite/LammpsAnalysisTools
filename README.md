# lat — LAMMPS Analysis Tool

> 面向 LAMMPS `data` / `lammpstrj` 的统一后处理工具。
> 命令风格对齐 GROMACS（`-s` 别名 topo、`-f` 别名 traj，单横线，无 `--topo`），
> 但语义、输出契约与文档体系完全自有。

## 一句话定位

把「分子动力学后处理」中高频、可复用、可自动化的分析，收敛成一组**边界清晰、自解释、可回归**的子命令，
让小白与 AI 都能 `lat <module> -h` 秒上手，让工程师只关心「选哪组原子、用什么判据」。

## 主旨（不可动摇的四条）

1. **完全兼容** LAMMPS 的 `data` / `lammpstrj`，并通过 GROMACS Ground Truth 做两级回归。
2. **自解释**：`lat -h` / `lat <module> -h` 输出用法、参数、默认值与手册入口；原理与判读见
   《算法原理手册》《测试人员手册》——文档即教材。
3. **数据与逻辑分离**：所有默认值、单位、输出模板、模块文档集中在唯一 `config/lat.yaml`；
   代码中**不得**存在任何 backup 默认值，配置错误直接崩溃是预期行为。
4. **格式无关的 IO**：内部基于 MDAnalysis 读取器，可直接吃 GROMACS 轨迹（**不对外暴露**，仅用于回归与内部复用）。

## 不在职责内（交给专门工具）

固体/金属力学（应力应变、位错、DXA、压痕、切削）、相互作用能/自由能、SDF 图、谱（VDOS/SAXS）、
反应网络、对接、PMF/WHAM、量子/介观。这些要么属于 LAMMPS on-the-fly，要么属于 OVITO/Zeo++ 等专用工具。

## 安装（一键，新设备可用）

```bash
cd lat && bash install.sh        # 建独立环境、装依赖、把 lat 放到 ~/.local/bin
export PATH="$HOME/.local/bin:$PATH"
lat --version                    # lat 0.1.0
```

无需 root、无需 `source`/激活、无需在源码目录下运行。系统 Python < 3.9 时安装脚本会自动用
`uv` 装一个私有 Python 3.11。环境变量：`LAT_CONFIG`（自定义配置）、`LAT_VENV`、`LAT_BIN`、`PYTHON`。

## 目录结构

```
lat/
├── README.md                 本文件
├── install.sh                一键安装脚本
├── config/
│   └── lat.yaml              唯一配置来源（默认值/单位/文档/输出模板）
├── docs/
│   ├── 00_provenance.md      来龙去脉与决策记录
│   ├── 01_meta_requirements.md 元要求（M1–M15）与范围
│   ├── 02_architecture.md    系统架构
│   ├── 03_io_selection.md    IO、拓扑、选择器、PBC、单位
│   ├── 04_output_contract.md 输出契约
│   ├── 05_regression.md      回归测试方案
│   ├── 06_reference.md       参考代码与权威来源
│   ├── 07_decision_log.md    已确认决策与优先级
│   ├── 08_handoff.md         compact 后的恢复/交接笔记（先读它）
│   ├── MANUAL_TESTER.md      **测试人员手册**（安装/逐模块操作/判读/排错）
│   ├── MANUAL_ALGORITHMS.md  **算法原理手册**（教材式，逐模块原理与精度）
│   ├── 09_test_campaign.md   测试战役报告（静态分析/幂等/文档/真实项目）
│   ├── 10_refdata_campaign.md 第三方 RefData 订单 + gmx 原生回归
│   ├── 11_examples_recheck.md 0.2.0 对已有项目产物的复跑/覆盖报告
│   ├── 12_0.2.0_release.md   0.2.0 发布说明（行为变更清单 + 验收）
│   ├── data/                 产物级归档（0.2.0_examples_mapping.tsv 等）
│   ├── figures/0.2.0/        修复前后对比图（16 张，见 11 号文档 §6）
│   └── BYPASS.md             旁路记录：讨论过但暂不纳入的功能及理由（非 roadmap）
├── src/lat/                  实现（P0 起逐步填充）
│   ├── core/                  IO/拓扑/选择/PBC/网格/相关/输出
│   ├── modules/               各子命令
│   └── plotting/              plot1d / plot2d
└── tests/
    ├── regression/
    │   ├── convert/          自研 gmx→lammps 转换器 + 校验（无需 InterMol/VMD）
    │   ├── examples_recheck/ 复跑/覆盖/对比图工具（行为变更后重算已有产物）
    │   ├── samples/          已落盘自包含回归样例（sys1_reduced）
    │   ├── level1_algo/      原生 gmx → 算法回归
    │   └── level2_format/    转换后 → 格式回归
    └── data/                 小型自造测试体系
```

## 当前状态

**全部 P0 模块已实现**。`lat <module> -s <data> [-f <traj>]` 可运行：

```
core    io / topology / selection / pbc / neighbor / grid / correlate / output / context / config
模块    rdf density hbond rms msd agg cluster pistack freevol sasa residence
        transport occupancy gyrate contact select describe
绘图    plot1d plot2d
```

回归：`pytest tests` **161 passed / 5 skipped**（含 level2；level1 逐值对拍为独立脚本，需本机 gmx）。四级证据：
- **单元/合成**（PBC、密度、MSD、RMSD、手算氢键/SSIP-CIP-AGG/ππ、FFV、滞留、占据、输运守卫，
  **键表 id 重排**）；
- **E2E 冒烟**（真实 6600 原子样例逐模块跑通）；
- **level1 对拍 gmx（GROMACS 原生输入）**（sys1–6 含 PMBA 聚合物：RDF/CN 相对误差 ≤0.2%，
  密度 <1%，MSD ~5%）；
- **level1 对拍 gmx（LAMMPS 输入，2026-10 新增）**：`convert/lammps_to_gmx.py` +
  `level1_algo/gmx_ref_lammps.py`，在 4 链 PET（6200 原子 / 501 帧）上逐帧比
  gyrate / rms / rmsf / sasa / rdf / hbond。结果：Rg 3e-4 Å、RMSD 1e-4 Å、RMSF 7e-4 Å、
  RDF 配位数 0.26%、hbond 0.21%（详见 `tests/regression/level1_algo/README.md`）。

**0.2.0 修掉的四个缺陷**（都会改数值，`docs/12_0.2.0_release.md` 是完整清单）：
1. `Bonds` 段写的是 **atom id** 而非 Atoms 段行号，Atoms 不按 id 排序时键图错乱 →
   `gyrate`/`rms`/`hbond` 静默出错（OVITO 写的 PET：Rg 高估 3.6×）；
2. `hbond` 把每一个 C–H 当给体 → 计数虚高（PET 体系 41.7 vs gmx 4.3/帧）；
3. `sasa` 未用 PBC → 系统性高估 +14.2%（`gmx sasa` 对拍）；
4. PSD 的 bin 宽/归一化/周期距离变换、connectivity 的环面与邻接全部重做。

同时把 `probe`（体素网格 + 曲面混在一起）拆成 **`lat freevol`**（ffv/psd/connectivity）
与 **`lat sasa`**（Shrake–Rupley）；`probe` 已删除，不留兼容壳。
所有周期几何统一到 `lat.core.periodic`，并由 `tests/unit/test_no_ad_hoc_pbc.py` 守卫。
尚待：`current` 与 gmx 逐值对拍；P1 审计。

命令示例：
```
./lat select -s sys.data -sel "type 2"
./lat rdf -s sys.data -f sys.lammpstrj -sel1 "type 2" -sel2 "type 5" -cutoff 10 -o out/sys
./lat plot1d out/sys.rdf.out.csv
```

剩余 P0 收尾（`current`、P1 审计）按 `docs/08_handoff.md` 的顺序继续。
P1 = 待 P0 验收后开发；P2 = 关键点验证后再开发（当前 `order`、`occupancy` 的分子体积自动 zone）。

## 快速约定

- 时间 = **帧号**，LAMMPS 时间字段一律忽略；寿命/滞留时间均以「帧」为单位输出，
  并在 header 注明「= k × 帧间距」。
- 不做任何分子/残基推断；计算目标必须由用户在运行前用静态 MDAnalysis 语义或 `ndx` 明确给出。
- 无交互式输入（保证幂等）。
