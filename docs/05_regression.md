# 05 · 回归测试方案

## 1. 两级输入回归

```
输入 A（原生 GROMACS: tpr/gro/xtc）
        │  MDAnalysis 直读（内部能力，不对外）
        └─> lat 模块 ──> 与 gmx 生成的 .xvg 数值比对     ← 验证【算法】
输入 B（同一体系转 LAMMPS data + lammpstrj）
        └─> lat 模块 ──> 与同一份 gmx .xvg 数值比对       ← 验证【格式解析 + 算法】
```

- 两级必须一致，才能证明「读 LAMMPS」与「读 GROMACS」得到同一物理量。
- 转换**已完成**，不依赖 InterMol/VMD：自研 MDAnalysis 转换器
  `tests/regression/convert/convert_gmx_to_lammps.py`（写 LAMMPS data + lammpstrj；
  MDAnalysis 本身只读不写 LAMMPS，故自实现 writer）。
- **不支持多文件/glob 拼接**：一轨迹一模拟；回归按单个 dump 组织。

## 2b. 已落盘自包含样例（`tests/regression/samples/sys1_reduced/`）

- 由 sys1 缩减（每 resname 200 分子 → 800 分子/6600 原子；帧 `[0,20,40,60,80]`）。
- 内容：`*.data`、`*.lammpstrj`、`*.gro`、`*.xtc`、`*.mapping.csv`、`*.selections.yaml`、`*.verify.json`。
- 自检全通过：坐标逐帧最大偏差 `9.5e-7 Å`；元素计数一致；`type 2`↔`resname LI` 各 200；
  净电荷 `2e-4 e`；键 6000。详见样例 `README.md`。
- 已知限制：`.gro` 不带质量/电荷（MDAnalysis 把 `LI` 质量猜为 0）；质量相关量以 `*.mapping.csv` 与 LAMMPS data 为准。
- 运行：`python convert/verify_lammps_roundtrip.py -p ...`（见 convert/）。

## 2. 主数据集：`/mnt/f/CRO/40_multi_electrolytes/11-systems`

- 6 体系，183k–188k 原子，正交盒。sys1–3 纯电解液（Li/FSI/FEC/FEMC/TTE），**sys4/sys6 含 PMBA 4 链**（聚合物 + 电解液）。
- 现成 Ground Truth：
  - `rdf/`：多 `ref/sel` 组合，bin=0.005 nm，含 `-cn` 配位数（`sys*_cn_fine_LI.xvg`）。
  - `density/`：`sys*_density.xvg`。
  - `msd/`：`sys*_npt_msd_{LI,FSI,FEC,FEMC,TTE,MBA,PMBA}.xvg`（含拟合 D 注释）。
- 覆盖模块：`rdf`、`density`、`msd`；`agg`/`cluster` 可用 Li–FSI 组合做一致性检查（无独立 GT，用自洽与转换一致性）。

## 3. 数据缺口与补测

| 模块 | 缺什么 | 补测方案 |
|---|---|---|
| `hbond` | 该数据集无 O–H 给体 | **自造计算对象**：几个带 A–H···B 的分子，手算期望点数；只验数值一致性，不考虑物理。**已有 gmx GT（2026-10）**：`level1_algo/gmx_ref_lammps.py` 在带端羟基的 PET 上与 `gmx hbond` 逐帧比，501 帧均值差 0.014、最大 1 条键 |
| `pistack` | 无 ππ 体系 | 移植代码，沿用 `/mnt/d/dev/83_pipi_stack_code/tests` 与苯液体 |
| `freevol` | 无多孔/膜 | 自造多孔网格或均相液体（FFV）。PSD 有**内部一致性**断言：`∫dV/dr dr = 1`（`test_periodic_core.py` 用暴力环面 EDT 校验距离变换）；connectivity 用环面闭合用例（`test_label_periodic_*`）。gmx 无对应工具 |
| `sasa` | 无多孔/膜 | **已有 gmx GT（2026-10）**：PET 体系逐帧对 `gmx sasa`，修正 PBC 后 frame 0 差 −0.17%、frame 1/2 差 −0.83%/−0.85%（0.1.0 为 +14.2%，见 `MANUAL_ALGORITHMS.md` §12.5）|
| `residence` | — | 自造「一中心 + 若干配体」体系，手算滞留区间 |
| `transport` | 需速度列 | dump 加 `vx vy vz`；与 msd 的 D 交叉验证 |
| `occupancy` | — | 自造两区域体系，手算穿越次数 |
| `gyrate`/`contact`/`order` | — | 用 40 体系（PMBA）做冒烟 + 手算小例。**gyrate 已有 gmx GT**：LAMMPS 路径逐链 Rg 差 2e-4 Å（`level1_algo/README.md`） |

### 3.1 LAMMPS 输入的 gmx Ground Truth（2026-10 新增）

上面第 2 节的 GT 全部来自 **GROMACS 原生输入**，LAMMPS 的 data/dump 路径此前没有原生参考。
现在补上：`convert/lammps_to_gmx.py`（LAMMPS → gmx 拓扑/轨迹）+ `level1_algo/gmx_ref_lammps.py`
（对拍 gyrate/rms/rmsf/sasa/rdf/hbond）。用法与结果见 `level1_algo/README.md` 与
`09_test_campaign.md` §7。

## 4. 容差

- RDF/密度/MSD：与 gmx xvg 相对误差 `≤ 1e-3`（同一轨迹、同 bin、同归一化）。
- CN：绝对误差 `≤ 1e-2`。
- 计数类（hbond/pistack/occupancy）：整数必须**完全相等**。
- 自相关/寿命：拟合窗口一致时相对误差 `≤ 1e-2`；帧数不足允许 FATAL/WARN 而非数值。

## 5. 口径提醒

- 回归**只考虑程序一致性，不考虑物理价值**。
- 固体/金属力学不在范围，不建回归。
- 转换后的 LAMMPS 体系若改变原子顺序/类型编号，需在 `tests/regression/convert/` 记录映射，保证选择器指向同一组原子。

## 6. 测试分层

```
tests/
├── conftest.py          定位 repo/样例、注入 src
├── regression/
│   ├── convert/         gmx→lammps 转换脚本 + 映射记录
│   ├── level1_algo/     原生 gmx 输入 → 比对 xvg（待 rdf 落地后补）
│   ├── level2_format/   转换后 lammps 输入 → 比对 xvg / 模块 E2E
│   │   ├── test_rdf_sys1_reduced.py
│   │   └── test_modules_sys1.py       所有模块在真实样例上的端到端冒烟
│   └── expected/        从 11-systems 复制的期望值（精简子集）
├── data/                自造小体系（hbond/residence/occupancy/freevol/sasa）
└── unit/                PBC/网格/归一化/自相关
    └── test_core.py
```

## 7. 早期回归结果（历史记录：仅 core/select/describe/rdf 阶段）

> **本节是 0.1.0 早期阶段的快照，不代表当前状态。** 当前状态见 §3.2、§8 与 `README.md`
> （`pytest tests` **161 passed / 5 skipped**）。

当时已实现 `core` + `select`/`describe`/`rdf`，回归状态：

| 检查 | 结果 |
|---|---|
| LAMMPS(data+lammpstrj) vs 内部 GROMACS(gro+xtc) 的 g(r) | **逐点完全相同**（max Δ = 0） |
| lat g(r) vs MDAnalysis InterRDF（同 dump，200 bins，0–10 Å） | `r` 轴一致；`max|Δg| = 4.7e-2`（首峰处 0.5%），`mean = 2.4e-3` |
| lat CN vs 由 MDAnalysis 对数重建的 CN | 完全相等（`atol 1e-12`） |
| 单测：正交/三斜最小镜像、元素质量匹配（含 reax_tools Li=7.0 必须失败）、选择器文法、备份轮转 | 全部通过 |

- 当时的命令：`/opt/miniconda3/bin/python -m pytest tests -q` → **6 passed**。
- 说明：该样例为**缩减稀释体系**（每 resname 200 分子），配位数绝对值不具物理意义，仅用于
  格式/算法一致性。
- 当时尚未接入 gmx `.xvg` GT；这一项后来由 §3.1/§3.2 与 `level1_algo/` 补齐。

## 8. 端到端（E2E）覆盖现状

`tests/regression/level2_format/test_modules_sys1.py` 在真实缩减样例（6600 原子/5 帧）上逐模块跑 CLI：

| 层 | 内容 | 状态 |
|---|---|---|
| **核心数值回归** | rdf：LAMMPS vs 内部 GROMACS 逐点相同；vs MDAnalysis InterRDF | ✅ 已接入 |
| **模块 E2E 冒烟** | `select describe rdf density(1d/2d) gyrate contact rms msd hbond agg cluster pistack freevol(ffv/connectivity) sasa residence occupancy` | ✅（CLI 退出 0 + 产出 CSV） |
| **守卫用例** | `transport` 无速度列必须 ERROR | ✅ |
| **几何边界** | hbond 在「部分帧 0 氢键」时不崩（E2E 抓到的真实 bug，已修） | ✅ |
| **单测（合成体系）** | 20 项：PBC/选择/元素/输出、密度、MSD、RMSD、hbond、SSIP-CIP-AGG、ππ、FFV、滞留、占据、输运守卫 | ✅ |
| **level1 gmx .xvg golden** | 全量 11-systems 的 `rdf/density/msd` 对 gmx | ⬜ 未接入（下一步） |

- 合计 `pytest tests` → **161 passed / 5 skipped**（level1 的逐值对拍是独立脚本，需要本机 gmx）。
- **level1 gmx 逐值回归已接入**：`tests/regression/level1_algo/compare_gmx.py` 运行时用 `gmx rdf` 生成 GT，
  对 sys1–sys6（含 PMBA 聚合物）比较。结果：max\|Δg\|≤0.184（首峰<1.5%）、rms≤0.024、mean~1e-5、
  **配位数相对误差≤0.2%**；密度 sys1 max 0.95%/rms 0.24%；MSD（mol_com）LI/FSI 曲线 ~3–7%。
  详见 `docs/09_test_campaign.md`。
- 仍属冒烟而非数值 GT 的模块：`density msd rms agg cluster pistack freevol residence transport occupancy`
  （已在合成体系上做数值断言，但未与 gmx/外部实现逐值对拍）。

### 3.2 PBC 平移不变性（2026-10 新增，0.2.0）

`tests/regression/level2_format/test_pbc_invariance.py` 把整个体系做一次**非格点平移**（+7.3 Å），
逐列比较模块输出：物理量必须**逐位相同**。覆盖
`contact`（min/center 两种模式）、`residence`（min/center）、`pistack`、`agg`、`cluster`、`rdf -intramol`。

在 0.1.0 代码上该测试的失败项：`contact -mode min/center` 的 `com_dist`（4.0 Å vs 6.0 Å，
真值 6.0 Å）、`residence -mode center` 的绑定计数、`pistack` 的堆积计数。
`rdf -intramol` 与 `agg/cluster` 通过（前者统计成对 MIC 距离，后者的小分子整体在参考原子半盒内）。

同时新增：
- `tests/unit/test_periodic_core.py`：环面距离变换（对 9×8×7 非均匀网距做 O(N²) 暴力校验）、
  环面连通分量（跨面合并/角落 26 邻接/贯通判定）、周期 SASA（单原子解析值、跨界不变性、
  跨面屏蔽）、整分子复原（跨面链、长于半盒的链、子集只遍历选中分子）、质心歧义。
- `tests/unit/test_no_ad_hoc_pbc.py`：架构守卫（模块不得自建 `cKDTree`/EDT/`ndimage.label`/`np.pad`/
  裸 `minimum_image`；PBC 敏感模块必须走 `core.periodic`；禁止裸 `mean(axis=0)` 求质心）。
