# 模块现状

面向开发与维护：记录每个模块的状态、契约、已固定的行为基线和已知限制。用户文档见 [README.md](../README.md)，逐模块手册见 [manual/](manual/README.md)，版本变更见 [CHANGELOG.md](CHANGELOG.md)。

状态口径：

- stable：接口与数值口径已固定，有回归测试；破坏性变更必须升 major 并写 CHANGELOG。
- beta：可用，但边界条件或数值口径仍可能调整。
- planned：已设计、未实现；不得在文档或 CLI 中宣称可用。

当前版本 1.1.0。全部已注册模块均为 stable。

## 模块总表

| 模块 | 状态 | 起始版本 | 拓扑 | 主要输出 | 关键测试 |
|---|---|---|---|---|---|
| rdf | stable | 0.1.0 | 必需 | rdf.out.csv、rdf.cn.out.csv | test_core、level2_format、level1 gmx 对拍 |
| density | stable | 0.1.0 | 必需 | density.out.csv、density.2d.out.csv | test_density_msd_rms、level2 |
| msd | stable | 0.1.0 | 必需 | msd.out.csv | test_density_msd_rms、test_core |
| gyrate | stable | 0.1.0 | 必需 | gyrate.out.csv | test_bond_atom_order、level1 gmx 对拍 |
| rms | stable | 0.1.0 | 必需 | rmsd.out.csv、rmsf.out.csv | test_density_msd_rms、level1 |
| contact | stable | 0.1.0 | 必需 | contact.out.csv | test_pbc_invariance、test_idempotency |
| hbond | stable | 0.1.0 | 必需 | hbond.out.csv、hbond.acf、hbond.dist_* | test_hbond、level1 |
| agg | stable | 0.1.0 | 必需 | agg.out.csv、agg.labels、agg.speciation | test_agg、test_pbc_invariance |
| cluster | stable | 0.1.0 | 必需 | cluster.cn.out.csv、cluster.size.out.csv | test_agg、test_pbc_invariance |
| pistack | stable | 0.1.0 | 必需 | pistack.out.csv、pistack.pairs.out.csv | test_pistack、test_pbc_invariance |
| freevol | stable | 0.2.0 | 必需 | freevol.out.csv、freevol.psd、freevol.connectivity | test_probe_residence_occupancy、test_periodic_core |
| sasa | stable | 0.2.0 | 必需 | sasa.out.csv | test_periodic_core、level1 gmx 对拍 |
| residence | stable | 0.1.0 | 必需 | residence.out.csv、residence.acf、residence.intervals | test_probe_residence_occupancy、test_pbc_invariance |
| transport | stable | 0.1.0 | 必需 | transport.d、transport.ne、transport.gk | test_density_msd_rms、test_idempotency |
| occupancy | stable | 0.1.0 | 必需 | occupancy、occupancy.crossings、occupancy.first_passage、occupancy.enrichment | test_probe_residence_occupancy |
| select | stable | 0.1.0 | 必需 | 屏幕输出；可选 .ndx | test_core、test_cli_contract |
| describe | stable | 0.1.0 | 不需要 | 屏幕输出 | test_cli_contract |
| plot1d | stable | 0.1.0 | 不需要 | 图片 | test_plot1d、test_cli_contract |
| plot2d | stable | 0.1.0 | 不需要 | 图片 | test_plot1d |
| update | stable | 1.1.0 | 不需要 | stderr 版本提示 + 本地缓存 | test_update |

## 已固定的行为基线

1. 时间 = 帧号；LAMMPS 时间字段忽略。寿命/滞留时间以帧为单位；配置 meta.frame_interval 后才给 ps 换算。
2. 所有周期距离走 lat.core.periodic；模块中不允许自建 cKDTree/distance_transform_edt/ndimage.label/np.pad/裸 minimum_image，由 tests/unit/test_no_ad_hoc_pbc.py 守卫。
3. Bonds 段按 atom id 解释；Topology 的逐原子数组按 id 升序存放；坐标与键索引必须同置换。
4. hbond 给体与受体都是 {N,O,F,S,Cl,Br,I}；C–H 不是给体；无给体直接报错。
5. freevol 网格只支持正交且固定的盒子；三斜直接报错，盒漂移 >2% 报错，0.1%–2% WARN 且网格钉在首帧。
6. sasa 的表面始终由全部原子构建；-sel 只决定 sasa_sel 报告哪部分。
7. 输出文件名 <base>.<tag>[.<mode>].out.csv；同名文件先备份为 .bak.N；CSV 头字段固定。
8. 默认值只来自 config/lat.yaml；代码不得为物理参数提供备用默认值。源码配置与包内配置逐字节一致（test_config_sync）。
9. 所有模块参数必须有 help 文本；docs/manual 的参数表由 lat.docs.render_param_table 生成，并由 test_manuals 静态校验。
10. 更新检查只读版本号，不下载/不执行远程内容；默认仅在交互终端、最多每 24 小时检查一次，缓存写 $XDG_CACHE_HOME/lat/update.json；任何网络失败都不得影响命令退出码。

## 已知限制

| 模块 | 限制 |
|---|---|
| rdf | cutoff 超过半盒长会失真；NPT 下 g 与 cn 只是同口径近似自洽；无分子 ID 时无法排除同分子对。 |
| density | 三斜/畸变盒的截面归一化是近似；质量密度需要质量；平滑只在绘图阶段。 |
| msd | 无 image flags/xu 时可能饱和；位移接近半盒时解包有歧义；少于 3 帧报错。 |
| gyrate | -sel all 在多分子周期体系不是良定义量，周期像歧义可达 17%；二维周期共价网络只能给近似并 WARN。 |
| rms | 参考帧必须在窗口内；二维周期网络给近似；少于 2 帧报错。 |
| contact | center 模式用几何质心（非质量加权），与 gmx 可差约 0.01 Å；非整分子跨盒会 WARN。 |
| hbond | 受体集合比 gmx 宽（含卤素/硫），含 Cl⁻ 体系计数可高 +13.8%；选择语义是给体组/受体组，与 gmx 的并集语义不同；ACF 在长轨迹上仍按 O(T²) 量级的滞后循环（已做增量优化）。 |
| agg | 需要残基/分子信息；全连通时 WARN；电荷不平衡只是 WARN。 |
| cluster | 需要残基/分子信息；-focus 不输出尺寸分布。 |
| pistack | 环定义必须显式给出；每个残基至少 3 个环原子；近退化环法向不稳。 |
| freevol | 体素网格近似，不是 Zeo++ 的 Voronoi 几何；网格可能因上限自动变粗；PSD 的 r 是 vdW 表面距离。 |
| sasa | 有限球面采样；需要元素；与 gmx 的残余差异约 0.6%（92 点）。 |
| residence | 轨迹短于 ACF 衰减时间时寿命是下界；-definition 目前只影响表头声明，ACF 始终同时给 inter/cont。 |
| transport | 必须含速度列；NE 忽略离子关联、通常高估电导；缺电荷/温度/frame_interval 时跳过 NE/GK。 |
| occupancy | 不支持按分子体积自动生成 zone；区域越界按交集计体积；少于 2 帧时穿越无意义。 |
| select | resname/name/resid 仅适用于 GROMACS 拓扑；mol 需要分子 ID。 |
| plot1d/plot2d | 依赖 matplotlib；-bar/-stack 互斥且不做平滑；plot2d 不接受 -columns/-bar/-stack。 |
| update | 依赖外部索引可达；内网镜像可能滞后，因此会同时查询多个源并取最高版本；只提示，不自动升级。 |

## 未提供 / 规划

| 名称 | 状态 | 说明 |
|---|---|---|
| order | planned | 取向序参数；参考轴定义未定，未注册，不能使用。 |
| current | planned | 电流分解；目前只有 transport 的 Green–Kubo 总电流。 |
| 自动 zone | planned | occupancy 按分子体积自动生成区域。 |
| -acceptors gmx | planned | 让 hbond 受体集合与 gmx 完全一致的可选开关。 |

## 测试与验收

    pytest tests/unit tests/regression/level2_format

层级：

- 单元/合成：PBC、密度均匀性、随机游走 MSD、刚体 RMSD、手算氢键/SSIP-CIP-AGG/ππ、区域穿越、FFV、键表 id 重排。
- 格式回归：LAMMPS data/lammpstrj 解析、平移不变性、模块冒烟。
- 原生 gmx 对拍：历史脚本未随 1.0.0 发布；关键结论已固化在 CHANGELOG 与模块手册的判读段落。

手册静态校验：tests/unit/test_manuals.py 会重新渲染每个模块的参数表，检查手册存在、help_doc 指向正确、参数集合与默认值一致、声明的输出文件与模块 outputs 一致。

## 维护规则

- 新增/修改参数：同时更新 Param.help、config/lat.yaml、模块手册的参数表（重新生成），然后跑 test_manuals。
- 修改输出文件名/tag：同步更新 cls.outputs、模块手册“输出”表与 CHANGELOG。
- 修改数值口径：必须在 CHANGELOG 记录影响面，并在 MODULE_STATUS 更新已知限制。
- 任何模块级重构后运行 pytest；架构守卫失败时不要绕过，回到 lat.core.periodic。
