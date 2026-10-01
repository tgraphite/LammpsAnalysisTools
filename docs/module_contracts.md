# 模块总表（输入 / 输出 / 警报）

> 命令行统一：`lat <module> -s <data> -f <traj> [...]`；`-s`=`-topo`，`-f`=`-traj`（单横线）。
> `frame_range` 一律为帧区间（可 `start:stop:stride`）。默认长度单位 Å。所有输出走统一 `Table` 契约（见 04）。
> 输出命名：`<basename>.<tag>.<mode>.out.csv`（命名即意图）。
> 除标注 **(P2)** 外，下表全部为 **P0**。详细算法与原理见 `docs/MANUAL_ALGORITHMS.md`。

| 模块 | 必选 | 可选 | 输出表 | 主要警报 |
|---|---|---|---|---|
| `rdf` | `-s` `-sel1` `-sel2` | `-f` `-nbins` `-cutoff` `-frame_range` `-intramol` `-pbc` | g(r) 表；CN 表 | 选择为空；帧不足；盒<2·cutoff；有键但无分子 ID 无法排除同分子对（默认排除，`-intramol` 包含） |
| `density` | `-s` `-sel` `-axis` | `-f` `-nbins` `-dim` `-frame_range` `-quantity` | 1D 曲线（X, ρ）；2D 长表（x,y,value） | 轴长/bins 退化；盒长偏离均值 >±2% |
| `hbond` | `-s` | `-f` `-sel1` `-sel2` `-distance`(默认3.5Å D-A) `-distance_mode` `-angle`(30°) `-nbins` `-frame_range` | 数量-帧；ACF+寿命(标量,单位帧)；3 个面分布 | 无给体；两组既不重合也不独立；帧少拟合不出 |
| `rms` | `-s` `-sel` | `-f` `-fit` `-ref` | RMSD-帧；RMSF-原子 | PBC 跳跃；选择非整分子 |
| `msd` | `-s` `-sel` | `-f` `-dim` `-fit_range` `-stride` | MSD-帧；标量 D、R²、拟合窗 | 帧少；R² 低；COM 漂移未去 |
| `agg` | `-s` `-sel1` `-sel2` `-cutoff`(必填) | `-f` `-mode`(min/center) `-solv` `-frame_range` | 团簇类型计数-帧；标签占比；AiBjCk 表 | 电荷不平衡；边数守恒偏差 |
| `cluster` | `-s` `-sel1` `-sel2` `-cutoff`(必填) | `-f` `-mode` `-focus` `-frame_range` | 配位数分布-帧；团簇尺寸 | 同 agg |
| `pistack` | `-s` `-sel` | `-f` `-distance` `-angle` `-frame_range` | 数量-帧；逐对几何 | 未检出芳环；需显式环定义 |
| `probe` | `-s` | `-f` `-radius`(默认1.445Å=H2) `-mode` `-grid_dr` `-sel` | FFV-帧 / SASA-帧 / PSD 表 / 孔隙率 / 连通分量 | 元素 vdW 缺失；网格过粗；周期连通未闭合 |
| `residence` | `-s` `-ref` `-lig` | `-f` `-cutoff` `-mode`(cont/inter) `-frame_range` | 存在函数 ACF；寿命(帧)；区间分布；占据数-帧 | 无配体命中；帧少；统计窗截断偏差 |
| `transport` | `-s` `-sel` | `-f` `-q` `-frame_range` | D 表；电导表；迁移数 | 缺速度列(ERROR)；缺电荷；帧少 |
| `occupancy` | `-s` `-sel` `-zone` | `-f` `-frame_range` | 区域占据-帧；穿越计数；首次通过分布 | 区域重叠/越界；无穿越 |
| `gyrate` | `-s` `-sel` | `-f` `-frame_range` | Rg-帧；标量均值/分布 | 选择非整分子 |
| `contact` | `-s` `-sel1` `-sel2` | `-f` `-cutoff`(必填) `-mode` `-frame_range` | 接触数-帧；最小/质心距离分布 | 同 rdf |
| `order` **(P2)** | `-s` `-sel` | `-f` `-axis` `-frame_range` | 序参数-帧；角度/二面角分布 | 取向无法定义 |
| `select` | `-s` `-sel` | `-f` `-o out.ndx` | ndx（可选）；屏显统计 | 选择为空 |
| `plot1d` | `-f`(CSV) | `-o` `-title` `-smooth` `-columns` `-bar` `-stack` | 图文件 | header 缺失；未知名；`-bar`/`-stack` 互斥且仅限 plot1d |
| `plot2d` | `-f`(CSV) | `-o` `-title` `-smooth` | 图文件 | header 缺失 |
| `describe` | `-module` 或全部 | — | 模块清单/文档 | — |

## 依赖关系（算法内核 → 模块）

| 内核 | 被哪些模块使用 |
|---|---|
| 近邻/距离 (K1) | rdf, hbond, agg, cluster, pistack, contact, residence, probe, order, msd |
| 直方图/网格 (K2/K9) | rdf, density, probe, hbond(分布), order |
| 时间相关 (K3/K4) | msd, hbond(ACF), residence, transport |
| 图/连通分量 (K5) | agg, cluster, probe(connectivity), pistack(网络可选) |
| 拟合/降维 (K6) | rms |
| 压力/速度 (K8) | transport |
| 区域/时间轴 (event-zone) | occupancy |
