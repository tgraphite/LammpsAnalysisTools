# 03 · IO、拓扑、选择器、PBC、单位

## 1. 输入契约

| 文件 | 别名 | 说明 |
|---|---|---|
| `data`（LAMMPS data file） | `-s` / `-topo` | 类型、质量、电荷、键/角、盒信息 |
| `lammpstrj`（dump） | `-f` / `-traj` | 逐帧坐标；必需列 `id type x y z`（或 `xs`） |
| `ndx` | `-index` / `-n` | 选择集；`-sel*` 可为其中组名 |
| `types` | `-types` | 显式 `type→element` 映射 |
| GROMACS（xtc/trr/tpr/gro） | 同 `-s/-f` | **内部支持、不对外文档化**，仅用于回归 |

**不支持多文件/glob 拼接**：一轨迹一模拟，用户自行保证（P0 决策）。

## 2. LAMMPS data 解析

- **data 文件不含 units/atom_style**（都是 input script 命令）。
- `atom_style` 由 `Atoms` 段列数推断（源码 `atom_vec_*.cpp: fields_data_atom`）：

| style | `Atoms` 列序 | 列数 |
|---|---|---|
| atomic | `id type x y z` | 5 |
| charge | `id type q x y z` | 6 |
| molecular / bond / angle | `id molecule type x y z` | 6 |
| full | `id molecule type q x y z` | 7 |

  → 5=atomic、7=full 可确定；**6 列歧义**（charge vs molecular）用启发式（第 2/3 列整型 mol/type vs 浮点 q），
  仍不确定 → **ERROR** 要求 `-style`。
- 读取 `Masses`、`Atoms`、`Bonds/Angles/Dihedrals`、盒（正交或 `xy xz yz` 三斜）。
- 不推断原子名/残基名；`type` 是主选择键。

## 3. lammpstrj 解析

- 解析 `ITEM: TIMESTEP`、`ITEM: BOX BOUNDS`、`ITEM: ATOMS ...`。
- 必需列：`id, type, x y z`（或 `xs ys zs`）。
- 可选列：`mol`、`q`、`ix iy iz`、`vx vy vz`（`transport` 必需）。
- 缺列 → 模块按需 ERROR（如 `transport` 无速度列）。
- **原子顺序按 `id` 对齐**：dump 的原子顺序是 LAMMPS 内部顺序，**不保证**与 data 文件
  （lat 按 `id` 升序排列）一致；逐帧按 `id` 重排到拓扑顺序，否则 `type` 会与错误坐标配对、
  g(r) 完全失真。缺 `id` 列时假定顺序一致并 WARN；`id` 集合与拓扑不一致 → ERROR。
- **时间语义**：只保留帧序号；`timestep` 值忽略（M2）。

## 4. 单位系统

- 默认长度 **Å**（real/metal 皆是）；未显式给 `-unit` → 报 **INFO**，按默认 real/Å。
- 内部统一：长度 Å、时间帧、质量 g/mol、电荷 e、密度 kg/m³。
- **`atom_style` 也不在 data 文件中**（见 §2 表格与推断规则）。
- 未知 unit → ERROR。

## 5. PBC 与三斜盒

- 支持正交与三斜，距离一律 minimum image；`-pbc` 隐藏开关，默认开。
- 分子内排除需用户提供分子定义（`mol` 列或 `Bonds`），**不自动推断**（M3）。

## 6. 坐标包装与解包

- `x y z` 可能 wrapped；`ix iy iz` 为映像计数。
- `rdf/density/contact/...` 用最小镜像，无需解包。
- `msd/transport` 需连续坐标：`xu yu zu` > `ix iy iz` > 最小镜像增量+WARN。
- NPT 变盒下做扩散 → WARN。

## 7. type → element 与 vdW

- 默认按**质量 ±0.2%** 匹配 `element_mass`；多解/无解 → **ERROR 中断**（不采用行尾注释）。
- 质量表来源 **IUPAC 标准原子量（`ase.data`）**，非 reax_tools：
  后者 `Li=7.0` 与真实 6.94 差 0.86%，会超出 ±0.2% 而误报（已在样例中验证）。
- 可用 `-types` 文件或 `lat.yaml:type_element` 显式覆盖。
- vdW 半径表来自 `reax_tools::ELEMENT_ATOMIC_RADII`（Bondi 型），整表在 `lat.yaml:vdw`；
  freevol/sasa 用纯 vdW，键感知才用 `0.6*(ri+rj)`。

## 11. 读取器细节（实现提醒）

- MDAnalysis **不按 `.lammpstrj` 扩展名自动选 reader**，必须显式 `format="LAMMPSDUMP"`。
- `data` 文件自动识别为 `DATA`；`atom_style full` 需能解析 `Atoms # full`。
- 已验证样例：`tests/regression/samples/sys1_reduced/`（坐标偏差 1e-6 Å，元素计数与选择一致）。

## 8. 选择器语义

- 静态 MDAnalysis 语义 + `ndx`；无交互输入（gmx 交互破坏幂等）。
- **LAMMPS 拓扑（data/dump）只有 `type`/`element`/`mol`/`index`/`id`**；`resname`/`name`/`resid`
  仅在 GROMACS 内部读入（回归）时可用。文档/示例不得在 LAMMPS 场景宣传 `resname`（UX 测试发现）。
- `element` 需 `type→element` 可推断；`mol` 需 data 的 `mol` 列或 Bonds。
- 读入 `-index`/`-n` 时，`-sel1/-sel2` 等**允许是 ndx 中的组名**；否则按 MDAnalysis 选择解析。
- 选择子支持 `and/or/not` 与括号。
- **selection 字面值只写入输出文件头**（可能为 `file.ndx -> water` 或 `type 1 2 3`）。
- `lat select`：打印命中数/type 分布/mass 映射；`-o out.ndx` 导出校验。

## 9. 自动网格

- 一维直方图 `nbins=200`；三维体素 `dr=1.0 Å`，单轴上限 500 bins。
- **绝对轴分布按 gmx 口径分数分箱**（NPT 盒变安全），轴按平均盒长标定。
- `smooth: none|savgol|gaussian` 只在 **plot1d/plot2d** 生效；计算命令始终输出原始值。

## 10. 内存与流式（M9）

- 永远逐帧访问，不缓存整条轨迹；模块声明复杂度与峰值内存；支持 `-stride`。
