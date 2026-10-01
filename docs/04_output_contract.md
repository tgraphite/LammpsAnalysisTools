# 04 · 输出契约

## 1. 原则

- **一文件一表**：允许同一 X 轴下多个字段（多条曲线）。
- **纯 CSV + `#` 注释头**：人类可读、机器可解析。
- **模块不写文件**：模块返回 `Table`，由 `output.py` 统一序列化。
- **原子写**：tmp + rename。

## 2. 命名（命名即意图）

```
<basename>.<tag>.<mode>.out.csv
```

- `basename`：用户给的系统名（`-o` 的基础名），多表共享。
- `tag`：模块与子项，如 `rdf`、`rdf.cn`、`rg`。
- `mode`：可选模式，如 `cont` / `inter`。
- 例：`mysystem.rdf.cn.out.csv`、`mysystem.rg.out.csv`、`mysystem.rdf_Li_FSI.out.csv`。
- **selection 字面值只存在于文件头**，不进文件名（其可能为 `file.ndx -> water` 或 `type 1 2 3`）。

## 3. 同名冲突

- 不覆盖、不报错：已存在的目标文件依次改名为 `*.bak.1`、`*.bak.2`…（`io.backup_on_conflict`）。
- 事件在审计中记录。

## 4. Table 对象

```python
Table(
  title   = "RDF Li-FSI",
  x       = Field(name="r", unit="Angstrom"),
  y       = [Field(name="g_Li_FSI", unit="-"), ...],
  scalars = {"r_first_min/A": ..., "cn_first_shell": ...},
  warnings= ["..."],
  meta    = {"module": "rdf", "seed": ...},
)
```

## 5. CSV 头（已冻结）

```
# program: LAT version 0.1.0
# module: rdf
# command: lat rdf -s sys.data -f sys.lammpstrj -sel1 type 1 -sel2 type 2
#
# selections: sel1=type 1; sel2=type 2
# definitions: cutoff=10 A; nbins=200
# frame_range: 0:100:1
# scalars: r_first_min/A=2.10; cn_first_shell=4.02
# warning: none
#
# title: RDF
# xlabel: r (Angstrom)
# ylabel: g(r)
#
# columns: r, g_Li_FSI
# --- data ---
2.000000,0.000000
...
```

- **绘图文本契约**：`title` / `xlabel` / `ylabel` / `columns` 四个字段供 `plot1d` 直接出图；
  **`columns` 的 y 列名即 legend**（用户改成 `Water-Oil` 之类即可得到明白的图例）。
- `selections`：与意图相关的选择输入（sel1/sel2/sel…）；`definitions`：**影响结果的实际使用的参数**
  （cutoff、angle、radius、nbins、temperature、mode…）。
- `scalars` / `warning` 保留在头部（WARN 同时屏显，M12）。
- 审计信息（`program`/`command`）默认写入文件（供复现），但不主动屏显（M5）。
- `-debug` 时额外落一份完整 JSON 审计（含开发统计字段）；**对外文档不暴露 debug 模式**。

## 6. 数据形态

| 形态 | 写法 |
|---|---|
| 1D 曲线 | `x, y1, y2, ...` |
| 2D 网格 | 长表 `x, y, value` |
| 时间序列 | `frame, y1, ...` |
| 标量 | 写 header `scalars`，不单独成表 |

## 7. 日志与退出

- 不用退出码；分级 **INFO / WARNING / ERROR**（ERROR 中断）。
- 全部输出到 **stderr**，绝不污染 stdout / CSV。
- `WARN` 同时进 header 与屏显；只提示，不修改结果（M12）。

## 8. 绘图

- 输入用 `-f`：`lat plot1d -f <csv> [-o out.png] [-title ...]`；**不再支持位置参数**。
- 只读 CSV，从 header 取 `title`/`xlabel`/`ylabel`；`columns` 的 y 列名作 legend，legend 默认**无边框**。
- **`-columns NAME[,NAME...]`**：只画给定名称的 y 列（逗号或空格分隔，保留用户给定顺序、自动去重）；
  **x 轴恒为第一列且总是绘制**，第一列名不可作为 y 选择（只选到第一列会报错；未知名报错并列出可用列）。
- **`-bar`**：分组柱状图，每个 y 列一根柱，按组并排（柱宽按 x 间距自适应）。
- **`-stack`**：堆积柱状图，**逐 x 把各 y 列归一化成占据比**（和为 1，零和组记 0），y 轴固定 0–1。
- `-bar` 与 `-stack` **互斥**；两者都只属于 `plot1d`；柱状模式**不做平滑**（`-smooth` 非 none 时报错）。
- 单文件单图；默认输出 `<csv stem>.png`。`plot2d` 同理，2D 用 `MaxNLocator` 配色条。
- 计算/绘图解耦（M14）。
