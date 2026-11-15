# select — 选择检查与 .ndx 导出

> 状态：稳定（1.0.0） · 需要拓扑：是 · 需要轨迹：否（可选）

## 用途

检查选择表达式命中的原子、类型分布、电荷、分子数和残基名，并可导出 GROMACS .ndx。

## 用法

    lat select -s system.data -sel "type 2"
lat select -s system.data -sel "type 2" -export-ndx group.ndx -name LIG

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-sel` | 字符串 | — | selection expression (omit to list all types) |
| `-export-ndx` | 字符串 | — | write the selection to this .ndx |
| `-name` | 字符串 | — | group name in the exported .ndx |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| （屏幕） | 原子数、类型分布、电荷和、分子数、残基名 |
| <name>.ndx (optional) | 可选：导出为 GROMACS 索引组 |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

选择表达式使用 MDAnalysis 风格子集：all、type、element、name、resname、resid、mol、index、id，支持 and/or/not 与括号；相邻关键字隐含 and。LAMMPS 拓扑通常只有 type/element/mol/id/index；resname/name/resid 需要 GROMACS 拓扑。导出时原子编号按 GROMACS 的 1-based 规则写入。

## 判读

先确认命中数、类型分布和电荷和是否符合预期，再把同一表达式用于分析模块。选错关键字通常会得到 0 原子并报错。

## 告警与限制

- 不提供轨迹也能工作；给了 -f 也不会读帧。
- resname/name/resid 只适用于带这些字段的拓扑。
- 选择 'mol' 需要拓扑中的分子 ID。
- 导出 .ndx 的组名默认是 lat。

## 相关

选择语法汇总见 [common.md](common.md#选择语法)。
