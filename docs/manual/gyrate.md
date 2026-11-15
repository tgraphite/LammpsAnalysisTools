# gyrate — 回转半径 Rg

> 状态：稳定（1.0.0） · 需要拓扑：是 · 需要轨迹：是（≥1 帧）

## 用途

逐帧计算选区分子的回转半径 Rg，用于表征聚合物链尺寸、溶胀与折叠。

## 用法

    lat gyrate -s system.data -f run.lammpstrj -sel "mol 1" -o out

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-sel` | 字符串 | 必填 | atoms making up the molecule/chain (use one molecule per run) |
| `-frame_range` | 字符串 | — | frame window start:stop[:stride] |
| `-stride` | 整数 | — | sample every N frames (default: 1) |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| <base>.gyrate.out.csv | frame 与 Rg（Å） |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

Rg=sqrt(Σ m_i|r_i−r_cm|²/Σ m_i)。计算前按键连接做整分子复原（等价 gmx rmpbc），否则跨盒聚合物会被系统性低估。无质量时使用等权并 WARN。键表按 atom id 解释，data 文件 Atoms 段不按 id 排序时也会正确重排。

## 判读

链状分子 Rg 应在 nm 量级并随链长/溶剂条件变化；跨盒大分子不应异常缩小。多分子体系请按分子分别计算，-sel all 的均值不是良定义量：每个分子可以放在任意周期像里，分子间相对位置依赖锚点约定，实测歧义可达 17%。

## 告警与限制

- 选择不是整分子时会 WARN。
- 二维周期共价网络（电极片层、凝胶）无法整分子复原，会 WARN 并给出近似值；gmx 在等价场景直接报错。
- 没有键或分子信息时跨盒分子可能被低估。

## 相关

RMSD/RMSF 见 [rms.md](rms.md)。
