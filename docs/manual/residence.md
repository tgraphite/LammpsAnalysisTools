# residence — 配体滞留时间

> 状态：稳定（1.0.0） · 需要拓扑：是 · 需要轨迹：是（≥2 帧）

## 用途

统计配体在参考组周围的结合/解离动力学：结合数、间断与连续自相关、寿命和连续结合区间分布。

## 用法

    lat residence -s system.data -f run.lammpstrj -ref "type 2" -lig "type 5" -cutoff 3.5 -o out

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-ref` | 字符串 | 必填 | reference group |
| `-lig` | 字符串 | 必填 | ligand group |
| `-cutoff` | 浮点数 | 必填 | binding cutoff in Angstrom (required; no default) |
| `-mode` | 字符串 | min | distance criterion: min atom-atom or residue COM distance |
| `-definition` | 字符串 | both | lifetime definition(s) to report |
| `-frame_range` | 字符串 | — | frame window start:stop[:stride] |
| `-stride` | 整数 | — | sample every N frames (default: 1) |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| <base>.residence.out.csv | frame 与结合配体数 |
| <base>.residence.acf.out.csv | lag、间断自相关 C_inter、连续自相关 C_cont |
| <base>.residence.intervals.out.csv | 连续结合区间长度（帧）与计数 |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

把每个配体分子视为一个单元：到参考组的距离（-mode min 用最近原子对，-mode center 用残基质心）小于 cutoff 记 h=1。间断自相关允许离开后再回来；连续自相关要求整个区间持续结合。寿命同时给积分法 τ_int=Σ_{τ≥1}C 和指数拟合 τ_fit=−1/斜率，单位是帧；配置 meta.frame_interval 后额外给出以 ps 表示的标量。区间分布统计每个配体连续 h=1 的游程长度。

## 判读

结合数应与体系配位情况一致；ACF 应从 1 单调衰减；τ_int 与 τ_fit 应同量级。若轨迹长度不足半个衰减时间，寿命会被截断并 WARN。

## 告警与限制

- cutoff 必填。
- 没有配体在 cutoff 内结合时会 WARN。
- 轨迹短于自相关衰减时间时寿命是下界，不能直接外推。
- -definition 目前只影响表头声明；ACF 表始终同时给出 C_inter 与 C_cont。

## 相关

氢键寿命见 [hbond.md](hbond.md)。
