# rms — RMSD 与 RMSF

> 状态：稳定（1.0.0） · 需要拓扑：是 · 需要轨迹：是（≥2 帧）

## 用途

用 Kabsch 最优旋转对齐后，逐帧计算 RMSD，并逐原子计算 RMSF。

## 用法

    lat rms -s system.data -f run.lammpstrj -sel "type 1" -fit "type 1" -ref first -o out

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-sel` | 字符串 | 必填 | atoms used for RMSD/RMSF |
| `-fit` | 字符串 | — | atoms used for alignment (default: -sel) |
| `-ref` | 字符串 | first | first, average, or an explicit frame index |
| `-rmsf` | 开关 | true | also write per-atom RMSF |
| `-frame_range` | 字符串 | — | frame window start:stop[:stride] |
| `-stride` | 整数 | — | sample every N frames (default: 1) |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| <base>.rmsd.out.csv | frame 与 RMSD |
| <base>.rmsf.out.csv | 原子 id 与 RMSF |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

先做整分子复原，再用 fit 组求最优旋转（SVD，含镜像修正），在 sel 组上计算 RMSD。参考结构可以是首帧、指定帧号或逐帧坐标的算术平均（-ref average）。RMSF 是对齐后各帧坐标相对其时间平均位置的涨落，逐原子输出。

## 判读

RMSD 用于判断整体构象漂移或平衡；RMSF 用于定位柔性区域。对比不同轨迹时保持 sel/fit 一致。

## 告警与限制

- 少于 2 帧直接报错。
- 参考帧必须在所选 frame_range 内。
- 跨盒且非整分子会 WARN；二维周期网络会给出近似值。
- 未请求 RMSF（-rmsf 关闭）时只写 RMSD 表并 WARN。

## 相关

回转半径见 [gyrate.md](gyrate.md)。
