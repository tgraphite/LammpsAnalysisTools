# pistack — π–π 堆积

> 状态：稳定（1.0.0） · 需要拓扑：是 · 需要轨迹：是（≥1 帧）

## 用途

以残基为单位定义芳环，用环心距离和环法向夹角判定 π–π 堆积，输出逐帧计数与逐对几何。

## 用法

    lat pistack -s system.data -f run.lammpstrj -sel1 "type 3" -o out
lat pistack -s system.data -f run.lammpstrj -sel1 "type 3" -sel2 "type 4" -distance 5.5 -angle 30 -o out

## 参数

<!-- lat:params:begin -->
| 参数 | 类型 | 默认 | 说明 |
|---|---|---|---|
| `-sel1` | 字符串 | 必填 | ring atoms of group A |
| `-sel2` | 字符串 | — | ring atoms of group B (default: self-stacking) |
| `-distance` | 浮点数 | 5.5 | ring-centre cutoff in Angstrom |
| `-angle` | 浮点数 | 30 | maximum angle between ring normals in degrees |
| `-frame_range` | 字符串 | — | frame window start:stop[:stride] |
| `-stride` | 整数 | — | sample every N frames (default: 1) |
<!-- lat:params:end -->

公共参数（-s、-f、-o、-frame_range、-stride 等）见 [common.md](common.md)。

## 输出

| 文件 | 内容 |
|---|---|
| <base>.pistack.out.csv | frame 与堆积对数 |
| <base>.pistack.pairs.out.csv | frame、环编号、环心距离、法向夹角 |

所有 CSV 都带统一文件头；标量结果在 # scalars: 行，见 [common.md](common.md#输出约定)。

## 原理

每个残基在 sel1（以及可选的 sel2）中贡献一个环：先用键连接整分子复原，再取环原子质心与 SVD 最小奇异向量作为法向。两个环判定为堆积需同时满足环心距离 < distance（默认 5.5 Å）和法向夹角 < angle（默认 30°）。sel2 省略时做自堆积。

## 判读

逐帧计数应与芳香堆积的预期一致；pairs 表可用于检查具体是哪两个环、距离和角度。对同一体系改变 distance/angle 时计数应平滑变化。

## 告警与限制

- sel1 中每个残基至少要有 3 个环原子，否则报错。
- 环定义必须由用户显式给出；lat 不猜测芳香性。
- 环法向由 SVD 得到，近共线/退化环会不稳定。
- 没有找到任何堆积对时会 WARN。

## 相关

氢键见 [hbond.md](hbond.md)。
