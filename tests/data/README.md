# 自造小体系（用于无 GT 的模块）

- `hbond/`：构造含 A–H···B 的少量分子，手算氢键数（只验数值一致）。
- `residence/`：一个中心 + 若干配体，构造已知停留区间。
- `occupancy/`：两/三区域 + 已知穿越次数的粒子。
- `freevol/`、`sasa/`：均相液体（FFV 基线）、人造多孔网格（PSD/连通性）、表面可及面积。
- `transport/`：带 `vx vy vz` 的小 dump。

每个体系附 `expected.md` 记录手算期望值与口径。
