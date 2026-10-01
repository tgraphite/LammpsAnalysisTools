# level1_algo — lat vs GROMACS 逐值回归

`compare_gmx.py` 在运行时用 `gmx` 生成 Ground Truth，再用 lat 的 Python API（拓扑只加载一次）
逐对比较 RDF 与配位数。

```
python compare_gmx.py --system sys1 \
  --base /mnt/f/CRO/40_multi_electrolytes/11-systems \
  --rmax 2.0 --bin 0.002 --pairs LI:LI,LI:FSI,LI:FEC,LI:FEMC --frames 90:101:1
```

- 真实体系：sys1–sys6（178k–188k 原子，sys4/sys6 含 PMBA 聚合物）。
- 参数：0.002 nm 网格，rmax 2 nm，帧 9000–10000 ps（11 帧）。
- 指标：`max|Δg|`、`rms`、`mean`、峰位/峰高、`cn_max|d|`、`cn_rel`。

## 结果（2026-09-18）

| 体系 | max\|Δg\| | rms | cn_rel | 备注 |
|---|---|---|---|---|
| sys1 | 0.134 | 0.021 | 0.16% | 纯电解液 |
| sys2 | 0.137 | 0.022 | 0.21% | +TTE |
| sys3 | 0.155 | 0.020 | 0.16% | +MBA |
| sys4 | 0.151 | 0.022 | 0.16% | **+PMBA 聚合物** |
| sys5 | 0.151 | 0.022 | 0.20% | MBA+TTE |
| sys6 | 0.184 | 0.024 | 0.20% | **PMBA+TTE** |

首峰高度可达 g≈12.7，故 max\|Δg\|≤0.184 ≈ **首峰 <1.5%**；mean 偏置 ~1e-5。结论：RDF/CN 与 `gmx rdf`
在分子混合物 + 大分子体系上事实上一致。密度剖面 <1%、MSD（mol_com）曲线 ~3–7%。

详见 `docs/09_test_campaign.md`。
