# tests/compat — 真实项目兼容性战役

`compat_run.py` 浅扫描真实项目目录，把每个 `*.lammpstrj`/`*.dump` 与**同目录且原子数一致**的
拓扑配对，然后对每个体系跑 15–16 个模块（rdf/density/gyrate/rms/msd/contact/occupancy/freevol/
transport/hbond/agg/cluster/residence），记录 PASS/FAIL/SKIP 与错误摘要。

```bash
python tests/compat/compat_run.py \
  --root /mnt/d/CurrentProjects --maxdepth 3 --timeout 120 \
  --out /tmp/compat.md --json /tmp/compat.json
```

注意：F 盘是 USB 机械硬盘，务必只指定**具体项目目录**并限制 `--maxdepth`，不要对 `/mnt/f` 全盘扫描。

要点：
- 按原子数配对，自动跳过 data/dump 不属于同一阶段的情况（GCMC 快照、ETH vs 混合体系）。
- 类型从 dump 首帧采样（`dump_types`），避免拓扑 Masses 列出的类型在 dump 中不存在。
- 拓扑无 `Masses` 时密度用 `-quantity number`；无质量且 dump 无 `element` 列时跳过 freevol/sasa。
- 真实战役结果与修复清单见 `docs/09_test_campaign.md` 第 5b 节。
