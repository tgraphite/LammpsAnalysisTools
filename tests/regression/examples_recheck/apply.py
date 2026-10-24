#!/usr/bin/env python3
"""Apply the rechecked results back into the example directories.

Nothing is touched until you pass `--apply`: `recheck.py` writes only into
`_lat_v020_recheck/`, and MAPPING.tsv records what would go where.

    python3 apply.py                  # dry run: print every action (mutates nothing)
    python3 apply.py --apply          # copy regenerated files over the originals
    python3 apply.py --apply --delete # ... and remove invalid/scaffolding files
    python3 apply.py --apply --restore  # undo: copy backup_before_overwrite/ back

Every overwritten or deleted file is first copied to
`_lat_v020_recheck/backup_before_overwrite/<original relative path>`.

Deletion needs **both** `--apply` and `--delete`; a dry run never writes or
removes anything (an earlier version of this script deleted during a dry run).
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXAMPLES = HERE.parent
BACKUP = HERE / "backup_before_overwrite"

# Round-1 scratch that lives in the examples tree: deleted only with --delete.
SCAFFOLD_DIRS = [
    "31_linear_polymer/analysis/lat_stage",
]


def rows():
    with open(HERE / "MAPPING.tsv") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def backup(path: Path) -> Path:
    dest = BACKUP / path.relative_to(EXAMPLES)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        if path.is_dir():
            shutil.copytree(path, dest)
        else:
            shutil.copy2(path, dest)
    return dest


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true", help="actually write/remove")
    ap.add_argument("--delete", action="store_true",
                    help="also remove files whose status is invalid/scaffolding")
    ap.add_argument("--restore", action="store_true",
                    help="copy backup_before_overwrite/ back over the examples")
    a = ap.parse_args()
    write = a.apply and not a.restore

    if a.restore:
        n = 0
        for src in sorted(BACKUP.rglob("*")):
            if src.is_file():
                dst = EXAMPLES / src.relative_to(BACKUP)
                print(f"restore {dst.relative_to(EXAMPLES)}")
                if write:
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(src, dst)
                n += 1
        print(f"{'restored' if write else 'would restore'} {n} file(s)"
              + ("" if write else "  (dry run -- add --apply)"))
        return 0

    n_copy = n_del = n_skip = 0
    for r in rows():
        status, old_rel, new_rel = r["status"], r["old_file"], r["new_file"]
        if status in ("affected", "verify") and old_rel and new_rel:
            src, dst = HERE / new_rel, EXAMPLES / old_rel
            if not src.is_file():
                print(f"SKIP (missing mirror file): {new_rel}")
                n_skip += 1
                continue
            print(f"copy   {new_rel}\n    -> {old_rel}")
            if write:
                if dst.exists():
                    backup(dst)
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)
            n_copy += 1

    if a.delete:
        to_remove = [r["old_file"] for r in rows()
                     if r["status"] in ("invalid", "scaffolding") and r["old_file"]]
        to_remove += SCAFFOLD_DIRS
        for rel in to_remove:
            dst = EXAMPLES / rel
            if not dst.exists():
                continue
            what = "dir " if dst.is_dir() else "file"
            print(f"delete {what} {rel}")
            if write:
                backup(dst)
                if dst.is_dir():
                    shutil.rmtree(dst)
                else:
                    dst.unlink()
            n_del += 1
    else:
        n_pending = sum(1 for r in rows()
                        if r["status"] in ("invalid", "scaffolding") and r["old_file"]
                        and (EXAMPLES / r["old_file"]).exists())
        n_pending += sum(1 for d in SCAFFOLD_DIRS if (EXAMPLES / d).exists())
        if n_pending:
            print(f"({n_pending} invalid/scaffolding path(s) still present; "
                  f"add --delete to remove them)")

    verb = "applied" if write else "would apply"
    print(f"\n{verb}: {n_copy} copy, {n_del} delete"
          + (f", {n_skip} skipped" if n_skip else "")
          + ("" if write else "  (dry run -- add --apply to write)"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
