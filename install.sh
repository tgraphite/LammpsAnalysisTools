#!/usr/bin/env bash
# One-click installer for lat (LAMMPS Analysis Tool).
#
#   bash install.sh
#   LAT_VENV=/opt/lat-venv bash install.sh
#   PYTHON=/usr/bin/python3.11 bash install.sh
#
# Works on machines whose system Python is older than 3.9: it bootstraps a
# private interpreter via `uv` (single static binary, no root needed).
# Afterwards run `lat <module> -s <data> -f <traj> ...` from anywhere.
set -euo pipefail

export PATH="$HOME/.local/bin:$PATH"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="${LAT_VENV:-$ROOT/.venv}"
BIN="${LAT_BIN:-$HOME/.local/bin}"

find_python() {
    if [[ -n "${PYTHON:-}" ]]; then echo "$PYTHON"; return; fi
    for p in python3 python3.13 python3.12 python3.11 python3.10 python3.9; do
        if command -v "$p" >/dev/null 2>&1; then
            if "$p" -c 'import sys; raise SystemExit(0 if sys.version_info>=(3,9) else 1)' 2>/dev/null; then
                echo "$p"; return
            fi
        fi
    done
}

PY="$(find_python || true)"
echo "[lat] repo    : $ROOT"
echo "[lat] venv    : $VENV"
echo "[lat] bin dir : $BIN"

if [[ -n "$PY" ]]; then
    echo "[lat] using system interpreter: $($PY --version 2>&1) ($(command -v "$PY"))"
    "$PY" -m venv --clear "$VENV"
    "$VENV/bin/python" -m pip install --upgrade pip wheel >/dev/null
    echo "[lat] installing lat + dependencies (numpy scipy pyyaml matplotlib MDAnalysis ase)..."
    "$VENV/bin/python" -m pip install "$ROOT"
else
    echo "[lat] no Python >= 3.9 found; bootstrapping a private one with uv..."
    if ! command -v uv >/dev/null 2>&1; then
        curl -LsSf https://astral.sh/uv/install.sh | sh
        export PATH="$HOME/.local/bin:$PATH"
    fi
    uv venv "$VENV" --python 3.11 --clear
    echo "[lat] installing lat + dependencies via uv..."
    uv pip install --python "$VENV/bin/python" "$ROOT"
fi

mkdir -p "$BIN"
ln -sf "$VENV/bin/lat" "$BIN/lat"

echo
echo "[lat] installed: $("$BIN/lat" --version)"
if ! echo ":$PATH:" | grep -q ":$BIN:"; then
    echo "[lat] NOTE: $BIN is not on your PATH. Add to your shell profile:"
    echo "      export PATH=\"$BIN:\$PATH\""
fi
echo "[lat] try:  lat --help"
