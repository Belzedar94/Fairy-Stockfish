#!/usr/bin/env bash
# Dedicated FoW C++ unit suites (visibility, belief, CFR, selection, KLUSS)
set -euo pipefail
cd "$(dirname "$0")/.."

if [[ ! -x src/stockfish ]]; then
  (cd src && make build ARCH=x86-64-modern)
fi

src/stockfish --fow-unittests
