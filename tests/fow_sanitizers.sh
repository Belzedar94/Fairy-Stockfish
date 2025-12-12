#!/usr/bin/env bash
# ASan/TSan-backed FoW init/teardown probes
set -euo pipefail
cd "$(dirname "$0")/.."

build_with_sanitize() {
  local mode=$1
  (cd src && make build ARCH=x86-64-modern debug=yes sanitize=$mode EXE="stockfish-${mode}")
}

for mode in address thread; do
  if [[ ! -x src/stockfish-${mode} ]]; then
    build_with_sanitize "$mode"
  fi
  src/stockfish-${mode} --fow-sanitizers || true
  src/stockfish-${mode} --fow-unittests || true
done
