#!/usr/bin/env bash
# Benchmark/comparison harness for FoW vs baseline play
set -euo pipefail
cd "$(dirname "$0")/.."

if [[ ! -x src/stockfish ]]; then
  (cd src && make build ARCH=x86-64-modern)
fi

# Baseline chess short bench
src/stockfish bench 1 1 1 default depth 9 nodes 200000 || true

# FoW micro-benchmark parity
src/stockfish --fow-benchmarks || true
