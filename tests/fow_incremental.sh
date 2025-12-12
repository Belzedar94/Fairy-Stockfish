#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
ENGINE="$ROOT_DIR/src/stockfish"

if [[ ! -x "$ENGINE" ]]; then
  echo "Error: build the engine at src/stockfish before running this test." >&2
  exit 1
fi

# Simple integration check for fog-of-war incremental belief handling
OUTPUT_FILE="$(mktemp)"

cat <<'SCRIPT' | "$ENGINE" >"$OUTPUT_FILE" 2>&1
uci
setoption name UCI_Variant value fogofwar
setoption name UCI_FoW value true
setoption name UCI_IISearch value true
position fog_fen ????????/??????pp/1?????1P/?1??p1?1/8/1P2P3/PB1P1PP1/NQ1NRBKR b KQk - 0 8
go movetime 50
stop
quit
SCRIPT

grep -q "info string FoW search" "$OUTPUT_FILE"
rm -f "$OUTPUT_FILE"

echo "FoW incremental belief test passed."
