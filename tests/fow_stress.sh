#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
ENGINE="$ROOT_DIR/src/stockfish"

if [[ ! -x "$ENGINE" ]]; then
  echo "Error: build the engine at src/stockfish before running this test." >&2
  exit 1
fi

run_stress_cycle() {
  local cycle=$1
  local output_file
  output_file="$(mktemp)"

  echo "Running FoW stress cycle ${cycle}..."
  cat <<SCRIPT | "$ENGINE" >"${output_file}" 2>&1
uci
setoption name UCI_Variant value fogofwar
setoption name UCI_FoW value true
setoption name UCI_IISearch value true
setoption name UCI_MinInfosetSize value 64
position fog_fen ????????/??????pp/1?????1P/?1??p1?1/8/1P2P3/PB1P1PP1/NQ1NRBKR b KQk - 0 8
go movetime 25
position startpos
setoption name UCI_FoW value false
position fen r1bqkbnr/pp1ppppp/2n5/2p5/3P4/5N2/PPP1PPPP/RNBQKB1R w KQkq - 2 3
go depth 6
quit
SCRIPT

  if ! grep -Eq "bestmove" "${output_file}"; then
    echo "[FAIL] FoW stress cycle ${cycle}: bestmove missing" >&2
    echo "--- Engine output ---" >&2
    cat "${output_file}" >&2
    echo "---------------------" >&2
    rm -f "${output_file}"
    exit 1
  fi

  rm -f "${output_file}"
  echo "[PASS] FoW stress cycle ${cycle}"
}

for cycle in $(seq 1 5); do
  run_stress_cycle "$cycle"
done

echo "FoW stress/regression cycles completed successfully."
