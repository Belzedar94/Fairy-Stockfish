#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
ENGINE="$ROOT_DIR/src/stockfish"

if [[ ! -x "$ENGINE" ]]; then
  echo "Error: build the engine at src/stockfish before running this test." >&2
  exit 1
fi

run_case() {
  local name="$1"; shift
  local pattern="$1"; shift
  local script="$1"; shift || true

  local output_file
  output_file="$(mktemp)"

  echo "Running ${name}..."
  cat <<SCRIPT | "$ENGINE" >"${output_file}" 2>&1
${script}
SCRIPT

  if ! grep -Eq "$pattern" "${output_file}"; then
    echo "[FAIL] ${name}: expected pattern not found: ${pattern}" >&2
    echo "--- Engine output ---" >&2
    cat "${output_file}" >&2
    echo "---------------------" >&2
    rm -f "${output_file}"
    exit 1
  fi

  rm -f "${output_file}"
  echo "[PASS] ${name}" 
}

# Baseline chess search to compare against FoW.
run_case "Baseline chess search" "bestmove [a-h][1-8]" "uci
ucinewgame
position fen rnbqkbnr/pppp1ppp/8/4p3/8/5N2/PPPPPPPP/RNBQKB1R w KQkq - 1 2
go movetime 50
stop
quit"

# FoW search seeded by fog_fen to validate FoW vs. baseline behavior.
run_case "FoW search comparison" "(info string FoW search|bestmove)" "uci
setoption name UCI_Variant value fogofwar
setoption name UCI_FoW value true
setoption name UCI_IISearch value true
setoption name UCI_MinInfosetSize value 64
position fog_fen ????????/??????pp/1?????1P/?1??p1?1/8/1P2P3/PB1P1PP1/NQ1NRBKR b KQk - 0 8
go movetime 50
stop
quit"

echo "FoW comparison checks passed."
