#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
ENGINE="$ROOT_DIR/src/stockfish"
VALGRIND_BIN="${VALGRIND:-valgrind}"

if [[ ! -x "$ENGINE" ]]; then
  echo "Error: build the engine at src/stockfish before running this test." >&2
  exit 1
fi

if ! command -v "$VALGRIND_BIN" >/dev/null 2>&1; then
  echo "Warning: valgrind not found; skipping FoW leak check." >&2
  exit 0
fi

run_under_valgrind() {
  local name="$1"; shift
  local script="$1"; shift || true
  local output_file
  output_file="$(mktemp)"

  echo "Running ${name} under valgrind..."
  cat <<SCRIPT | "$VALGRIND_BIN" --error-exitcode=42 --leak-check=full --show-leak-kinds=definite,indirect --quiet "$ENGINE" >"${output_file}" 2>&1
${script}
SCRIPT

  local status=$?
  if [[ $status -ne 0 ]]; then
    echo "[FAIL] ${name}: valgrind reported errors (exit ${status})." >&2
    echo "--- Valgrind/engine output ---" >&2
    cat "${output_file}" >&2
    echo "--------------------------------" >&2
    rm -f "${output_file}"
    exit $status
  fi

  rm -f "${output_file}"
  echo "[PASS] ${name}" 
}

# Short FoW search path.
run_under_valgrind "FoW leak check" "uci
setoption name UCI_Variant value fogofwar
setoption name UCI_FoW value true
setoption name UCI_IISearch value true
setoption name UCI_MinInfosetSize value 64
position fog_fen ????????/??????pp/1?????1P/?1??p1?1/8/1P2P3/PB1P1PP1/NQ1NRBKR b KQk - 0 8
go movetime 10
stop
quit"

# Followed by a non-FoW search to validate teardown.
run_under_valgrind "Baseline leak check" "uci
ucinewgame
position startpos
go movetime 10
stop
quit"

echo "FoW valgrind leak checks passed."
