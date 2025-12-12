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

# Castling rights must remain visible even when the board view is partially hidden via fog_fen.
run_case "FoW castling rights preserved" "Fen:.*KQkq" "uci
setoption name UCI_Variant value fogofwar
setoption name UCI_FoW value true
setoption name UCI_IISearch value true
position fog_fen r3k2r/8/8/8/8/8/8/R3K2R w KQkq - 0 1
d
quit"

# En-passant targets should be announced even when the capturing pawn is hidden behind fog.
run_case "FoW en-passant marker visible" "En passant: e6" "uci
setoption name UCI_Variant value fogofwar
setoption name UCI_FoW value true
position fog_fen rnbqkbnr/ppp1pppp/8/3pP3/8/8/PPPP1PPP/RNBQKBNR b KQkq e6 0 2

d
quit"

# Crazyhouse pockets must be exposed in darkcrazyhouse while FoW is active.
run_case "FoW crazyhouse pockets visible" "\\[.*\\]" "uci
setoption name UCI_Variant value darkcrazyhouse
setoption name UCI_FoW value true
position startpos

d
quit"

echo "All FoW API visibility checks passed."
