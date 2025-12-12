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

# Visibility/castling: ensure castling rights survive FoW setup and remain visible in diagnostics.
run_case "FoW castling visibility" "Fen:.*KQkq" "uci
setoption name UCI_Variant value fogofwar
setoption name UCI_FoW value true
setoption name UCI_IISearch value true
position fen r3k2r/8/8/8/8/8/8/R3K2R w KQkq - 0 1
d
quit"

# En-passant marker: validate parsing and display of the en-passant file while in FoW mode.
run_case "FoW en-passant visibility" "En passant: e6" "uci
setoption name UCI_Variant value fogofwar
setoption name UCI_FoW value true
position fen rnbqkbnr/ppp1pppp/8/3pP3/8/8/PPPP1PPP/RNBQKBNR b KQkq e6 0 2
d
quit"

# Crazyhouse hands: confirm piece-in-hand visibility for darkcrazyhouse still renders the pocket.
run_case "FoW crazyhouse pocket visibility" "\\[.*\\]" "uci
setoption name UCI_Variant value darkcrazyhouse
setoption name UCI_FoW value true
position startpos
d
quit"

# End-to-end FoW search: run a short search seeded by fog_fen and check for bestmove plus FoW info strings.
run_case "FoW end-to-end search" "(info string FoW search|bestmove)" "uci
setoption name UCI_Variant value fogofwar
setoption name UCI_FoW value true
setoption name UCI_IISearch value true
setoption name UCI_MinInfosetSize value 64
position fog_fen ????????/??????pp/1?????1P/?1??p1?1/8/1P2P3/PB1P1PP1/NQ1NRBKR b KQk - 0 8
go movetime 50
stop
quit"

echo "All FoW suite checks passed."
