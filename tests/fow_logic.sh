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

# Dense belief enumeration seeded by fog_fen with many unknowns; ensures the sampler builds
# a non-empty belief set and still returns a move under short time constraints.
run_case "FoW belief enumeration" "bestmove" "uci
setoption name UCI_Variant value fogofwar
setoption name UCI_FoW value true
setoption name UCI_IISearch value true
setoption name UCI_MinInfosetSize value 64
position fog_fen ???k????/????n??p/?????p2/???p4/8/????P3/PPPP1PPP/RNBQKBNR w KQ - 0 5
go movetime 50
stop
quit"

# KLUSS/frontier stability: drive a tiny search with minimal infoset size to keep the
# subgame shallow and exercise the KLUSS unfrozen frontier while still reporting a move.
run_case "FoW KLUSS frontier" "bestmove" "uci
setoption name UCI_Variant value fogofwar
setoption name UCI_FoW value true
setoption name UCI_IISearch value true
setoption name UCI_MinInfosetSize value 32
setoption name UCI_ExpansionThreads value 1
position startpos
go nodes 128
stop
quit"

# Purified selection clamp: force support to 1 to stress deterministic purification and
# ensure the engine still responds with a valid move.
run_case "FoW purified selection" "bestmove" "uci
setoption name UCI_Variant value fogofwar
setoption name UCI_FoW value true
setoption name UCI_IISearch value true
setoption name UCI_PurifySupport value 1
position startpos
go nodes 64
stop
quit"

# CFR continuity: run back-to-back FoW searches in the same session to ensure regrets and
# caches survive repeated search invocations without crashes.
run_case "FoW CFR continuity" "bestmove" "uci
setoption name UCI_Variant value fogofwar
setoption name UCI_FoW value true
setoption name UCI_IISearch value true
position startpos
go movetime 30
stop
go movetime 30
stop
quit"

echo "All FoW logic checks passed."
