# KOTH-Stockfish runtime contract

KOTH-Stockfish implements the `KOTH_LICHESS_V1` King of the Hill profile on an
official Stockfish development base. The goal squares are `d4`, `e4`, `d5`, and
`e5`. A legal king move onto one of those squares is terminal at the accepted
move boundary. Checkmate takes precedence over a simultaneous hill win; a hill
win takes precedence over stalemate and automatic draws.

## Required legacy network

The engineering runtime does not embed, download, or bundle an evaluation
network. Supply the external owner-authored legacy network with one of these
case-sensitive filenames:

- `kingofthehill-978b86d0e6a4.nnue` (canonical and default)
- `KOTH_v1.nnue` (byte-identical compatibility alias only)

The required file is 47,721,371 bytes with SHA-256
`978B86D0E6A45E05F9F1375DCED129CEA0ACEA13041EA65960691632EDC47AF7`.
The engine authenticates the complete file before activation. A missing,
renamed, truncated, extended, corrupt, or incompatible file aborts rather than
falling back to another evaluator.

Configure an absolute path before `isready` or a nonterminal search:

```text
setoption name EvalFile value C:\absolute\path\kingofthehill-978b86d0e6a4.nnue
isready
```

`KOTH_v1.nnue` does not select a different model and is never the default.

## Protocol and referee

`UCI_Variant` has exactly one supported value: `kingofthehill`. The engine emits
`info string koth ...` records for authoritative terminal states and
`info string koth pv_terminal ...` when a principal variation reaches a KOTH
terminal. Those records preserve the terminal predicate bitmask, primary
reason, and winner; numeric mate scores alone do not distinguish checkmate from
a hill win.

The upstream `bench`, `speedtest`, `go mate`, network export, Chess960, Syzygy,
and WDL interfaces are unsupported for this profile. Use `tools/koth_bench.py`
for the fixed-work semantic digest. The certified independent referee is
`tools/koth_referee.py`; install its exact dependency with:

```text
python -m pip install -r tools/koth-referee-requirements.txt
```

## T0 engineering scope

The T0 search baseline disables transposition-table authority and makes Syzygy
unreachable. This is a correctness and reproducibility baseline, not evidence
of playing strength, model selection, release readiness, or OpenBench routing.
The network remains an external runtime input and is not included in engine
packages.
