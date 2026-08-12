# Independent executable rules reference

Status: implementation candidate; differential certification is pending.

The project reference is `tools/koth_reference.py`. It uses standard physical
move legality from `chess==1.11.2` and implements the frozen
`KOTH_LICHESS_V1` trajectory, terminal-predicate, precedence, notation, and
loaded-state policies in project-owned code.

It is intentionally independent from the engine KOTH modules and from the
tournament referee. It does not evaluate positions, search moves, load NNUE
weights, manage processes, or adjudicate clocks.

## Pinned dependency

| Field | Identity |
| --- | --- |
| Package | `chess==1.11.2` |
| Upstream source revision | `9c24454dcea4f8a30259d811a2f10b26e911deb4` |
| Source archive | `chess-1.11.2.tar.gz` |
| Source archive SHA-256 | `A8B43E5678FDB3000695BDAA573117AD683761E5CA38E591C4826EBA6D25BB39` |
| License | GPL-3.0+ |

The exact runtime interpreter, script bytes, installed package inventory, and
their hashes belong in each execution receipt. A package version string alone
is not an executable identity.

## Input contract

The CLI accepts one UTF-8 JSON object on standard input:

```json
{
  "schema": "koth-reference-input-v1",
  "rule_profile": "KOTH_LICHESS_V1",
  "root_fen": "7k/8/8/8/8/2K5/8/8 w - - 0 1",
  "moves": ["c3d4"]
}
```

The root must be a canonical six-field standard-chess FEN. A physically
invalid root, an already-terminal root, or either king already occupying a
goal is rejected. Every UCI move token must be legal in sequence, and any
token after the first terminal transition is rejected.

Successful output includes each physical transition, SAN, post-FEN, sorted
physical legal-move digest, all terminal predicates, primary reason, winner,
and public result. Physical legal moves after a hill transition may exist; the
game-legal continuation set is still empty.

## Usage

```sh
python tools/koth_reference.py identity
python tools/koth_reference.py replay < input.json
python -m unittest tests/koth_reference_test.py
```

The reference must disagree loudly. It is not permitted to normalize a WLD,
SAN, terminal-ply, predicate, or loaded-state mismatch into agreement.
