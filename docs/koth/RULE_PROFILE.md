# KOTH_LICHESS_V1 Rule Profile

Status: owner-approved engineering contract. Referee certification is still
pending.

## Identity

- Project profile: `KOTH_LICHESS_V1`
- Public variant name: `King of the Hill`
- UCI variant value: `kingofthehill`
- Goal squares: `d4`, `e4`, `d5`, `e5`
- Board and orthodox move-legality basis: standard chess

`KOTH_LICHESS_V1` is a KOTH-Stockfish version identifier. It is not a name
assigned by Lichess.

## Legal hill entry

A king may enter or remain on a goal square only through a move that is legal
under standard king-safety rules. The destination must not be attacked after
the move. The goal set does not override check, pins, occupancy, castling,
en-passant, or any other legality constraint.

In a normal trajectory, a hill win becomes terminal immediately after an
accepted legal move leaves the previous mover's king on a goal square. A king
approaching a goal, or a non-king piece occupying one, has no terminal meaning.

## Terminal ordering

When more than one terminal predicate is true after the same accepted move, the
primary status order is:

1. checkmate;
2. hill/variant end;
3. stalemate;
4. automatic draw.

The lossless game record must retain every true predicate as well as the primary
reason. In particular, a legal hill transition takes precedence over the
100-halfmove and repetition draw predicates when they arise together, unless
the position is checkmate under the ordering above.

For deterministic headless engine games:

- threefold repetition is an automatic claim;
- fivefold repetition is an automatic draw;
- a halfmove clock of at least 100 is an automatic draw;
- orthodox insufficient-material adjudication is disabled for this variant.

Clock adjudication precedes board progress: a move classified as out of time
after clock and lag handling is a time loss, and its board transition is not
committed as the game result.

Certified local matches use the exact `KOTH_CLOCK_V1` profile in
[`CLOCK_PROFILE.md`](CLOCK_PROFILE.md). Online-server lag compensation is
evidence about Lichess operation, not an implicit input to local adjudication.

## Loaded positions

FEN does not encode which player made the preceding move or authenticate a
terminal transition. A loaded position with a king already on a goal square is
therefore rejected as ambiguous unless authenticated predecessor data proves
the accepted move, mover, and resulting terminal predicate set.

This fail-closed rule applies to engine tests, referee inputs, data generation,
and match runners. A position must not infer a winner merely from a king's
presence on a goal square without the required trajectory evidence.

## Protocol, notation, and persistence

- UCI selects the variant with `UCI_Variant=kingofthehill`.
- FEN/EPD describes physical position state but is not a lossless terminal
  record for an already-completed game.
- SAN for a pure hill-winning move ends in `#`.
- The public machine result uses `variantEnd` with the correct winner.
- Lichess-compatible PGN uses the normal game result and `Termination "Normal"`
  for a pure hill win. That PGN field alone is not a lossless terminal reason.
- A project game record must persist the raw PGN, physical board, side to move,
  castling rights, en-passant/last-move information, halfmove clock, repetition
  history, clocks, all terminal predicates, primary reason, rule-profile
  version, and referee identity.

## Evaluation and data boundary

Goal membership is derived from the physical board. “Distance to the center”
and other search features are hypotheses, not rules or physical labels.

Training labels use the terminal result from the recorded player's perspective.
No board symmetry is authorized by geometry alone: pawn direction, castling
rights, en-passant state, side to move, and trajectory history must all survive
an exact round trip before a transformation can be enabled.

## Certification boundary

This document freezes the intended dialect but is not a referee certificate.
Certification requires differential fixtures for every goal, legal and illegal
entries, terminal precedence, clocks, notation, loaded states, and replay
persistence across the engine, executable reference, and tournament referee.
