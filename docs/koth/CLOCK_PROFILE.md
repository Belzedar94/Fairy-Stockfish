# KOTH_CLOCK_V1

Status: frozen engineering contract; referee certification is pending.

`KOTH_CLOCK_V1` defines the deterministic clock boundary for certified local
KOTH-Stockfish engine matches. It is deliberately narrower than an online
chess server. It has no network-lag compensation, human input allowance, or
implementation-dependent grace period.

## Supported time-control form

The certified lane accepts sudden-death time controls of the form
`base_ms + increment_ms`, with non-negative signed 64-bit integer values.
Moves-per-period, delay, hourglass, byoyomi, and node-based clocks are rejected
by the certified referee profile. Fixed-node work remains a test mode, not a
clocked game.

Each side starts its first engine-controlled turn with exactly `base_ms`.
Opening-book or setup moves are an authenticated root prefix: they consume no
clock and receive no increment. The first timed position, complete opening
prefix, side to move, and both initial clocks are persisted before `go`.

## Timing boundary

The referee uses one process-wide monotonic clock with nanosecond readings.
Wall time, timezone changes, and system-clock adjustment are irrelevant.

For each proposal:

1. `go_sent_ns` is sampled immediately after the complete newline-terminated
   UCI `go` command has been written and flushed to the engine pipe.
2. `proposal_received_ns` is sampled when the referee has received the final
   byte of a complete newline-terminated `bestmove` line.
3. `elapsed_ns = max(0, proposal_received_ns - go_sent_ns)`.
4. `elapsed_ms = elapsed_ns / 1_000_000`, using integer truncation toward zero.
5. `remaining_before_increment_ms = pre_clock_ms - elapsed_ms`, using checked
   signed 64-bit arithmetic.

The proposal is out of time when
`remaining_before_increment_ms <= 0`. Equality is therefore a timeout. The
profile has `lag_allowance_ms = 0` and `expiry_margin_ms = 0`.

If the proposal is timely and legal, the move is committed once and then
`increment_ms` is added with checked arithmetic. Increment is never granted to
a timed-out, malformed, illegal, missing, or rejected proposal.

## Transaction order

The live referee must use this order without make/undo/make persistence:

1. capture the pre-position, legal-move digest, clocks, and identities;
2. send and flush `go`, then sample `go_sent_ns`;
3. receive a complete proposal and sample `proposal_received_ns`;
4. compute the clock decision;
5. if timed out, persist a rejected-proposal receipt and end the game without
   changing board state, UCI history, SAN, PGN moves, or repetition history;
6. otherwise validate the proposal against the pre-position legal moves;
7. compute SAN from the pre-position;
8. commit the accepted move exactly once;
9. add increment, compute terminal status, and persist the accepted record;
10. either terminate or send the updated position to the engines.

A timely hill move can win only after step 8. A late hill proposal loses on
time at step 5 and never becomes a physical or recorded board transition.

For a valid KOTH trajectory, the opponent retains winning potential while its
king exists. Orthodox insufficient-material logic must not convert a KOTH
timeout into a draw.

## Failure and overflow policy

- A monotonic-clock regression is `KCLOCK_E_MONOTONIC` and aborts the game.
- Arithmetic overflow is `KCLOCK_E_RANGE` and aborts the game.
- Unsupported time-control syntax is `KCLOCK_E_UNSUPPORTED` and rejects game
  creation.
- A partial or unterminated `bestmove` line is not a received proposal.
- Engine crash, protocol EOF, malformed move, and illegal move are distinct
  forfeits and are never rewritten as timeouts.
- Timeout, crash, and protocol-forfeit records preserve the unchanged
  pre-position and do not contain a committed `post_fen`.

Canonical diagnostics contain stable codes and numeric fields, never localized
OS text. Raw process transcripts are retained separately.

## Required receipt fields

Every proposal receipt must include at least:

- `clock_profile = KOTH_CLOCK_V1`;
- `pre_clock_ms`, `go_sent_ns`, `proposal_received_ns`, `elapsed_ns`, and
  `elapsed_ms`;
- `lag_allowance_ms = 0`, `expiry_margin_ms = 0`, and `increment_ms`;
- the clock decision and stable reason code;
- whether the proposal was accepted;
- pre-position identity and, only for an accepted move, post-position identity;
- engine, referee, rule-profile, and time-control identities.

## Certification boundary

This contract is not a clock certificate. Certification requires injected
boundary tests below, at, and above expiry; real process-level tests; timeout
state immutability; increment tests; book-prefix tests; crash/protocol tests;
and independent recomputation from the recorded monotonic timestamps.
