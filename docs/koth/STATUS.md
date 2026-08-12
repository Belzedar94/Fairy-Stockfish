# Project Status

Last updated: 2026-08-12.

Current state: `GO_NEXT_PHASE` for bounded local engineering. No release or
strength claim exists.

## Frozen inputs

- source base: Fairy-Stockfish
  `c19b5f6c66894fdb0e88d0dd100e3885f744760a`;
- rule profile: `KOTH_LICHESS_V1`;
- exact goal set: `d4`, `e4`, `d5`, `e5`;
- legacy network identity: SHA-256
  `978B86D0E6A45E05F9F1375DCED129CEA0ACEA13041EA65960691632EDC47AF7`;
- public network distribution: held pending an explicit license;
- local engineering ceiling: two CPU threads, no GPU, no shared workers.

## Gate boundary

| Gate | State | Meaning |
| --- | --- | --- |
| Discovery | Passed for engineering | Primary rule and source identities are frozen. |
| Source/build | In progress | Reproducible project builds and deterministic digests are not yet certified. |
| Referee | Open | The unmodified tournament referee has known result and notation conflicts. |
| Network loader | Failing diagnostic | Alias routing and full-byte authentication require engineering. |
| CI/artifacts | Open | No project artifact is currently a release candidate. |
| Search/strength | Not authorized | Correctness and asset gates must pass first. |
| Data/model | Not authorized | No KOTH dataset or V2 model has been selected. |
| Release | Not authorized | A complete draft and owner approval are required. |

Smoke tests, canaries, loss investigations, and short matches remain in their
own evidence classes. They are not Elo, correctness, or release evidence.
