# KOTH-Stockfish

KOTH-Stockfish is a correctness-first King of the Hill specialization of
[Fairy-Stockfish](https://github.com/fairy-stockfish/Fairy-Stockfish). The project
starts from upstream commit
`c19b5f6c66894fdb0e88d0dd100e3885f744760a` and targets the explicitly versioned
`KOTH_LICHESS_V1` rule profile.

This repository is under active engineering. There is no stable KOTH-Stockfish
release yet, and the current default branch must not be treated as a strength
claim or a certified tournament binary.

## Rule contract

The project does not treat “distance to the center” as a game rule. A King of
the Hill win is a physical, legal transition in which the previous mover's king
occupies one of exactly four goal squares: `d4`, `e4`, `d5`, or `e5`.

The complete versioned contract, including terminal ordering, draw handling,
loaded-position policy, notation, and persistence requirements, is in
[`docs/koth/RULE_PROFILE.md`](docs/koth/RULE_PROFILE.md).

## Legacy network

The official legacy compatibility network is an external asset identified by
its complete SHA-256 digest. Its bytes are not stored in this repository and
must not be redistributed from this repository. `KOTH_v1.nnue` is reserved as a
byte-identical alias; it is not a new network and will not change the engine's
default until its loader contract is certified.

See [`docs/koth/LEGACY_NETWORK.md`](docs/koth/LEGACY_NETWORK.md) and
[`networks/manifest.json`](networks/manifest.json) for the exact identity and
current compatibility status.

## Current status

The public engineering state and gate boundary are recorded in
[`docs/koth/STATUS.md`](docs/koth/STATUS.md). Results from smoke tests, canaries,
or short matches are not promoted to Elo, correctness, or release evidence.

## Build

Until the KOTH build profile is frozen, use the upstream Fairy-Stockfish build
interface from `src`:

```sh
make -j2 build ARCH=x86-64
```

Build outputs are not release artifacts unless they are accompanied by the
project's reproducibility, correctness, referee, network-loader, and runtime
receipts.

## License and provenance

Engine source code is distributed under the GNU General Public License v3.0.
See [`Copying.txt`](Copying.txt). This fork retains the Fairy-Stockfish history
and attribution in [`AUTHORS`](AUTHORS).

Neural-network weights are separate assets. The GPL license for the engine does
not by itself grant redistribution rights for an external network.
