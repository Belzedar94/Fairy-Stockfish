# Frozen Fairy-Stockfish Atomic baseline on OpenBench

This branch is a non-playing harness pinned to upstream Fairy-Stockfish commit
`fb78cb561aa01708338e35b3dc3b65a42149a3c4`. It exists only to make the
historical Atomic baseline reproducible on public OpenBench workers; the oracle
checkout remains untouched.

The worker contract from `src/` is:

```text
make -j EXE=<output> CXX=<compiler> EVALFILE=<absolute-network-path>
```

When `EVALFILE` is present, the default goal builds `x86-64-bmi2` with
`all=no`, `largeboards=no`, and `nnue=yes`. It uses MinGW on Windows and GCC on
Linux, while preserving the worker's `EXE` and `CXX`. A normal developer
`make` without `EVALFILE` keeps Fairy-Stockfish's original `help` default.
Fairy's historical native-MinGW path omits LTO, so the shim supplies the exact
`-flto -flto-partition=one` compile/link flags used by the frozen optimized
baseline, plus `-save-temps -Wl,--no-insert-timestamp` for deterministic PE
linking. These additions are scoped to the Windows OpenBench target; ordinary
builds and Linux builds retain their upstream flags.

OpenBench authenticates the assigned network before invoking Make. The shim
copies that file to Fairy-Stockfish's canonical embedded name
`atomic_run3b_e202_l05.nnue` and calls the private optimized `all` target
directly. It intentionally does not call `build nnue=yes`, whose `net`
prerequisite would reject the Atomic filename. The harness macro changes
Fairy's compiled EvalFile default to that Atomic-prefixed canonical name; this
is required by Fairy's variant-to-network safety check and leaves non-harness
builds unchanged.

The frozen Atomic network is `atomic_run3b_e202_l05.nnue`, SHA-256
`99DC67EABF26A64FAEECA3A88B4C38597A840B8D4A874B9F2CF658C6F92A04A6`.
The public binary contains that network and does not depend on the worker's
temporary download path.

`OPENBENCH_ATOMIC_BASELINE` changes only missing defaults in the benchmark
command. Bare `bench` is therefore equivalent to:

```text
bench atomic 16 1 13 default depth NNUE
```

Its deterministic signature with the frozen network is:

```text
Nodes searched  : 97362
```

Explicit benchmark arguments still take precedence. Search, evaluation,
variant rules, protocols, and all game-time options are unchanged.
