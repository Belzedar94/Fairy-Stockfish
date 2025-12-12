# Obscuro-Style Fog-of-War Chess Implementation for Fairy-Stockfish

## Implementation Status

This implementation follows the Obscuro paper's algorithms for imperfect-information search in Fog-of-War chess.

### Completed Modules

#### Phase 1: Foundations
- [x] **Visibility Module** (`src/imperfect/Visibility.{h,cpp}`)
  - Computes visible squares under FoW rules (Appendix A)
  - Special handling for blocked pawns, en-passant visibility
  - Implements all paper's visibility rules

- [x] **Belief State** (`src/imperfect/Belief.{h,cpp}`)
  - Maintains set P of consistent positions
  - Observation history tracking
  - From-scratch enumeration (Figure 9, lines 2-4)
  - Enumerates hidden opponent piece permutations consistent with `fog_fen` observations (capped at 1024 states)

- [x] **Evaluator Hook** (`src/imperfect/Evaluator.{h,cpp}`)
  - MultiPV depth-1 evaluation for all children
  - Normalized to [-1, +1] range (Appendix B.3.4)
  - Used as leaf heuristic for expansion

#### Phase 2: Core CFR Engine
- [x] **Subgame & KLUSS** (`src/imperfect/Subgame.{h,cpp}`)
  - 2-KLUSS implementation (order-2 knowledge region, unfrozen at distance 1)
  - Infoset representation with sequence IDs
  - Resolve and Maxmargin gadget games (Appendix B.3.1, B.3.2)
  - Gift and alternative value calculations

- [x] **PCFR+ Solver** (`src/imperfect/CFR.{h,cpp}`)
  - Predictive CFR+ with PRM+ (Positive Regret Matching Plus)
  - Last-iterate play (no strategy averaging at runtime)
  - Gadget switching (Resolve ↔ Maxmargin) (Figure 10, lines 6-13)
  - Alternative value handling for Resolve (Figure 10, lines 18-19)

- [x] **GT-CFR Expander** (`src/imperfect/Expander.{h,cpp}`)
  - One-sided GT-CFR expansion (Appendix B.3.3)
  - PUCT selection with variance (C=1.0)
  - Alternating exploring side
  - Initialize to best child (Appendix B.3.4, Figure 12 lines 24-27)

- [x] **Action Selection & Purification** (`src/imperfect/Selection.{h,cpp}`)
  - Purified strategy with MaxSupport=3 (Appendix B.3.7)
  - Deterministic play in Resolve gadget
  - Stable action filtering (non-negative margins)

#### Phase 3: Coordination & Integration
- [x] **Planner** (`src/imperfect/Planner.{h,cpp}`)
  - Main coordinator implementing Figure 8 (Move loop)
  - Threading: 1 CFR solver + 2 expanders
  - Time management
  - Statistics collection

- [x] **UCI Integration**
  - UCI options for FoW configuration
  - Hook into `go` command in `uci.cpp`
  - Option defaults match paper parameters

### UCI Options Added

```
UCI_FoW               = false      // Enable Fog-of-War mode
UCI_IISearch          = true       // Enable imperfect information search
UCI_MinInfosetSize    = 256        // Sample size for root infoset (Figure 9)
UCI_ExpansionThreads  = 2          // Number of expander threads
UCI_CFRThreads        = 1          // Number of CFR solver threads
UCI_PurifySupport     = 3          // Max actions in purified strategy
UCI_PUCT_C            = 100        // PUCT constant C (x100 for precision)
UCI_FoW_TimeMs        = 5000       // Time budget per move in milliseconds
```

### Paper Algorithm Mapping

| Paper Reference | Implementation File | Status |
|----------------|---------------------|--------|
| Figure 8 (Move loop) | `Planner.cpp::plan_move()` | ✓ |
| Figure 9 (ConstructSubgame) | `Planner.cpp::construct_subgame()` | ✓ |
| Figure 10 (RunSolverThread) | `CFR.cpp::run_continuous()` | ✓ |
| Figure 11 (MakeUtilities) | `CFR.cpp::compute_cfv()` | ✓ |
| Figure 12 (RunExpanderThread) | `Expander.cpp::run_continuous()` | ✓ |
| Appendix A (FoW Rules) | `Visibility.cpp` | ✓ |
| Appendix B.3.1 (Gadgets) | `Subgame.cpp::compute_alternative_value()` | ✓ |
| Appendix B.3.2 (Resolve prior) | `Subgame.cpp::construct()` | ✓ |
| Appendix B.3.3 (One-sided GT-CFR) | `Expander.cpp::select_leaf()` | ✓ |
| Appendix B.3.4 (Leaf init) | `Expander.cpp::initialize_to_best_child()` | ✓ |
| Appendix B.3.6 (PCFR+) | `CFR.cpp::run_iteration()` | ✓ |
| Appendix B.3.7 (Purification) | `Selection.cpp::purify_strategy()` | ✓ |

### Recent Bug Fixes (December 2025)

 

The following critical bugs were fixed to make the FoW search functional:

 

1. **Selection.cpp crash**: `compute_margins()` called `std::max_element` on empty `qValues` vector

   - Fixed by adding empty checks and returning zero margins when qValues is empty

 

2. **Position.cpp crash**: `do_move()` accessed `thisThread->nodes` when `thisThread` was null (FoW context)

   - Fixed by adding null check before incrementing node counter

 

3. **CFR.cpp race condition**: `compute_cfv()` accessed `node->children` while Expander was modifying it

   - Fixed by checking `expanded` flag and taking snapshot of children size

 

4. **Expander.cpp wrong variant**: Used hardcoded "chess" variant instead of actual variant

   - Fixed by storing variant pointer in Subgame and retrieving it

 

5. **Uninitialized root infoset**: Root infoset had empty actions vector

   - Fixed by populating with legal moves from position in `construct_subgame()`

 

6. **fog_fen support**: Added basic parsing and storage of fog_fen (partial observation FEN)

   - New commands: `position fog_fen <partial_fen>` stores the observation

   - Functions `get_fog_fen()` and `clear_fog_fen()` added to uci.h

 

---

 

## What's Missing for Complete Implementation

 

This section provides a detailed analysis of what remains to be implemented to achieve a production-quality Obscuro-style FoW search engine. Items are organized by priority and complexity.

 

### Critical Missing Features

 

#### Addressed Critical Gaps

The latest update implements the two previously missing foundation pieces:

- **Full belief state enumeration**: `BeliefState::enumerate_candidates()` now permutes hidden opponent pieces across unseen squares (masked by visibility) and keeps every state consistent with the latest observation. Enumeration is capped at 1024 states to prevent combinatorial blowups, and illegal or king-capturable positions are filtered out before sampling.

- **fog_fen integration**: `BeliefState::parse_fog_fen()` converts partial FoW FEN strings (supports `*` or `?` for unknown squares) into observations, and `Planner::construct_subgame()` seeds the belief state directly from a supplied `fog_fen` before running the solver.

 

#### 1. Proper KLUSS Order-2 Neighborhood (DONE)

`Subgame::compute_kluss_region()` now walks the tree and marks nodes within depth 2 as inside the KLUSS region, with infosets inside distance 1 flagged as "unfrozen". Frozen infosets capture and reuse their trunk strategies so CFR only updates regrets for the unfrozen frontier. Sequence IDs are extended per-move, keeping trunk/frozen boundaries aligned with move depth.

 

#### 2. Thread Synchronization Improvements (DONE)

Tree access now relies on a `shared_mutex` to separate read and write paths: expanders take shared locks while selecting leaves and upgrade to unique locks during expansion, and counters/strategies are guarded by infoset-level mutexes. The planner already stops expanders before the solver to drain in-flight work, matching the paper's shutdown order.

 

### Important Missing Features

 

#### 3. Action Purification (DONE)

`Selection::compute_margins()` now derives margins from CFR regrets (falling back to Q-values) and `purify_strategy()` filters to the top MaxSupport non-negative actions. Resolve gadget play is deterministic, and selection respects the purified distribution when collapsing to a single action.

 

#### 4. Gadget Implementation (DONE)

Resolve gadgets now build the paper's prior α(J) via `compute_resolve_prior()` (50/50 uniform/opponent mix) and fold the alternative value back into counterfactuals. Gift computation uses the alternative estimate so switching to Maxmargin preserves the safety guarantee once `resolveEntered` is raised by the solver.

 

#### 5. Leaf Evaluation Integration (DONE)

Depth-1 child evaluation now uses a shallow Stockfish search (depth = 1) with terminal detection and normalization to [-1, +1]. `Evaluator::evaluate()` caches results by Zobrist key to avoid re-running searches on duplicate states, and `evaluate_belief_state()` walks the belief set to return the averaged score. Terminal states (mate/draw/stalemate and variant-specific endings) are intercepted before search to keep heuristics stable.

 

### Lower Priority Enhancements

 

#### 6. Instrumentation (Appendix B.4) (DONE)

- Planner statistics now record exploitability approximation using the mean positive regret across infosets, action entropy at the root, peak node counts, and a timeline of node sizes.
- Time is broken down across construction, search, and selection so the UCI `info string` can report where FoW time was spent.


#### 7. Memory Management (DONE)

- Subgame nodes are pulled from and returned to a recycling pool, keeping allocations bounded.
- Out-of-KLUSS branches are pruned eagerly and a soft node limit triggers pruning of deep or frozen leaves.
- Belief states are compressed to a configurable cap before sampling, shrinking both memory footprint and sampling cost.

 

#### 8. Incremental Belief Updates (DONE)

- `BeliefState::update_incrementally()` now filters the current belief set against the
  latest observation, reuses cached positions when visibility shrinks, and triggers
  a full rebuild only when the observation expands or filtering collapses the set.
- The planner wires incremental updates through `enableIncrementalBelief`, falling
  back to full reconstruction when disabled or when the belief set underflows.

 

### Testing Status

#### Implemented FoW Test Coverage

- `tests/fow_suite.sh` – lightweight suite that now covers:
  - Castling-rights visibility in FoW positions
  - En-passant target visibility
  - Crazyhouse pockets (piece-in-hand) rendering for Dark Crazyhouse
  - End-to-end FoW search seeded via `fog_fen` with IISearch enabled
- `tests/fow_api_checks.sh` – targeted API/visibility validator that confirms castling
  rights, en-passant markers, and darkcrazyhouse pockets remain exposed when positions
  are supplied via `fog_fen`.
- `tests/fow_incremental.sh` – smoke/integration run for fog_fen parsing, incremental
  belief filtering, and the FoW planner pipeline end-to-end.
- `tests/fow_stress.sh` – short multi-cycle regression loop that alternates FoW and
  non-FoW searches to exercise teardown/initialization paths and catch stability
  regressions.
- `tests/fow_compare.sh` – baseline vs FoW comparison harness to validate that standard
  chess searches and FoW searches both produce bestmoves and share a stable setup path.
- `tests/fow_leakcheck.sh` – optional valgrind-driven sanity pass that runs short FoW and
  baseline searches to flag obvious leaks in initialization/teardown.
- `tests/fow_logic.sh` – focused integration checks for belief enumeration on heavy
  fog_fen inputs, KLUSS frontier stability at small infoset sizes, deterministic
  purification, and repeated CFR invocations in a single session.
- `tests/fow_units.sh` – C++ unit suites linked into the engine via `--fow-unittests`
  (visibility, belief enumeration/sampling caps, CFR lifecycle, purification/selection,
  and KLUSS sequence hashing).
- `tests/fow_sanitizers.sh` – ASan/TSan-backed initialization and teardown loops that
  reuse the unit suites under sanitizer builds.
- `tests/fow_benchmarks.sh` – baseline vs FoW micro-benchmark harness that exercises the
  evaluator parity checks alongside a short chess bench.

#### Remaining Gaps

- The new C++ suites are intentionally lightweight; expand with deeper fixture coverage
  (edge-case fog_fen parsing, larger KLUSS frontiers, deeper CFR run loops) as time
  permits.
- Stress loops remain short to keep CI reasonable; add overnight-duration runs locally
  when changing allocation, threading, or belief management paths.

 

### Implementation Roadmap

 

**Phase 1: Core Functionality** (Essential for correct play)

1. ✅ Full belief state enumeration (hidden-piece permutations, capped at 1024)

2. ✅ fog_fen integration (parse + seed belief state)

3. Action purification

4. Thread synchronization fixes

 

**Phase 2: Quality Improvements** (Better play quality)

1. KLUSS order-2 neighborhood

2. Gadget implementation

3. Leaf evaluation integration

4. Instrumentation

 

**Phase 3: Polish** (Production readiness)

1. Memory management

2. Incremental belief updates

3. Comprehensive testing

4. Performance optimization

 

### Known Issues & TODOs

 

#### Remaining API Issues

- Primary API visibility checks have been automated via `tests/fow_api_checks.sh`.
  Continue to watch for edge cases in complex `fog_fen` setups and add targeted C++
  unit cases as the harness matures.

 

#### Performance Concerns

1. Belief enumeration will be slow for complex positions

2. CFR convergence rate unknown for FoW chess

3. Memory usage with large belief states

### Build Instructions

```bash
cd src
make build ARCH=x86-64-modern
```

### Usage

To enable Fog-of-War Obscuro search:

```
setoption name UCI_FoW value true
setoption name UCI_IISearch value true
setoption name UCI_FoW_TimeMs value 5000
go
```

### Architecture

```
src/imperfect/
├── Visibility.{h,cpp}    # FoW visibility computation
├── Belief.{h,cpp}        # Belief state (set P of positions)
├── Evaluator.{h,cpp}     # Leaf evaluation hook
├── Subgame.{h,cpp}       # Game tree, infosets, KLUSS
├── CFR.{h,cpp}           # PCFR+ solver
├── Expander.{h,cpp}      # One-sided GT-CFR expansion
├── Selection.{h,cpp}     # Purification
└── Planner.{h,cpp}       # Main coordinator
```

### Paper Parameters Implemented

| Parameter | Value | Paper Reference |
|-----------|-------|-----------------|
| MinInfosetSize | 256 | Figure 9, line 10 |
| PUCT Constant C | 1.0 | Appendix B.3.3 |
| MaxSupport | 3 | Appendix B.3.7 |
| Solver Threads | 1 | Section 3.4 |
| Expander Threads | 2 | Section 3.4 |
| Variance Prior | {-1, +1} | Appendix B.3.3 |
| Resolve α(J) | 0.5·uniform + 0.5·y(J) | Appendix B.3.2 |

### References

Implementation follows:
- "Combining Private and Public Information for Imperfect-Information Games" (Obscuro paper)
- Fairy-Stockfish engine architecture
- UCT/PUCT Monte Carlo tree search
- Counterfactual Regret Minimization (CFR+)

### License

GPLv3 (consistent with Fairy-Stockfish)
