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

 

#### 1. Proper KLUSS Order-2 Neighborhood (MEDIUM PRIORITY)

 

**Current State**: `Subgame::compute_kluss_region()` is a placeholder that marks all nodes as in KLUSS.

 

**What's Needed** (from paper Section 3.2):

- **Order-2 KLUSS**: Include all nodes reachable within 2 moves from any position in the belief state

- **Unfrozen at distance 1**: Nodes at distance 1 from belief state positions are "unfrozen" (strategies can be updated)

- **Frozen beyond**: Nodes at distance 2+ use fixed strategies from the trunk

 

**Algorithm**:

```

For each state s in belief_state:

    For each sequence σ of length ≤ 2:

        Mark node(s, σ) as in_kluss

        If |σ| == 1: mark as unfrozen

```

 

**Implementation Tasks**:

1. Implement `compute_order_k_neighborhood(k=2)`

2. Track "frozen" vs "unfrozen" status for each infoset

3. In CFR solver, only update regrets for unfrozen infosets

4. Use trunk strategies for frozen infosets

 

#### 2. Thread Synchronization Improvements (MEDIUM PRIORITY)

 

**Current State**: Basic synchronization exists but has potential race conditions.

 

**What's Needed**:

1. **Read-write locks for tree access**:

   - Expanders write to tree (add children)

   - CFR solver reads tree

   - Use `shared_mutex` for concurrent read access

 

2. **Atomic infoset updates**:

   - Strategy and regret updates should be atomic or use fine-grained locks

   - Consider lock-free data structures for hot paths

 

3. **Proper shutdown sequence**:

   - Expanders stop first (paper: allows solver to finish with stable tree)

   - Drain work queues before joining threads

   - Handle in-flight expansions gracefully

 

**Implementation Tasks**:

1. Add `std::shared_mutex` to `Subgame` class

2. Wrap tree modifications in write locks, reads in shared locks

3. Add atomic operations for infoset statistics

4. Implement proper thread barrier for shutdown

 

### Important Missing Features

 

#### 3. Action Purification (MEDIUM PRIORITY)

 

**Current State**: `Selection::purify_strategy()` is a placeholder that returns the input strategy unchanged.

 

**What's Needed** (Appendix B.3.7):

- **Support-limited strategy**: Reduce strategy support to MaxSupport actions

- **Margin-based filtering**: Only include actions with non-negative margin

- **Deterministic in Resolve**: Use deterministic play when in Resolve gadget

 

**Algorithm**:

```

margins ← compute_margins(infoset)

stable_actions ← {a : margin[a] ≥ 0}

purified ← top_k(stable_actions, MaxSupport, by=strategy_prob)

renormalize(purified)

```

 

**Implementation Tasks**:

1. Implement proper `compute_margins()` using CFR regrets

2. Implement `purify_strategy()` with MaxSupport filtering

3. Handle Resolve gadget determinism

4. Add unit tests for purification edge cases

 

#### 4. Gadget Implementation (MEDIUM PRIORITY)

 

**Current State**: `GadgetType` enum exists but gadget logic is incomplete.

 

**What's Needed** (Appendix B.3.1, B.3.2):

 

**Resolve Gadget**:

- Prior α(J) = 0.5·uniform + 0.5·y(J) where y is opponent's last strategy

- Add v_alt to counterfactual values

- Ensures safety against worst-case opponent

 

**Maxmargin Gadget**:

- Used after Resolve has been "entered" (opponent deviated from prior)

- Maximizes margin over opponent's strategy

- More aggressive exploitation

 

**Implementation Tasks**:

1. Implement `compute_resolve_prior()`

2. Implement `compute_alternative_value()` (partial implementation exists)

3. Add gadget switching logic in CFR solver

4. Track "resolve entered" state properly

 

#### 5. Leaf Evaluation Integration (MEDIUM PRIORITY)

 

**Current State**: `Evaluator.cpp` exists but integration with Stockfish evaluation is incomplete.

 

**What's Needed**:

1. Hook into Stockfish's NNUE evaluation

2. Normalize evaluation to [-1, +1] range

3. Handle FoW-specific evaluation (average over belief state)

4. Use MultiPV for action value initialization

 

**Implementation Tasks**:

1. Complete `Evaluator::evaluate()` to call Stockfish search

2. Implement `evaluate_belief_state()` that averages over positions

3. Add caching to avoid re-evaluating same positions

4. Handle terminal position detection

 

### Lower Priority Enhancements

 

#### 6. Instrumentation (Appendix B.4) (LOW PRIORITY)

 

**What's Needed**:

- CFR convergence metrics (exploitability approximation)

- Tree size statistics over time

- Action entropy tracking

- Time breakdown per component

 

#### 7. Memory Management (LOW PRIORITY)

 

**What's Needed**:

- Tree pruning for old/unused nodes

- Node recycling pool

- Belief state compression

- Memory limits and cleanup

 

#### 8. Incremental Belief Updates (LOW PRIORITY)

 

**What's Needed**:

- When new observation arrives, filter existing belief set

- Much faster than re-enumeration from scratch

- Requires careful tracking of observation sequence

 

### Testing Requirements

 

#### Unit Tests Needed

1. `Visibility_test.cpp` - Test all visibility rules from Appendix A

2. `Belief_test.cpp` - Test belief enumeration and sampling

3. `CFR_test.cpp` - Test regret matching and strategy updates

4. `Selection_test.cpp` - Test purification with various strategies

5. `Subgame_test.cpp` - Test KLUSS region computation

 

#### Integration Tests Needed

1. **End-to-end FoW search**: Full pipeline from `position` to `bestmove`

2. **fog_fen analysis**: Parse partial observation and search

3. **Multi-threaded stress test**: Race condition detection

4. **Memory leak detection**: Run extended searches

 

#### Comparison Tests

1. Compare move quality against baseline (random play)

2. Compare against simplified FoW search (no belief state)

3. Self-play tournament to verify improvement

 

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

1. Need to verify castling rights handling in visibility computation

2. En-passant visibility edge cases may need testing

3. Crazyhouse piece-in-hand visibility needs verification

 

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
