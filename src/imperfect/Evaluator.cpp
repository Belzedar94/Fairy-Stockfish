/*
  Fairy-Stockfish, a UCI chess variant playing engine derived from Stockfish
  Copyright (C) 2018-2024 Fabian Fichter

  Fairy-Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Fairy-Stockfish is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Evaluator.h"
#include "../evaluate.h"
#include "../movegen.h"
#include "../position.h"
#include "../search.h"
#include "../thread.h"
#include <algorithm>
#include <cmath>
#include <mutex>
#include <shared_mutex>
#include <unordered_map>

namespace Stockfish {
namespace FogOfWar {

namespace {

std::mutex searchMutex;
std::shared_mutex cacheMutex;
std::unordered_map<Key, float> evaluationCache;

} // namespace

/// normalize_value() converts centipawn evaluation to [-1, +1]
/// Mate scores map to +/-1, material scores are clamped
float normalize_value(Value v) {
    // Handle mate scores
    if (v >= VALUE_MATE_IN_MAX_PLY)
        return 1.0f;
    if (v <= VALUE_MATED_IN_MAX_PLY)
        return -1.0f;

    // Normalize material scores using sigmoid-like curve
    // Map roughly [-1000, +1000] centipawns to [-1, +1]
    constexpr float scale = 1000.0f;
    float normalized = float(v) / scale;

    // Clamp to [-1, +1]
    return std::max(-1.0f, std::min(1.0f, normalized));
}

/// evaluate() runs a shallow search (depth=1) and normalizes the result
/// Also detects terminal positions before launching a search
float evaluate(Position& pos) {
    // Terminal detection
    Value terminalValue;
    if (pos.is_game_end(terminalValue, 0))
        return normalize_value(terminalValue);

    // Cache lookup
    Key key = pos.key();
    {
        std::shared_lock<std::shared_mutex> lock(cacheMutex);
        auto it = evaluationCache.find(key);
        if (it != evaluationCache.end())
            return it->second;
    }

    std::lock_guard<std::mutex> guard(searchMutex);

    // Re-check cache after acquiring the search lock to avoid duplicate work
    {
        std::shared_lock<std::shared_mutex> lock(cacheMutex);
        auto it = evaluationCache.find(key);
        if (it != evaluationCache.end())
            return it->second;
    }

    Search::LimitsType limits;
    limits.depth = 1; // Depth-1 as described in the paper

    StateListPtr states(new std::deque<StateInfo>(1));
    if (pos.state())
        (*states)[0] = *pos.state();

    Threads.start_thinking(pos, states, limits, false);
    Threads.main()->wait_for_search_finished();

    Value searchValue = VALUE_ZERO;
    const auto& rootMoves = Threads.main()->rootMoves;
    if (!rootMoves.empty())
        searchValue = rootMoves.front().score;

    float normalized = normalize_value(searchValue);

    {
        std::unique_lock<std::shared_mutex> lock(cacheMutex);
        evaluationCache[key] = normalized;
    }

    return normalized;
}

/// evaluate_children() evaluates all legal child positions
/// Implements the depth-1 MultiPV evaluation described in Appendix B.3.4
std::vector<ChildEvaluation> evaluate_children(Position& pos) {
    std::vector<ChildEvaluation> evaluations;

    // Generate all legal moves
    for (const auto& m : MoveList<LEGAL>(pos)) {
        StateInfo st;
        // Make the move
        pos.do_move(m, st);

        // Evaluate the resulting position using shallow search
        float eval = -evaluate(pos);

        // Undo the move
        pos.undo_move(m);

        // Store normalized evaluation
        ChildEvaluation ce{m, eval};
        evaluations.push_back(ce);
    }

    return evaluations;
}

float evaluate_belief_state(const BeliefState& beliefState, const Variant* variant) {
    if (!variant || beliefState.empty())
        return 0.0f;

    const auto& states = beliefState.all_states();
    if (states.empty())
        return 0.0f;

    float total = 0.0f;
    size_t count = 0;

    for (const auto& fen : states) {
        StateInfo st;
        Position pos;
        pos.set(variant, fen, variant->chess960, &st, nullptr, true);

        total += evaluate(pos);
        ++count;
    }

    return count ? total / static_cast<float>(count) : 0.0f;
}

/// get_best_child() returns the move with the highest evaluation
Move get_best_child(const std::vector<ChildEvaluation>& evals) {
    if (evals.empty())
        return MOVE_NONE;

    auto best = std::max_element(evals.begin(), evals.end(),
        [](const ChildEvaluation& a, const ChildEvaluation& b) {
            return a.value < b.value;
        });

    return best->move;
}

} // namespace FogOfWar
} // namespace Stockfish
