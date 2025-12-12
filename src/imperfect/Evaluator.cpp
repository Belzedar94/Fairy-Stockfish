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

/// simple_material_eval() computes a basic material-only evaluation
/// Used when pos.this_thread() is null and Eval::evaluate() cannot be called
Value simple_material_eval(const Position& pos) {
    // Simple material count: pawns=100, knights/bishops=300, rooks=500, queens=900
    constexpr int PawnValue = 100;
    constexpr int KnightValue = 300;
    constexpr int BishopValue = 300;
    constexpr int RookValue = 500;
    constexpr int QueenValue = 900;

    int score = 0;
    Color us = pos.side_to_move();
    Color them = ~us;

    score += (pos.count<PAWN>(us) - pos.count<PAWN>(them)) * PawnValue;
    score += (pos.count<KNIGHT>(us) - pos.count<KNIGHT>(them)) * KnightValue;
    score += (pos.count<BISHOP>(us) - pos.count<BISHOP>(them)) * BishopValue;
    score += (pos.count<ROOK>(us) - pos.count<ROOK>(them)) * RookValue;
    score += (pos.count<QUEEN>(us) - pos.count<QUEEN>(them)) * QueenValue;

    return Value(score);
}

/// evaluate() returns a static evaluation of the position, normalized to [-1, +1]
/// Uses Eval::evaluate() for static evaluation instead of running a search
/// to avoid threading issues when called from FoW planner threads.
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

    // Use static evaluation
    // Note: Eval::evaluate() requires pos.this_thread() to be non-null
    // because it accesses thread-specific tables. When called from FoW
    // threads where Position is created with nullptr thread, we fall back
    // to a simple material-only evaluation.
    Value evalValue;
    if (pos.this_thread())
        evalValue = Eval::evaluate(pos);
    else
        evalValue = simple_material_eval(pos);

    float normalized = normalize_value(evalValue);

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
