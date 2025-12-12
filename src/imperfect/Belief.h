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

#ifndef BELIEF_H_INCLUDED
#define BELIEF_H_INCLUDED

#include <vector>
#include <unordered_set>
#include "../types.h"
#include "../position.h"
#include "Visibility.h"

namespace Stockfish {
namespace FogOfWar {

/// StateKey uniquely identifies a position using Zobrist hashing
using StateKey = uint64_t;

/// Observation represents what the player can see at a given point
struct Observation {
    Bitboard visible = 0;            // Squares that are visible
    Bitboard myPieces = 0;           // Our pieces (always known exactly)
    Bitboard seenOpponentPieces = 0; // Opponent pieces that we can see
    Color sideToMove = WHITE;        // Who moves next
    Bitboard epSquares = 0;          // En-passant squares if visible
    int castlingRights = 0;          // Our known castling rights

    // For reconstruction
    int halfmoveClock = 0;           // 50-move counter
    int fullmoveNumber = 1;          // Full move number
};

/// ObservationHistory maintains the sequence of observations
class ObservationHistory {
public:
    void add_observation(const Observation& obs);
    void clear();
    const std::vector<Observation>& observations() const { return history; }
    size_t size() const { return history.size(); }
    const Observation& last() const { return history.back(); }
    bool empty() const { return history.empty(); }

private:
    std::vector<Observation> history;
};

/// BeliefState maintains the set P of all positions consistent with observations
/// Implements from-scratch enumeration (Figure 9, lines 2-4 of the paper)
class BeliefState {
public:
    BeliefState() = default;

    /// rebuild_from_observations() reconstructs P from scratch given observation history
    /// This is the baseline implementation (Figure 9 of paper)
    void rebuild_from_observations(const ObservationHistory& obsHist,
                                    const Position& truePos);

    /// update_incrementally() attempts incremental update; falls back to rebuild if needed
    void update_incrementally(const Observation& newObs);

    /// Sample a subset of states for building the subgame
    /// Returns FEN strings of sampled positions
    std::vector<std::string> sample_states(size_t n, uint64_t seed = 0) const;

    /// Accessors
    size_t size() const { return stateFens.size(); }
    const std::vector<std::string>& all_states() const { return stateFens; }
    bool empty() const { return stateFens.empty(); }

    /// Check if a position is consistent with observations
    static bool is_consistent(const Position& pos, const Observation& obs);

    /// Parse a partial FoW FEN (fog_fen) into an Observation structure
    static Observation parse_fog_fen(const std::string& fogFen, const Variant* variant);

private:
    std::vector<std::string> stateFens; // FEN strings instead of Position objects
    std::unordered_set<StateKey> stateKeys; // For deduplication

    const Variant* variant = nullptr;
    bool isChess960 = false;
    Thread* owningThread = nullptr;

    /// Generate candidate positions from observations
    void enumerate_candidates(const ObservationHistory& obsHist,
                              const Position& truePos);

    /// Filter candidates by legality and consistency
    void filter_illegal_states();

    /// Helper: Check if king is capturable (game would have ended)
    bool is_king_capturable(const Position& pos) const;

    /// Helper: create a Position from a FEN using the current variant context
    bool set_position_from_fen(Position& pos, StateInfo& st, const std::string& fen) const;
};

/// create_observation() creates an observation from current position
Observation create_observation(const Position& pos);

} // namespace FogOfWar
} // namespace Stockfish

#endif // #ifndef BELIEF_H_INCLUDED
