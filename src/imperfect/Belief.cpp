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

#include <algorithm>
#include <functional>
#include <numeric>
#include <random>
#include <sstream>
#include <unordered_map>
#include "Belief.h"
#include "../movegen.h"
#include "../bitboard.h"

namespace Stockfish {
namespace FogOfWar {

/// ObservationHistory implementation

void ObservationHistory::add_observation(const Observation& obs) {
    history.push_back(obs);
}

void ObservationHistory::clear() {
    history.clear();
}

/// create_observation() creates an observation from the current position
Observation create_observation(const Position& pos) {
    Observation obs;
    VisibilityInfo vi = compute_visibility(pos);

    obs.visible = vi.visible;
    obs.myPieces = vi.myPieces;
    obs.seenOpponentPieces = vi.seenOpponentPieces;
    obs.sideToMove = pos.side_to_move();
    obs.epSquares = pos.ep_squares();

    // Query castling rights for our side
    Color us = pos.side_to_move();
    CastlingRights ourCastlingRights = us == WHITE ? WHITE_CASTLING : BLACK_CASTLING;
    obs.castlingRights = pos.can_castle(ourCastlingRights) ? 1 : 0;

    obs.halfmoveClock = pos.rule50_count();
    obs.fullmoveNumber = pos.game_ply();

    return obs;
}

/// BeliefState implementation

bool BeliefState::is_consistent(const Position& pos, const Observation& obs) {
    // Check side to move
    if (pos.side_to_move() != obs.sideToMove)
        return false;

    Color us = obs.sideToMove;
    Color them = ~us;

    // Check that our pieces match exactly
    if (pos.pieces(us) != obs.myPieces)
        return false;

    // Check that visible opponent pieces match
    Bitboard theirPieces = pos.pieces(them);
    Bitboard visibleTheirPieces = theirPieces & obs.visible;
    if (visibleTheirPieces != obs.seenOpponentPieces)
        return false;

    // Check that opponent pieces are not on squares we see as empty
    Bitboard visibleEmpty = obs.visible & ~obs.myPieces & ~obs.seenOpponentPieces;
    if (theirPieces & visibleEmpty)
        return false;

    // Check en-passant consistency
    if (obs.epSquares && pos.ep_squares() != obs.epSquares)
        return false;

    // Castling rights consistency (we know our own castling rights)
    CastlingRights ourCastlingRights = us == WHITE ? WHITE_CASTLING : BLACK_CASTLING;
    int posCastling = pos.can_castle(ourCastlingRights) ? 1 : 0;
    if (posCastling != obs.castlingRights)
        return false;

    return true;
}

Observation BeliefState::parse_fog_fen(const std::string& fogFen, const Variant* variant) {
    Observation obs;

    if (fogFen.empty())
        return obs;

    std::istringstream ss(fogFen);
    std::string boardToken, stmToken, castlingToken, epToken;

    ss >> boardToken >> stmToken >> castlingToken >> epToken;
    if (!(ss >> obs.halfmoveClock))
        obs.halfmoveClock = 0;
    if (!(ss >> obs.fullmoveNumber))
        obs.fullmoveNumber = 1;

    obs.sideToMove = (stmToken == "b" ? BLACK : WHITE);

    Color us = obs.sideToMove;

    auto files = variant ? int(variant->maxFile) + 1 : FILE_NB;
    auto ranks = variant ? int(variant->maxRank) + 1 : RANK_NB;

    int rankIdx = 0;
    int fileIdx = 0;

    auto mark_square = [&](int f, int r) {
        if (f >= files || r >= ranks)
            return Square(SQ_NONE);
        return make_square(File(f), Rank(ranks - 1 - r));
    };

    for (char c : boardToken) {
        if (c == '/') {
            rankIdx++;
            fileIdx = 0;
            continue;
        }

        if (std::isdigit(static_cast<unsigned char>(c))) {
            int emptyCount = c - '0';
            for (int i = 0; i < emptyCount && fileIdx < files; ++i, ++fileIdx) {
                Square sq = mark_square(fileIdx, rankIdx);
                if (sq != SQ_NONE)
                    obs.visible |= sq;
            }
            continue;
        }

        Square sq = mark_square(fileIdx, rankIdx);
        ++fileIdx;

        if (sq == SQ_NONE)
            continue;

        if (c == '*' || c == '?')
            continue; // Unknown square, remains invisible

        obs.visible |= sq;

        if (c >= 'A' && c <= 'Z') {
            if (us == WHITE)
                obs.myPieces |= sq;
            else
                obs.seenOpponentPieces |= sq;
        } else if (c >= 'a' && c <= 'z') {
            if (us == BLACK)
                obs.myPieces |= sq;
            else
                obs.seenOpponentPieces |= sq;
        }
    }

    // Castling rights we know for our side only
    if (castlingToken != "-") {
        if (us == WHITE && (castlingToken.find('K') != std::string::npos ||
                            castlingToken.find('Q') != std::string::npos))
            obs.castlingRights = 1;
        else if (us == BLACK && (castlingToken.find('k') != std::string::npos ||
                                 castlingToken.find('q') != std::string::npos))
            obs.castlingRights = 1;
    }

    // En-passant square visibility
    if (epToken != "-") {
        if (epToken.size() >= 2) {
            File f = File(epToken[0] - 'a');
            Rank r = Rank(epToken[1] - '1');
            if (f >= FILE_A && f < FILE_NB && r >= RANK_1 && r < RANK_NB) {
                Square sq = make_square(f, r);
                obs.epSquares |= sq;
                obs.visible |= sq;
            }
        }
    }

    return obs;
}

bool BeliefState::is_king_capturable(const Position& pos) const {
    // Check if the side-to-move can capture the opponent's king
    // If so, the game would have already ended, so this state is illegal
    Color us = pos.side_to_move();
    Color them = ~us;

    Square theirKing = pos.square<KING>(them);
    if (theirKing == SQ_NONE)
        return true; // No king = capturable (illegal state)

    // Check if any of our pieces attack their king
    return bool(pos.attackers_to(theirKing) & pos.pieces(us));
}

bool BeliefState::set_position_from_fen(Position& pos, StateInfo& st, const std::string& fen) const {
    if (!variant)
        return false;

#if defined(__cpp_exceptions) || defined(__EXCEPTIONS)
    try {
        pos.set(variant, fen, isChess960, &st, owningThread);
    } catch (...) {
        return false;
    }
#else
    // Exceptions are disabled (e.g., when built with -fno-exceptions)
    pos.set(variant, fen, isChess960, &st, owningThread);
#endif

    return true;
}

void BeliefState::enumerate_candidates(const ObservationHistory& obsHist,
                                        const Position& truePos) {
    if (obsHist.empty())
        return;

    const Observation& obs = obsHist.last();
    Color us = obs.sideToMove;
    Color them = ~us;

    Bitboard visible = obs.visible;
    Bitboard myPieces = obs.myPieces;

    Bitboard oppPieces = truePos.pieces(them);
    Bitboard hiddenOpp = oppPieces & ~visible;

    std::vector<Square> hiddenSquares;
    Bitboard unseen = (~visible & ~myPieces) & AllSquares;
    Bitboard tmpUnseen = unseen;
    while (tmpUnseen) {
        hiddenSquares.push_back(pop_lsb(tmpUnseen));
    }

    std::vector<Square> hiddenOppSquares;
    std::vector<Piece> hiddenOppPieces;
    Bitboard tmpHiddenOpp = hiddenOpp;
    while (tmpHiddenOpp) {
        Square sq = pop_lsb(tmpHiddenOpp);
        hiddenOppSquares.push_back(sq);
        hiddenOppPieces.push_back(truePos.piece_on(sq));
    }

    // If there are no hidden opponent pieces, the true position is the only consistent state
    if (hiddenOppPieces.empty()) {
        stateFens.push_back(truePos.fen());
        stateKeys.insert(truePos.key());
        return;
    }

    if (hiddenSquares.size() < hiddenOppPieces.size()) {
        // Not enough locations to hide opponent pieces, keep true position only
        stateFens.push_back(truePos.fen());
        stateKeys.insert(truePos.key());
        return;
    }

    // Limit enumeration to avoid combinatorial explosion
    constexpr size_t kMaxEnumeratedStates = 1024;
    size_t generated = 0;
    std::vector<Square> assignment(hiddenOppPieces.size());

    std::string baseFen = truePos.fen();

    std::function<void(size_t, Bitboard)> dfs = [&](size_t idx, Bitboard used) {
        if (generated >= kMaxEnumeratedStates)
            return;

        if (idx == hiddenOppPieces.size()) {
            StateInfo st;
            Position candidate;
            if (!set_position_from_fen(candidate, st, baseFen))
                return;

            for (Square sq : hiddenOppSquares)
                candidate.remove_piece(sq);

            for (size_t i = 0; i < hiddenOppPieces.size(); ++i)
                candidate.put_piece(hiddenOppPieces[i], assignment[i], false, NO_PIECE);

            if (is_consistent(candidate, obs) && !is_king_capturable(candidate)) {
                StateKey key = candidate.key();
                if (stateKeys.insert(key).second) {
                    stateFens.push_back(candidate.fen());
                    ++generated;
                }
            }
            return;
        }

        for (Square sq : hiddenSquares) {
            Bitboard sqBB = square_bb(sq);
            if (used & sqBB)
                continue;
            assignment[idx] = sq;
            dfs(idx + 1, used | sqBB);
            if (generated >= kMaxEnumeratedStates)
                break;
        }
    };

    dfs(0, 0);

    if (stateFens.empty()) {
        stateFens.push_back(truePos.fen());
        stateKeys.insert(truePos.key());
    }
}

void BeliefState::filter_illegal_states() {
    if (!variant)
        return;

    auto it = stateFens.begin();
    while (it != stateFens.end()) {
        StateInfo st;
        Position tempPos;
        if (!set_position_from_fen(tempPos, st, *it)) {
            it = stateFens.erase(it);
            continue;
        }

        if (is_king_capturable(tempPos)) {
            stateKeys.erase(tempPos.key());
            it = stateFens.erase(it);
        } else {
            ++it;
        }
    }
}

void BeliefState::rebuild_from_observations(const ObservationHistory& obsHist,
                                             const Position& truePos) {
    stateFens.clear();
    stateKeys.clear();

    variant = truePos.variant();
    isChess960 = truePos.is_chess960();
    owningThread = truePos.this_thread();

    if (obsHist.empty())
        return;

    // Enumerate all candidate positions
    enumerate_candidates(obsHist, truePos);

    // Filter out illegal states
    filter_illegal_states();

    // Verify all states are consistent with latest observation
    // For now, simplified - we only have the true position
    // TODO: Implement full consistency checking with FEN parsing
}

void BeliefState::update_incrementally(const ObservationHistory& obsHist, const Position& truePos) {
    if (!variant) {
        variant = truePos.variant();
        isChess960 = truePos.is_chess960();
        owningThread = truePos.this_thread();
    }

    if (obsHist.empty()) {
        stateFens.clear();
        stateKeys.clear();
        return;
    }

    const Observation& newObs = obsHist.last();

    if (stateFens.empty()) {
        rebuild_from_observations(obsHist, truePos);
        return;
    }

    size_t beforeSize = stateFens.size();
    auto it = stateFens.begin();
    while (it != stateFens.end()) {
        StateInfo st;
        Position pos;
        if (!set_position_from_fen(pos, st, *it)) {
            it = stateFens.erase(it);
            continue;
        }

        if (!is_consistent(pos, newObs) || is_king_capturable(pos)) {
            stateKeys.erase(pos.key());
            it = stateFens.erase(it);
        } else {
            ++it;
        }
    }

    bool observationExpanded = obsHist.size() >= 2
        && (   obsHist.last().visible != obsHist.observations()[obsHist.size() - 2].visible
            || obsHist.last().seenOpponentPieces != obsHist.observations()[obsHist.size() - 2].seenOpponentPieces
            || obsHist.last().epSquares != obsHist.observations()[obsHist.size() - 2].epSquares
            || obsHist.last().castlingRights != obsHist.observations()[obsHist.size() - 2].castlingRights);

    if (stateFens.empty() || observationExpanded || stateFens.size() < beforeSize / 4) {
        rebuild_from_observations(obsHist, truePos);
    }
}

std::vector<std::string> BeliefState::sample_states(size_t n, uint64_t seed) const {
    if (stateFens.empty())
        return {};

    if (stateFens.size() <= n)
        return stateFens;

    // Random sampling without replacement
    std::vector<std::string> sampled;
    std::vector<size_t> indices(stateFens.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::mt19937_64 rng(seed);
    std::shuffle(indices.begin(), indices.end(), rng);

    for (size_t i = 0; i < n && i < indices.size(); ++i)
        sampled.push_back(stateFens[indices[i]]);

    return sampled;
}

void BeliefState::compress(size_t maxStates) {
    if (!maxStates || stateFens.size() <= maxStates)
        return;

    stateFens.resize(maxStates);
    stateKeys.clear();

    if (!variant) {
        stateFens.shrink_to_fit();
        return;
    }

    for (const auto& fen : stateFens) {
        Position pos;
        StateInfo st;
        if (set_position_from_fen(pos, st, fen))
            stateKeys.insert(pos.key());
    }

    stateFens.shrink_to_fit();
}

} // namespace FogOfWar
} // namespace Stockfish
