/*
  KOTH-Stockfish, a King of the Hill specialization of Stockfish
  Copyright (C) 2026 The KOTH-Stockfish developers

  KOTH-Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
*/

#ifndef KOTH_H_INCLUDED
#define KOTH_H_INCLUDED

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "types.h"

namespace Stockfish {

class Position;

namespace Koth {

enum class TransitionKind : std::uint8_t {
    ROOT,
    ACCEPTED_MOVE,
    SEARCH_NULL
};

struct TransitionEvent {
    Move           move;
    TransitionKind kind;
    bool           kingEnteredHill;
};

enum class TerminalPredicate : std::uint16_t {
    CHECKMATE              = 1u << 0,
    HILL                   = 1u << 1,
    STALEMATE              = 1u << 2,
    REPETITION3_AUTO       = 1u << 3,
    REPETITION5_DIAGNOSTIC = 1u << 4,
    RULE50_AUTO            = 1u << 5
};

enum class PrimaryReason : std::uint8_t {
    NONE,
    CHECKMATE,
    HILL,
    STALEMATE,
    AUTOMATIC_DRAW
};

enum class AdjudicationContext : std::uint8_t {
    ACCEPTED_TRAJECTORY,
    RAW_ROOT_VALIDATION,
    SEARCH_NULL,
    DIAGNOSTIC_PHYSICAL_REPLAY
};

struct GameStatus {
    std::uint16_t        predicates = 0;
    PrimaryReason        primary    = PrimaryReason::NONE;
    std::optional<Color> winner;

    constexpr bool terminal() const { return primary != PrimaryReason::NONE; }
    constexpr bool has(TerminalPredicate predicate) const {
        return predicates & static_cast<std::uint16_t>(predicate);
    }
};

constexpr Bitboard HillSquares =
  (Bitboard(1) << SQ_D4) | (Bitboard(1) << SQ_E4) | (Bitboard(1) << SQ_D5) | (Bitboard(1) << SQ_E5);

constexpr bool is_hill(Square square) {
    return is_ok(square) && (HillSquares & (Bitboard(1) << square));
}

bool any_king_on_hill(const Position& position);
bool is_canonical_uci_move(std::string_view token);

GameStatus classify(const Position&     position,
                    AdjudicationContext context = AdjudicationContext::ACCEPTED_TRAJECTORY);

std::vector<Move> physical_moves(const Position& position);
std::vector<Move> game_moves(const Position& position);
bool              is_immediate_hill_move(const Position& position, Move move);
std::vector<Move> immediate_hill_moves(const Position& position);

std::string primary_reason_name(PrimaryReason reason);
std::string serialize(const GameStatus& status);
std::string state_selftest();

}  // namespace Koth
}  // namespace Stockfish

#endif  // KOTH_H_INCLUDED
