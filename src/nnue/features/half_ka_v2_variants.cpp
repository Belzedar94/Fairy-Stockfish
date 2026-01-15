/*
  Stockfish, a UCI chess playing engine derived from Glaurung 2.1
  Copyright (C) 2004-2022 The Stockfish developers (see AUTHORS file)

  Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  Stockfish is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

//Definition of input features HalfKAv2 of NNUE evaluation function

#include "half_ka_v2_variants.h"

#include "../../position.h"

namespace Stockfish::Eval::NNUE::Features {

  // Map square to numbering on variant board
  inline Square to_variant_square(Square s, const Position& pos) {
    return Square(s - rank_of(s) * (FILE_MAX - pos.max_file()));
  }

  // Orient a square according to perspective (rotates by 180 for black)
  // Missing kings map to index 0 (SQ_A1)
  inline Square HalfKAv2Variants::orient(Color perspective, Square s, const Position& pos) {
    return s != SQ_NONE ? to_variant_square(  perspective == WHITE || (pos.flag_region(BLACK) & Rank8BB) ? s
                                            : flip_rank(s, pos.max_rank()), pos) : SQ_A1;
  }

  // Index of a feature for a given king position and another piece on some square
  inline IndexType HalfKAv2Variants::make_index(Color perspective, Square s, Piece pc, Square ksq, const Position& pos) {
    return IndexType(orient(perspective, s, pos) + pos.nnue_piece_square_index(perspective, pc) + pos.nnue_king_square_index(ksq));
  }

  // Index of a feature for a given king position and another piece on some square
  inline IndexType HalfKAv2Variants::make_index(Color perspective, int handCount, Piece pc, Square ksq, const Position& pos) {
    return IndexType(handCount + pos.nnue_piece_hand_index(perspective, pc) + pos.nnue_king_square_index(ksq));
  }

  inline IndexType HalfKAv2Variants::make_points_index(Color perspective, int plane, Square ksq, const Position& pos) {
    (void)perspective;
    return IndexType(plane + pos.nnue_points_index_base() + pos.nnue_king_square_index(ksq));
  }

  // Get a list of indices for active features
  void HalfKAv2Variants::append_active_indices(
    const Position& pos,
    Color perspective,
    ValueListInserter<IndexType> active
  ) {
    Square oriented_ksq = orient(perspective, pos.nnue_king_square(perspective), pos);
    Bitboard bb = pos.pieces(WHITE) | pos.pieces(BLACK);
    while (bb)
    {
      Square s = pop_lsb(bb);
      active.push_back(make_index(perspective, s, pos.piece_on(s), oriented_ksq, pos));
    }

    // Indices for pieces in hand
    if (pos.nnue_use_pockets())
      for (Color c : {WHITE, BLACK})
          for (PieceSet ps = pos.piece_types(); ps;)
          {
              PieceType pt = pop_lsb(ps);
              for (int i = 0; i < pos.count_in_hand(c, pt); i++)
                  active.push_back(make_index(perspective, i, make_piece(c, pt), oriented_ksq, pos));
          }

    if (pos.nnue_points_index_base() >= 0)
    {
      int planeOffset = 0;
      if (pos.nnue_points_score_planes())
      {
        int usScore = pos.points_score_clamped(perspective);
        int themScore = pos.points_score_clamped(~perspective);
        for (int bit = 0; bit < POINTS_SCORE_BITS; ++bit)
        {
          int mask = 1 << bit;
          if (usScore & mask)
            active.push_back(make_points_index(perspective, planeOffset + bit, oriented_ksq, pos));
          if (themScore & mask)
            active.push_back(make_points_index(perspective, planeOffset + POINTS_SCORE_BITS + bit, oriented_ksq, pos));
        }
        planeOffset += pos.nnue_points_score_planes();
      }
      if (pos.nnue_points_check_planes())
      {
        int usChecks = std::min<int>(std::max(0, int(pos.checks_remaining(perspective))), CHECKS_MAX);
        int themChecks = std::min<int>(std::max(0, int(pos.checks_remaining(~perspective))), CHECKS_MAX);
        for (int bit = 0; bit < CHECKS_BITS; ++bit)
        {
          int mask = 1 << bit;
          if (usChecks & mask)
            active.push_back(make_points_index(perspective, planeOffset + bit, oriented_ksq, pos));
          if (themChecks & mask)
            active.push_back(make_points_index(perspective, planeOffset + CHECKS_BITS + bit, oriented_ksq, pos));
        }
      }
    }

  }

  // append_changed_indices() : get a list of indices for recently changed features

  void HalfKAv2Variants::append_changed_indices(
    Square ksq,
    StateInfo* st,
    Color perspective,
    ValueListInserter<IndexType> removed,
    ValueListInserter<IndexType> added,
    const Position& pos
  ) {
    const auto& dp = st->dirtyPiece;
    Square oriented_ksq = orient(perspective, ksq, pos);
    for (int i = 0; i < dp.dirty_num; ++i) {
      Piece pc = dp.piece[i];
      if (dp.from[i] != SQ_NONE)
        removed.push_back(make_index(perspective, dp.from[i], pc, oriented_ksq, pos));
      else if (dp.handPiece[i] != NO_PIECE)
        removed.push_back(make_index(perspective, dp.handCount[i] - 1, dp.handPiece[i], oriented_ksq, pos));
      if (dp.to[i] != SQ_NONE)
        added.push_back(make_index(perspective, dp.to[i], pc, oriented_ksq, pos));
      else if (dp.handPiece[i] != NO_PIECE)
        added.push_back(make_index(perspective, dp.handCount[i] - 1, dp.handPiece[i], oriented_ksq, pos));
    }

    if (pos.nnue_points_index_base() >= 0)
    {
      auto clamp_points = [](int value) {
        return std::max(0, std::min(value, POINTS_SCORE_MAX));
      };
      auto clamp_checks = [](int value) {
        return std::max(0, std::min(value, CHECKS_MAX));
      };
      auto add_changed_bits = [&](int oldValue, int newValue, int baseOffset, int bits) {
        for (int bit = 0; bit < bits; ++bit)
        {
          int mask = 1 << bit;
          if ((oldValue ^ newValue) & mask)
          {
            if (oldValue & mask)
              removed.push_back(make_points_index(perspective, baseOffset + bit, oriented_ksq, pos));
            else
              added.push_back(make_points_index(perspective, baseOffset + bit, oriented_ksq, pos));
          }
        }
      };

      int planeOffset = 0;
      if (pos.nnue_points_score_planes())
      {
        int oldUs = st->previous ? clamp_points(st->previous->pointsScore[perspective]) : 0;
        int newUs = clamp_points(st->pointsScore[perspective]);
        int oldThem = st->previous ? clamp_points(st->previous->pointsScore[~perspective]) : 0;
        int newThem = clamp_points(st->pointsScore[~perspective]);
        add_changed_bits(oldUs, newUs, planeOffset, POINTS_SCORE_BITS);
        add_changed_bits(oldThem, newThem, planeOffset + POINTS_SCORE_BITS, POINTS_SCORE_BITS);
        planeOffset += pos.nnue_points_score_planes();
      }
      if (pos.nnue_points_check_planes())
      {
        int oldUs = st->previous ? clamp_checks(st->previous->checksRemaining[perspective]) : 0;
        int newUs = clamp_checks(st->checksRemaining[perspective]);
        int oldThem = st->previous ? clamp_checks(st->previous->checksRemaining[~perspective]) : 0;
        int newThem = clamp_checks(st->checksRemaining[~perspective]);
        add_changed_bits(oldUs, newUs, planeOffset, CHECKS_BITS);
        add_changed_bits(oldThem, newThem, planeOffset + CHECKS_BITS, CHECKS_BITS);
      }
    }
  }

  int HalfKAv2Variants::update_cost(StateInfo* st) {
    int cost = st->dirtyPiece.dirty_num;
    if (!currentNnueVariant || currentNnueVariant->nnuePointsIndexBase < 0)
      return cost;

    auto clamp_points = [](int value) {
      return std::max(0, std::min(value, POINTS_SCORE_MAX));
    };
    auto clamp_checks = [](int value) {
      return std::max(0, std::min(value, CHECKS_MAX));
    };
    auto prev = st->previous;
    if (currentNnueVariant->nnuePointsScorePlanes)
    {
      int oldW = prev ? clamp_points(prev->pointsScore[WHITE]) : 0;
      int newW = clamp_points(st->pointsScore[WHITE]);
      int oldB = prev ? clamp_points(prev->pointsScore[BLACK]) : 0;
      int newB = clamp_points(st->pointsScore[BLACK]);
      cost += popcount(Bitboard(oldW ^ newW)) + popcount(Bitboard(oldB ^ newB));
    }
    if (currentNnueVariant->nnuePointsCheckPlanes)
    {
      int oldW = prev ? clamp_checks(prev->checksRemaining[WHITE]) : 0;
      int newW = clamp_checks(st->checksRemaining[WHITE]);
      int oldB = prev ? clamp_checks(prev->checksRemaining[BLACK]) : 0;
      int newB = clamp_checks(st->checksRemaining[BLACK]);
      cost += popcount(Bitboard(oldW ^ newW)) + popcount(Bitboard(oldB ^ newB));
    }
    return cost;
  }

  int HalfKAv2Variants::refresh_cost(const Position& pos) {
    int cost = pos.count<ALL_PIECES>();
    if (pos.nnue_points_index_base() >= 0)
    {
      if (pos.nnue_points_score_planes())
        cost += popcount(Bitboard(pos.points_score_clamped(WHITE)))
             + popcount(Bitboard(pos.points_score_clamped(BLACK)));
      if (pos.nnue_points_check_planes())
      {
        int checksW = std::min<int>(std::max(0, int(pos.checks_remaining(WHITE))), CHECKS_MAX);
        int checksB = std::min<int>(std::max(0, int(pos.checks_remaining(BLACK))), CHECKS_MAX);
        cost += popcount(Bitboard(checksW)) + popcount(Bitboard(checksB));
      }
    }
    return cost;
  }

  bool HalfKAv2Variants::requires_refresh(StateInfo* st, Color perspective, const Position& pos) {
    return st->dirtyPiece.piece[0] == make_piece(perspective, pos.nnue_king()) || pos.flip_enclosed_pieces();
  }

}  // namespace Stockfish::Eval::NNUE::Features
