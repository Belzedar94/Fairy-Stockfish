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

#ifndef NNUE_FEATURES_HALF_KA_V2_VARIANTS_H_INCLUDED
#define NNUE_FEATURES_HALF_KA_V2_VARIANTS_H_INCLUDED

#include "../nnue_common.h"

#include "../../evaluate.h"
#include "../../misc.h"

namespace Stockfish {
  struct StateInfo;
}

namespace Stockfish::Eval::NNUE::Features {

  // Feature HalfKAv2: Combination of the position of own king
  // and the position of pieces
  class HalfKAv2Variants {

    // Orient a square according to perspective (rotates by 180 for black)
    static Square orient(Color perspective, Square s, const Position& pos);

    // Index of a feature for a given king position and another piece on some square
    static IndexType make_index(Color perspective, Square s, Piece pc, Square ksq, const Position& pos);

    // Index of a feature for a given king position and another piece in hand   
    static IndexType make_index(Color perspective, int handCount, Piece pc, Square ksq, const Position& pos);

    // Index of a feature for a given king position and potion zone square
    static IndexType make_potion_zone_index(Color perspective, Color potionColor,
                                            Variant::PotionType potion, Square s,
                                            Square ksq, const Position& pos);

    // Index of a feature for a given king position and potion cooldown bit
    static IndexType make_potion_cooldown_index(Color perspective, Color potionColor,
                                                Variant::PotionType potion, int bit,
                                                Square ksq, const Position& pos);

   public:
    // Feature name
    static constexpr const char* Name = "HalfKAv2(Friend)";

    // Hash value embedded in the evaluation file
    static constexpr std::uint32_t HashValueNoPotions = 0x5f234cb8u;
    static constexpr std::uint32_t HashValueWithPotions = 0x6a8f3c12u;

    static std::uint32_t get_hash_value() {
      return currentNnueVariant && currentNnueVariant->nnuePotionZoneIndexBase >= 0
             ? HashValueWithPotions
             : HashValueNoPotions;
    }

    // Number of feature dimensions
    static constexpr IndexType PotionZonePlanes =
        COLOR_NB * Variant::POTION_TYPE_NB;
    static constexpr IndexType PotionCooldownFeatures =
        COLOR_NB * Variant::POTION_TYPE_NB * POTION_COOLDOWN_BITS;
    static constexpr IndexType Dimensions =
        static_cast<IndexType>(SQUARE_NB)
        * (static_cast<IndexType>(SQUARE_NB) * (19 + PotionZonePlanes)
           + PotionCooldownFeatures);

    static IndexType get_dimensions() {
      return currentNnueVariant->nnueDimensions;
    }

    // Maximum number of simultaneously active features.
    static constexpr IndexType MaxActiveDimensions =
        4 * static_cast<IndexType>(SQUARE_NB);

    // Get a list of indices for active features
    static void append_active_indices(
      const Position& pos,
      Color perspective,
      ValueListInserter<IndexType> active);

    // Get a list of indices for recently changed features
    static void append_changed_indices(
      Square ksq,
      StateInfo* st,
      Color perspective,
      ValueListInserter<IndexType> removed,
      ValueListInserter<IndexType> added,
      const Position& pos);

    // Returns the cost of updating one perspective, the most costly one.
    // Assumes no refresh needed.
    static int update_cost(StateInfo* st);
    static int refresh_cost(const Position& pos);

    // Returns whether the change stored in this StateInfo means that
    // a full accumulator refresh is required.
    static bool requires_refresh(StateInfo* st, Color perspective, const Position& pos);
  };

}  // namespace Stockfish::Eval::NNUE::Features

#endif // #ifndef NNUE_FEATURES_HALF_KA_V2_VARIANTS_H_INCLUDED
