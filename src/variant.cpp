/*
  Fairy-Stockfish, a UCI chess variant playing engine derived from Stockfish
  Copyright (C) 2018-2022 Fabian Fichter

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

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "piece.h"
#include "variant.h"

using std::string;

namespace Stockfish {

VariantMap variants; // Global object

namespace {
    // Base variant
    Variant* variant_base() {
        Variant* v = new Variant();
        return v;
    }
    // Base for all fairy variants
    Variant* chess_variant_base() {
        Variant* v = variant_base()->init();
        v->pieceToCharTable = "PNBRQ................Kpnbrq................k";
        return v;
    }
    // Standard chess (no potions)
    Variant* chess_variant() {
        Variant* v = chess_variant_base()->init();
        v->nnueAlias = "nn-";
        v->variantTemplate = "chess";
        v->startFen = "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1";
        return v;
    }
    // Spell-chess
    Variant* spell_chess_variant() {
        Variant* v = chess_variant_base()->init();
        v->nnueAlias = "nn-";
        v->variantTemplate = "spell-chess";
        v->potions = true;
        v->potionPiece[Variant::POTION_FREEZE] = CUSTOM_PIECE_1;
        v->potionPiece[Variant::POTION_JUMP] = CUSTOM_PIECE_2;
        v->potionCooldown[Variant::POTION_FREEZE] = 3;
        v->potionCooldown[Variant::POTION_JUMP] = 3;
        v->potionDropOnOccupied = true;
        v->remove_piece(KING);
        v->add_piece(COMMONER, 'k');
        v->royalPiece = COMMONER;
        v->castlingKingPiece[WHITE] = v->castlingKingPiece[BLACK] = COMMONER;
        v->pieceToChar[make_piece(WHITE, CUSTOM_PIECE_1)] = 'F';
        v->pieceToChar[make_piece(BLACK, CUSTOM_PIECE_1)] = 'f';
        v->pieceToChar[make_piece(WHITE, CUSTOM_PIECE_2)] = 'J';
        v->pieceToChar[make_piece(BLACK, CUSTOM_PIECE_2)] = 'j';
        v->extinctionValue = -VALUE_MATE;
        v->extinctionPieceTypes = piece_set(COMMONER);
        v->extinctionPieceCount = 0;
        v->startFen = "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR[JJFFFFFjjfffff] w KQkq - 0 1";
        return v;
    }
} // namespace


/// VariantMap::init() is called at startup to initialize all predefined variants

void VariantMap::init() {
    // Add to UCI_Variant option
    add("chess", chess_variant());
    add("spell-chess", spell_chess_variant());
}

// Spell-only build: keep variant parsing entry points as no-ops for link
// compatibility with UCI VariantPath handling in other build configurations.
template <bool DoCheck>
void VariantMap::parse_istream(std::istream& file) {
    (void)file;
    // No-op: variants are hardcoded in this build.
    // If DoCheck is true, callers expect validation side effects; none apply here.
}

template <bool DoCheck>
void VariantMap::parse(const std::string& data) {
    std::istringstream ss(data);
    parse_istream<DoCheck>(ss);
}

template void VariantMap::parse_istream<true>(std::istream& file);
template void VariantMap::parse_istream<false>(std::istream& file);
template void VariantMap::parse<true>(const std::string& data);
template void VariantMap::parse<false>(const std::string& data);


// Pre-calculate derived properties
Variant* Variant::conclude() {
    // Enforce consistency to allow runtime optimizations
    if (!doubleStep)
        doubleStepRegion[WHITE] = doubleStepRegion[BLACK] = 0;
    if (!doubleStepRegion[WHITE] && !doubleStepRegion[BLACK])
        doubleStep = false;

    PieceSet originalPieceTypes = pieceTypes;
    PieceSet potionPieces = NO_PIECE_SET;
    if (potions)
        for (int idx = 0; idx < Variant::POTION_TYPE_NB; ++idx)
        {
            PieceType potion = potionPiece[idx];
            if (potion != NO_PIECE_TYPE)
            {
                pieceTypes |= piece_set(potion);
                if (!(originalPieceTypes & piece_set(potion)))
                    potionPieces |= piece_set(potion);
            }
        }

    // Determine optimizations
    bool restrictedMobility = false;
    for (PieceSet ps = pieceTypes; !restrictedMobility && ps;)
    {
        PieceType pt = pop_lsb(ps);
        if (mobilityRegion[WHITE][pt] || mobilityRegion[BLACK][pt])
          restrictedMobility = true;
    }
    PieceSet boardPieceTypes = pieceTypes & ~potionPieces;

    fastAttacks =  !(boardPieceTypes & ~(CHESS_PIECES | COMMON_FAIRY_PIECES))
                  && kingType == KING
                  && !restrictedMobility
                  && !cambodianMoves
                  && !diagonalLines;
    fastAttacks2 =  !(boardPieceTypes & ~(SHOGI_PIECES | COMMON_STEP_PIECES))
                  && kingType == KING
                  && !restrictedMobility
                  && !cambodianMoves
                  && !diagonalLines;

    // Initialize calculated NNUE properties
    nnueKing =  pieceTypes & KING ? KING
              : extinctionPieceCount == 0 && (extinctionPieceTypes & COMMONER) ? COMMONER
              : NO_PIECE_TYPE;
    // The nnueKing has to present exactly once and must not change in count
    if (nnueKing != NO_PIECE_TYPE)
    {
        // If the nnueKing is involved in promotion, count might change
        if (   ((promotionPawnTypes[WHITE] | promotionPawnTypes[BLACK]) & nnueKing)
            || ((promotionPieceTypes[WHITE] | promotionPieceTypes[BLACK]) & nnueKing)
            || std::find(std::begin(promotedPieceType), std::end(promotedPieceType), nnueKing) != std::end(promotedPieceType))
            nnueKing = NO_PIECE_TYPE;
    }
    if (nnueKing != NO_PIECE_TYPE)
    {
        std::string fenBoard = startFen.substr(0, startFen.find(' '));
        // Switch NNUE from KA to A if there is no unique piece
        if (   std::count(fenBoard.begin(), fenBoard.end(), pieceToChar[make_piece(WHITE, nnueKing)]) != 1
            || std::count(fenBoard.begin(), fenBoard.end(), pieceToChar[make_piece(BLACK, nnueKing)]) != 1)
            nnueKing = NO_PIECE_TYPE;
    }
    // We can not use popcount here yet, as the lookup tables are initialized after the variants
    int nnueSquares = (maxRank + 1) * (maxFile + 1);
    nnueUsePockets = (pieceDrops && (capturesToHand || (!mustDrop && std::bitset<64>(pieceTypes).count() != 1)))
                     || seirawanGating
                     || potions;
    int nnuePockets = nnueUsePockets ? 2 * int(maxFile + 1) : 0;
    int nnueNonDropPieceIndices = (2 * std::bitset<64>(pieceTypes).count() - (nnueKing != NO_PIECE_TYPE)) * nnueSquares;
    int nnuePieceIndices = nnueNonDropPieceIndices + 2 * (std::bitset<64>(pieceTypes).count() - (nnueKing != NO_PIECE_TYPE)) * nnuePockets;
    bool nnueHasPotions = potions;
    nnuePotionZoneIndexBase = nnueHasPotions ? nnuePieceIndices : -1;
    if (nnueHasPotions)
        nnuePieceIndices += nnueSquares * COLOR_NB * Variant::POTION_TYPE_NB;
    nnuePotionCooldownIndexBase = nnueHasPotions ? nnuePieceIndices : -1;
    if (nnueHasPotions)
        nnuePieceIndices += COLOR_NB * Variant::POTION_TYPE_NB * POTION_COOLDOWN_BITS;
    int i = 0;
    for (PieceSet ps = pieceTypes; ps;)
    {
        // Make sure that the nnueKing type gets the last index, since the NNUE architecture relies on that
        PieceType pt = lsb(ps != piece_set(nnueKing) ? ps & ~piece_set(nnueKing) : ps);
        ps ^= pt;
        assert(pt != nnueKing || !ps);

        for (Color c : { WHITE, BLACK})
        {
            pieceSquareIndex[c][make_piece(c, pt)] = 2 * i * nnueSquares;
            pieceSquareIndex[c][make_piece(~c, pt)] = (2 * i + (pt != nnueKing)) * nnueSquares;
            pieceHandIndex[c][make_piece(c, pt)] = 2 * i * nnuePockets + nnueNonDropPieceIndices;
            pieceHandIndex[c][make_piece(~c, pt)] = (2 * i + 1) * nnuePockets + nnueNonDropPieceIndices;
        }
        i++;
    }

    // Map king squares to enumeration of actually available squares.
    // E.g., for xiangqi map from 0-89 to 0-8.
    // Variants might be initialized before bitboards, so do not rely on precomputed bitboards (like SquareBB).
    // Furthermore conclude() might be called on invalid configuration during validation,
    // therefore skip proper initialization in case of invalid board size.
    int nnueKingSquare = 0;
    if (nnueKing && nnueSquares <= SQUARE_NB)
        for (Square s = SQ_A1; s < nnueSquares; ++s)
        {
            Square bitboardSquare = Square(s + s / (maxFile + 1) * (FILE_MAX - maxFile));
            if (   !mobilityRegion[WHITE][nnueKing] || !mobilityRegion[BLACK][nnueKing]
                || (mobilityRegion[WHITE][nnueKing] & make_bitboard(bitboardSquare))
                || (mobilityRegion[BLACK][nnueKing] & make_bitboard(relative_square(BLACK, bitboardSquare, maxRank))))
            {
                kingSquareIndex[s] = nnueKingSquare++ * nnuePieceIndices;
            }
        }
    else
        kingSquareIndex[SQ_A1] = nnueKingSquare++ * nnuePieceIndices;
    nnueDimensions = nnueKingSquare * nnuePieceIndices;

    // Determine maximum piece count
    std::istringstream ss(startFen);
    ss >> std::noskipws;
    unsigned char token;
    nnueMaxPieces = 0;
    while ((ss >> token) && !isspace(token))
    {
        if (pieceToChar.find(token) != std::string::npos || pieceToCharSynonyms.find(token) != std::string::npos)
            nnueMaxPieces++;
    }
    if (twoBoards)
        nnueMaxPieces *= 2;

    // For endgame evaluation to be applicable, no special win rules must apply.
    // Furthermore, rules significantly changing game mechanics also invalidate it.
    endgameEval =  endgameEval != EG_EVAL_CHESS
                 ||
                   (   endgameEval == EG_EVAL_CHESS
                    && extinctionValue == VALUE_NONE
                    && checkmateValue == -VALUE_MATE
                    && stalemateValue == VALUE_DRAW
                    && !materialCounting
                    && !(flagRegion[WHITE] || flagRegion[BLACK])
                    && !mustCapture
                    && !checkCounting
                    && !makpongRule
                    && !connectN
                    && !blastOnCapture
                    && !petrifyOnCaptureTypes
                    && !capturesToHand
                    && !twoBoards
                    && !restrictedMobility
                    && kingType == KING
                   )
                 ? endgameEval : NO_EG_EVAL;

    shogiStylePromotions = false;
    for (PieceType current: promotedPieceType)
        if (current != NO_PIECE_TYPE)
        {
            shogiStylePromotions = true;
            break;
        }

    connectDirections.clear();
    if (connectHorizontal)
    {
        connectDirections.push_back(EAST);
    }
    if (connectVertical)
    {
        connectDirections.push_back(NORTH);
    }
    if (connectDiagonal)
    {
        connectDirections.push_back(NORTH_EAST);
        connectDirections.push_back(SOUTH_EAST);
    }

    // If not a connect variant, set connectPieceTypesTrimmed to no pieces.
    // connectPieceTypesTrimmed is separated so that connectPieceTypes is left unchanged for inheritance.
    if ( !(connectRegion1[WHITE] || connectRegion1[BLACK] || connectN || connectNxN || collinearN) )
    {
          connectPieceTypesTrimmed = NO_PIECE_SET;
    }
    //Otherwise optimize to pieces actually in the game.
    else
    {
        connectPieceTypesTrimmed = connectPieceTypes & pieceTypes;
    };

    return this;
}


void VariantMap::add(std::string s, Variant* v) {
  insert(std::pair<std::string, const Variant*>(s, v->conclude()));
}

void VariantMap::clear_all() {
  for (auto const& element : *this)
      delete element.second;
  clear();
}

std::vector<std::string> VariantMap::get_keys() {
  std::vector<std::string> keys;
  for (auto const& element : *this)
      keys.push_back(element.first);
  return keys;
}

} // namespace Stockfish
