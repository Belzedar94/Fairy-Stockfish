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

#ifndef VARIANT_H_INCLUDED
#define VARIANT_H_INCLUDED

#include <bitset>
#include <set>
#include <map>
#include <vector>
#include <string>
#include <functional>
#include <sstream>
#include <iostream>
#include <cctype>

#include "types.h"
#include "bitboard.h"

namespace Stockfish {

/// Variant struct stores information needed to determine the rules of a variant.

struct Variant {
  std::string variantTemplate = "fairy";
  std::string pieceToCharTable = "-";
  int pocketSize = 0;
  Rank maxRank = RANK_8;
  File maxFile = FILE_H;
  bool chess960 = false;
  bool twoBoards = false;
  int pieceValue[PHASE_NB][PIECE_TYPE_NB] = {};
  std::string customPiece[CUSTOM_PIECES_NB] = {};
  PieceSet pieceTypes = CHESS_PIECES;
  std::string pieceToChar =  " PNBRQ" + std::string(KING - QUEEN - 1, ' ') + "K" + std::string(PIECE_TYPE_NB - KING - 1, ' ')
                           + " pnbrq" + std::string(KING - QUEEN - 1, ' ') + "k" + std::string(PIECE_TYPE_NB - KING - 1, ' ');
  std::string pieceToCharSynonyms = std::string(PIECE_NB, ' ');
  std::vector<std::string> pieceToSymbol = std::vector<std::string>(PIECE_NB, "");
  std::vector<std::string> pieceToSymbolSynonyms = std::vector<std::string>(PIECE_NB, "");
  std::map<std::string, Piece> symbolToPiece;
  std::map<std::string, PieceType> symbolToPieceType;
  std::string startFen = "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1";
  Bitboard mobilityRegion[COLOR_NB][PIECE_TYPE_NB] = {};
  Bitboard promotionRegion[COLOR_NB] = {Rank8BB, Rank1BB};
  PieceType mainPromotionPawnType[COLOR_NB] = {PAWN, PAWN};
  PieceSet promotionPawnTypes[COLOR_NB] = {piece_set(PAWN), piece_set(PAWN)};
  PieceSet promotionPieceTypes[COLOR_NB] = {piece_set(QUEEN) | ROOK | BISHOP | KNIGHT,
                                            piece_set(QUEEN) | ROOK | BISHOP | KNIGHT};
  bool sittuyinPromotion = false;
  int promotionLimit[PIECE_TYPE_NB] = {}; // 0 means unlimited
  PieceType promotedPieceType[PIECE_TYPE_NB] = {};
  bool piecePromotionOnCapture = false;
  bool mandatoryPawnPromotion = true;
  bool mandatoryPiecePromotion = false;
  bool pieceDemotion = false;
  bool blastOnCapture = false;
  PieceSet blastImmuneTypes = NO_PIECE_SET;
  PieceSet mutuallyImmuneTypes = NO_PIECE_SET;
  PieceSet petrifyOnCaptureTypes = NO_PIECE_SET;
  bool petrifyBlastPieces = false;
  bool doubleStep = true;
  Bitboard doubleStepRegion[COLOR_NB] = {Rank2BB, Rank7BB};
  Bitboard tripleStepRegion[COLOR_NB] = {};
  Bitboard enPassantRegion[COLOR_NB] = {AllSquares, AllSquares};
  PieceSet enPassantTypes[COLOR_NB] = {piece_set(PAWN), piece_set(PAWN)};
  bool castling = true;
  bool castlingDroppedPiece = false;
  File castlingKingsideFile = FILE_G;
  File castlingQueensideFile = FILE_C;
  Rank castlingRank = RANK_1;
  File castlingKingFile = FILE_E;
  PieceType castlingKingPiece[COLOR_NB] = {KING, KING};
  File castlingRookKingsideFile = FILE_MAX; // only has to match if rook is not in corner in non-960 variants
  File castlingRookQueensideFile = FILE_A; // only has to match if rook is not in corner in non-960 variants
  PieceSet castlingRookPieces[COLOR_NB] = {piece_set(ROOK), piece_set(ROOK)};
  bool oppositeCastling = false;
  PieceType kingType = KING;
  bool checking = true;
  bool dropChecks = true;
  bool mustCapture = false;
  bool mustDrop = false;
  PieceType mustDropType = ALL_PIECES;
  bool pieceDrops = false;
  bool dropLoop = false;
  bool capturesToHand = false;
  bool firstRankPawnDrops = false;
  bool promotionZonePawnDrops = false;
  EnclosingRule enclosingDrop = NO_ENCLOSING;
  Bitboard enclosingDropStart = 0;
  Bitboard dropRegion[COLOR_NB] = {AllSquares, AllSquares};
  bool sittuyinRookDrop = false;
  bool dropOppositeColoredBishop = false;
  bool dropPromoted = false;
  PieceType dropNoDoubled = NO_PIECE_TYPE;
  int dropNoDoubledCount = 1;
  bool immobilityIllegal = false;
  bool gating = false;
  WallingRule wallingRule = NO_WALLING;
  Bitboard wallingRegion[COLOR_NB] = {AllSquares, AllSquares};
  bool wallOrMove = false;
  bool seirawanGating = false;
  bool cambodianMoves = false;
  Bitboard diagonalLines = 0;
  bool pass[COLOR_NB] = {false, false};
  bool passOnStalemate[COLOR_NB] = {false, false};
  bool makpongRule = false;
  bool flyingGeneral = false;
  Rank soldierPromotionRank = RANK_1;
  EnclosingRule flipEnclosedPieces = NO_ENCLOSING;
  bool freeDrops = false;

  // game end
  PieceSet nMoveRuleTypes[COLOR_NB] = {piece_set(PAWN), piece_set(PAWN)};
  int nMoveRule = 50;
  int nFoldRule = 3;
  Value nFoldValue = VALUE_DRAW;
  bool nFoldValueAbsolute = false;
  bool perpetualCheckIllegal = false;
  bool moveRepetitionIllegal = false;
  ChasingRule chasingRule = NO_CHASING;
  Value stalemateValue = VALUE_DRAW;
  bool stalematePieceCount = false; // multiply stalemate value by sign(count(~stm) - count(stm))
  Value checkmateValue = -VALUE_MATE;
  bool shogiPawnDropMateIllegal = false;
  bool shatarMateRule = false;
  bool bikjangRule = false;
  Value extinctionValue = VALUE_NONE;
  bool extinctionClaim = false;
  bool extinctionPseudoRoyal = false;
  bool dupleCheck = false;
  PieceSet extinctionPieceTypes = NO_PIECE_SET;
  int extinctionPieceCount = 0;
  int extinctionOpponentPieceCount = 0;
  PieceType flagPiece[COLOR_NB] = {ALL_PIECES, ALL_PIECES};
  Bitboard flagRegion[COLOR_NB] = {};
  int flagPieceCount = 1;
  bool flagPieceBlockedWin = false;
  bool flagMove = false;
  bool flagPieceSafe = false;
  bool checkCounting = false;
  int connectN = 0;
  PieceSet connectPieceTypes = ~NO_PIECE_SET;
  bool connectHorizontal = true;
  bool connectVertical = true;
  bool connectDiagonal = true;
  Bitboard connectRegion1[COLOR_NB] = {};
  Bitboard connectRegion2[COLOR_NB] = {};
  int connectNxN = 0;
  int collinearN = 0;
  Value connectValue = VALUE_MATE;
  MaterialCounting materialCounting = NO_MATERIAL_COUNTING;
  bool adjudicateFullBoard = false;
  CountingRule countingRule = NO_COUNTING;
  CastlingRights castlingWins = NO_CASTLING;

  // Derived properties
  bool fastAttacks = true;
  bool fastAttacks2 = true;
  std::string nnueAlias = "";
  PieceType nnueKing = KING;
  int nnueDimensions;
  bool nnueUsePockets;
  int pieceSquareIndex[COLOR_NB][PIECE_NB];
  int pieceHandIndex[COLOR_NB][PIECE_NB];
  int kingSquareIndex[SQUARE_NB];
  int nnueMaxPieces;
  EndgameEval endgameEval = EG_EVAL_CHESS;
  bool shogiStylePromotions = false;
  std::vector<Direction> connectDirections;
  PieceSet connectPieceTypesTrimmed = ~NO_PIECE_SET;
  static bool is_piece_id_suffix(char c) {
      return c == '\'' || c == '"' || c == '!';
  }

  static bool is_piece_id_start(char c) {
      return std::isalpha(static_cast<unsigned char>(c));
  }

  static std::string normalize_piece_symbol(const std::string& token, Color c) {
      if (token.empty() || !is_piece_id_start(token[0]) || token.size() > 2)
          return "";

      std::string symbol(1, c == WHITE ? char(std::toupper(static_cast<unsigned char>(token[0])))
                                      : char(std::tolower(static_cast<unsigned char>(token[0]))));
      if (token.size() == 2) {
          if (!is_piece_id_suffix(token[1]))
              return "";
          symbol.push_back(token[1]);
      }
      return symbol;
  }

  void rebuild_piece_symbol_maps() {
      symbolToPiece.clear();
      symbolToPieceType.clear();
      for (Piece pc = W_PAWN; pc < PIECE_NB; ++pc)
      {
          if (pieceToSymbol[pc].empty() && pieceToChar[pc] != ' ')
              pieceToSymbol[pc] = std::string(1, pieceToChar[pc]);
          if (pieceToSymbolSynonyms[pc].empty() && pieceToCharSynonyms[pc] != ' ')
              pieceToSymbolSynonyms[pc] = std::string(1, pieceToCharSynonyms[pc]);
          if (!pieceToSymbol[pc].empty())
          {
              symbolToPiece[pieceToSymbol[pc]] = pc;
              std::string typeSymbol = pieceToSymbol[pc];
              typeSymbol[0] = char(std::toupper(static_cast<unsigned char>(typeSymbol[0])));
              symbolToPieceType[typeSymbol] = type_of(pc);
          }
          if (!pieceToSymbolSynonyms[pc].empty())
          {
              symbolToPiece[pieceToSymbolSynonyms[pc]] = pc;
              std::string typeSymbol = pieceToSymbolSynonyms[pc];
              typeSymbol[0] = char(std::toupper(static_cast<unsigned char>(typeSymbol[0])));
              symbolToPieceType[typeSymbol] = type_of(pc);
          }
      }
  }

  Piece piece_from_symbol(const std::string& token) const {
      auto it = symbolToPiece.find(token);
      return it == symbolToPiece.end() ? NO_PIECE : it->second;
  }

  PieceType piece_type_from_symbol(const std::string& token) const {
      if (token.empty())
          return NO_PIECE_TYPE;
      std::string typeToken = token;
      typeToken[0] = char(std::toupper(static_cast<unsigned char>(typeToken[0])));
      auto it = symbolToPieceType.find(typeToken);
      return it == symbolToPieceType.end() ? NO_PIECE_TYPE : it->second;
  }

  const std::string& piece_symbol(Piece pc) const {
      static const std::string empty;
      return pieceToSymbol[pc].empty() ? empty : pieceToSymbol[pc];
  }

  const std::string& piece_symbol_synonym(Piece pc) const {
      static const std::string empty;
      return pieceToSymbolSynonyms[pc].empty() ? empty : pieceToSymbolSynonyms[pc];
  }

  void add_piece(PieceType pt, const std::string& token, std::string betza = "", const std::string& token2 = "") {
      std::string whiteSymbol = normalize_piece_symbol(token, WHITE);
      std::string blackSymbol = normalize_piece_symbol(token, BLACK);
      std::string whiteSyn = normalize_piece_symbol(token2, WHITE);
      std::string blackSyn = normalize_piece_symbol(token2, BLACK);

      if (whiteSymbol.empty() || blackSymbol.empty())
      {
          remove_piece(pt);
          return;
      }

      auto remove_if_duplicate = [&](const std::string& symbol) {
          if (symbol.empty())
              return;
          for (Piece p = W_PAWN; p < PIECE_NB; ++p)
              if (pieceToSymbol[p] == symbol || pieceToSymbolSynonyms[p] == symbol)
                  remove_piece(type_of(p));
      };

      remove_if_duplicate(whiteSymbol);
      remove_if_duplicate(blackSymbol);
      remove_if_duplicate(whiteSyn);
      remove_if_duplicate(blackSyn);

      pieceToChar[make_piece(WHITE, pt)] = whiteSymbol[0];
      pieceToChar[make_piece(BLACK, pt)] = blackSymbol[0];
      pieceToCharSynonyms[make_piece(WHITE, pt)] = whiteSyn.empty() ? ' ' : whiteSyn[0];
      pieceToCharSynonyms[make_piece(BLACK, pt)] = blackSyn.empty() ? ' ' : blackSyn[0];
      pieceToSymbol[make_piece(WHITE, pt)] = whiteSymbol;
      pieceToSymbol[make_piece(BLACK, pt)] = blackSymbol;
      pieceToSymbolSynonyms[make_piece(WHITE, pt)] = whiteSyn;
      pieceToSymbolSynonyms[make_piece(BLACK, pt)] = blackSyn;
      pieceTypes |= pt;
      if (is_custom(pt))
          customPiece[pt - CUSTOM_PIECES] = betza;
  }

  void add_piece(PieceType pt, char c, std::string betza = "", char c2 = ' ') {
      std::string token(1, c);
      std::string token2;
      if (c2 != ' ')
          token2 = std::string(1, c2);
      add_piece(pt, token, betza, token2);
  }

  void add_piece(PieceType pt, char c, char c2) {
      add_piece(pt, c, "", c2);
  }

  void remove_piece(PieceType pt) {
      pieceToChar[make_piece(WHITE, pt)] = ' ';
      pieceToChar[make_piece(BLACK, pt)] = ' ';
      pieceToCharSynonyms[make_piece(WHITE, pt)] = ' ';
      pieceToCharSynonyms[make_piece(BLACK, pt)] = ' ';
      pieceToSymbol[make_piece(WHITE, pt)].clear();
      pieceToSymbol[make_piece(BLACK, pt)].clear();
      pieceToSymbolSynonyms[make_piece(WHITE, pt)].clear();
      pieceToSymbolSynonyms[make_piece(BLACK, pt)].clear();
      pieceTypes &= ~piece_set(pt);
      // erase from promotion types to ensure consistency
      promotionPieceTypes[WHITE] &= ~piece_set(pt);
      promotionPieceTypes[BLACK] &= ~piece_set(pt);
  }

  void reset_pieces() {
      pieceToChar = std::string(PIECE_NB, ' ');
      pieceToCharSynonyms = std::string(PIECE_NB, ' ');
      pieceToSymbol.assign(PIECE_NB, "");
      pieceToSymbolSynonyms.assign(PIECE_NB, "");
      symbolToPiece.clear();
      symbolToPieceType.clear();
      pieceTypes = NO_PIECE_SET;
      // clear promotion types to ensure consistency
      promotionPieceTypes[WHITE] = NO_PIECE_SET;
      promotionPieceTypes[BLACK] = NO_PIECE_SET;
  }

  // Reset values that always need to be redefined
  Variant* init() {
      nnueAlias = "";
      endgameEval = EG_EVAL_CHESS;
      return this;
  }

  Variant* conclude();
};

class VariantMap : public std::map<std::string, const Variant*> {
public:
  void init();
  template <bool DoCheck> void parse(std::string path);
  template <bool DoCheck> void parse_istream(std::istream& file);
  void clear_all();
  std::vector<std::string> get_keys();

private:
  void add(std::string s, Variant* v);
};

extern VariantMap variants;

} // namespace Stockfish

#endif // #ifndef VARIANT_H_INCLUDED
