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

#include <algorithm>
#include <array>
#include <bitset>
#include <map>
#include <vector>

#include "bitboard.h"
#include "magic.h"
#include "misc.h"
#include "piece.h"

namespace Stockfish {

#if defined(VERY_LARGE_BOARDS)
Bitboard AllSquares;
Bitboard DarkSquares;

Bitboard FileABB;
Bitboard FileBBB;
Bitboard FileCBB;
Bitboard FileDBB;
Bitboard FileEBB;
Bitboard FileFBB;
Bitboard FileGBB;
Bitboard FileHBB;
Bitboard FileIBB;
Bitboard FileJBB;
Bitboard FileKBB;
Bitboard FileLBB;
Bitboard FileMBB;
Bitboard FileNBB;
Bitboard FileOBB;
Bitboard FilePBB;

Bitboard Rank1BB;
Bitboard Rank2BB;
Bitboard Rank3BB;
Bitboard Rank4BB;
Bitboard Rank5BB;
Bitboard Rank6BB;
Bitboard Rank7BB;
Bitboard Rank8BB;
Bitboard Rank9BB;
Bitboard Rank10BB;
Bitboard Rank11BB;
Bitboard Rank12BB;
Bitboard Rank13BB;
Bitboard Rank14BB;
Bitboard Rank15BB;
Bitboard Rank16BB;

Bitboard QueenSide;
Bitboard CenterFiles;
Bitboard KingSide;
Bitboard Center;
Bitboard KingFlank[FILE_NB];
#endif

uint8_t PopCnt16[1 << 16];
uint8_t SquareDistance[SQUARE_NB][SQUARE_NB];

Bitboard SquareBB[SQUARE_NB];
Bitboard LineBB[SQUARE_NB][SQUARE_NB];
Bitboard BetweenBB[SQUARE_NB][SQUARE_NB];
Bitboard PseudoAttacks[COLOR_NB][PIECE_TYPE_NB][SQUARE_NB];
Bitboard PseudoMoves[2][COLOR_NB][PIECE_TYPE_NB][SQUARE_NB];
Bitboard LeaperAttacks[COLOR_NB][PIECE_TYPE_NB][SQUARE_NB];
Bitboard LeaperMoves[2][COLOR_NB][PIECE_TYPE_NB][SQUARE_NB];
Bitboard BoardSizeBB[FILE_NB][RANK_NB];
RiderType AttackRiderTypes[PIECE_TYPE_NB];
RiderType MoveRiderTypes[2][PIECE_TYPE_NB];

Magic RookMagicsH[SQUARE_NB];
Magic RookMagicsV[SQUARE_NB];
Magic BishopMagics[SQUARE_NB];
Magic CannonMagicsH[SQUARE_NB];
Magic CannonMagicsV[SQUARE_NB];
Magic LameDabbabaMagics[SQUARE_NB];
Magic HorseMagics[SQUARE_NB];
Magic ElephantMagics[SQUARE_NB];
Magic JanggiElephantMagics[SQUARE_NB];
Magic CannonDiagMagics[SQUARE_NB];
Magic NightriderMagics[SQUARE_NB];
Magic GrasshopperMagicsH[SQUARE_NB];
Magic GrasshopperMagicsV[SQUARE_NB];
Magic GrasshopperMagicsD[SQUARE_NB];

Magic* magics[] = {BishopMagics, RookMagicsH, RookMagicsV, CannonMagicsH, CannonMagicsV,
                   LameDabbabaMagics, HorseMagics, ElephantMagics, JanggiElephantMagics, CannonDiagMagics, NightriderMagics,
                   GrasshopperMagicsH, GrasshopperMagicsV, GrasshopperMagicsD};

namespace {

// Attack tables are allocated at runtime to match the board size.
  std::vector<Bitboard> RookTableH;  // To store horizontal rook attacks
  std::vector<Bitboard> RookTableV;  // To store vertical rook attacks
  std::vector<Bitboard> BishopTable; // To store bishop attacks
  std::vector<Bitboard> CannonTableH;  // To store horizontal cannon attacks
  std::vector<Bitboard> CannonTableV;  // To store vertical cannon attacks
  std::vector<Bitboard> LameDabbabaTable;  // To store lame dabbaba attacks
  std::vector<Bitboard> HorseTable;  // To store horse attacks
  std::vector<Bitboard> ElephantTable;  // To store elephant attacks
  std::vector<Bitboard> JanggiElephantTable;  // To store janggi elephant attacks
  std::vector<Bitboard> CannonDiagTable; // To store diagonal cannon attacks
  std::vector<Bitboard> NightriderTable; // To store nightrider attacks
  std::vector<Bitboard> GrasshopperTableH;  // To store horizontal grasshopper attacks
  std::vector<Bitboard> GrasshopperTableV;  // To store vertical grasshopper attacks
  std::vector<Bitboard> GrasshopperTableD; // To store diagonal grasshopper attacks

  // Rider directions
  const std::map<Direction, int> RookDirectionsV { {NORTH, 0}, {SOUTH, 0}};
  const std::map<Direction, int> RookDirectionsH { {EAST, 0}, {WEST, 0} };
  const std::map<Direction, int> BishopDirections { {NORTH_EAST, 0}, {SOUTH_EAST, 0}, {SOUTH_WEST, 0}, {NORTH_WEST, 0} };
  const std::map<Direction, int> LameDabbabaDirections { {2 * NORTH, 0}, {2 * EAST, 0}, {2 * SOUTH, 0}, {2 * WEST, 0} };
  const std::map<Direction, int> HorseDirections { {2 * SOUTH + WEST, 0}, {2 * SOUTH + EAST, 0}, {SOUTH + 2 * WEST, 0}, {SOUTH + 2 * EAST, 0},
                                                   {NORTH + 2 * WEST, 0}, {NORTH + 2 * EAST, 0}, {2 * NORTH + WEST, 0}, {2 * NORTH + EAST, 0} };
  const std::map<Direction, int> ElephantDirections { {2 * NORTH_EAST, 0}, {2 * SOUTH_EAST, 0}, {2 * SOUTH_WEST, 0}, {2 * NORTH_WEST, 0} };
  const std::map<Direction, int> JanggiElephantDirections { {NORTH + 2 * NORTH_EAST, 0}, {EAST  + 2 * NORTH_EAST, 0},
                                                            {EAST  + 2 * SOUTH_EAST, 0}, {SOUTH + 2 * SOUTH_EAST, 0},
                                                            {SOUTH + 2 * SOUTH_WEST, 0}, {WEST  + 2 * SOUTH_WEST, 0},
                                                            {WEST  + 2 * NORTH_WEST, 0}, {NORTH + 2 * NORTH_WEST, 0} };
  const std::map<Direction, int> GrasshopperDirectionsV { {NORTH, 1}, {SOUTH, 1}};
  const std::map<Direction, int> GrasshopperDirectionsH { {EAST, 1}, {WEST, 1} };
  const std::map<Direction, int> GrasshopperDirectionsD { {NORTH_EAST, 1}, {SOUTH_EAST, 1}, {SOUTH_WEST, 1}, {NORTH_WEST, 1} };

  enum MovementType { RIDER, HOPPER, LAME_LEAPER, HOPPER_RANGE };

  template <MovementType MT>
#ifdef PRECOMPUTED_MAGICS
  void init_magics(std::vector<Bitboard>& table, Magic magics[], std::map<Direction, int> directions, const Bitboard magicsInit[]);
#else
  void init_magics(std::vector<Bitboard>& table, Magic magics[], std::map<Direction, int> directions);
#endif

  template <MovementType MT>
  Bitboard sliding_attack(std::map<Direction, int> directions, Square sq, Bitboard occupied, Color c = WHITE) {
    assert(MT != LAME_LEAPER);

    Bitboard attack = 0;

    for (auto const& [d, limit] : directions)
    {
        int count = 0;
        bool hurdle = false;
        for (Square s = sq + (c == WHITE ? d : -d);
             is_ok(s) && distance(s, s - (c == WHITE ? d : -d)) <= 2;
             s += (c == WHITE ? d : -d))
        {
            if (MT != HOPPER || hurdle)
            {
                attack |= s;
                // For hoppers we consider limit == 1 as a grasshopper,
                // but limit > 1 as a limited distance hopper
                if (limit && !(MT == HOPPER_RANGE && limit == 1) && ++count >= limit)
                    break;
            }

            if (occupied & s)
            {
                if (MT == HOPPER && !hurdle)
                    hurdle = true;
                else
                    break;
            }
        }
    }

    return attack;
  }

  Bitboard lame_leaper_path(Direction d, Square s) {
    Direction dr = d > 0 ? NORTH : SOUTH;
    Direction df = (std::abs(d % NORTH) < NORTH / 2 ? d % NORTH : -(d % NORTH)) < 0 ? WEST : EAST;
    Square to = s + d;
    Bitboard b = 0;
    if (!is_ok(to) || distance(s, to) >= 4)
        return b;
    while (s != to)
    {
        int diff = std::abs(file_of(to) - file_of(s)) - std::abs(rank_of(to) - rank_of(s));
        if (diff > 0)
            s += df;
        else if (diff < 0)
            s += dr;
        else
            s += df + dr;

        if (s != to)
            b |= s;
    }
    return b;
  }

  Bitboard lame_leaper_path(std::map<Direction, int> directions, Square s) {
    Bitboard b = 0;
    for (const auto& i : directions)
        b |= lame_leaper_path(i.first, s);
    return b;
  }

  Bitboard lame_leaper_attack(std::map<Direction, int> directions, Square s, Bitboard occupied) {
    Bitboard b = 0;
    for (const auto& i : directions)
    {
        Square to = s + i.first;
        if (is_ok(to) && distance(s, to) < 4 && !(lame_leaper_path(i.first, s) & occupied))
            b |= to;
    }
    return b;
  }

#if defined(VERY_LARGE_BOARDS)
  Bitboard on_the_fly_rider_attacks_impl(RiderType R, Square s, Bitboard occupied) {
    switch (R)
    {
    case RIDER_BISHOP:
        return sliding_attack<RIDER>(BishopDirections, s, occupied);
    case RIDER_ROOK_H:
        return sliding_attack<RIDER>(RookDirectionsH, s, occupied);
    case RIDER_ROOK_V:
        return sliding_attack<RIDER>(RookDirectionsV, s, occupied);
    case RIDER_CANNON_H:
        return sliding_attack<HOPPER>(RookDirectionsH, s, occupied);
    case RIDER_CANNON_V:
        return sliding_attack<HOPPER>(RookDirectionsV, s, occupied);
    case RIDER_CANNON_DIAG:
        return sliding_attack<HOPPER>(BishopDirections, s, occupied);
    case RIDER_NIGHTRIDER:
        return sliding_attack<RIDER>(HorseDirections, s, occupied);
    case RIDER_GRASSHOPPER_H:
        return sliding_attack<HOPPER>(GrasshopperDirectionsH, s, occupied);
    case RIDER_GRASSHOPPER_V:
        return sliding_attack<HOPPER>(GrasshopperDirectionsV, s, occupied);
    case RIDER_GRASSHOPPER_D:
        return sliding_attack<HOPPER>(GrasshopperDirectionsD, s, occupied);
    case RIDER_LAME_DABBABA:
        return lame_leaper_attack(LameDabbabaDirections, s, occupied);
    case RIDER_HORSE:
        return lame_leaper_attack(HorseDirections, s, occupied);
    case RIDER_ELEPHANT:
        return lame_leaper_attack(ElephantDirections, s, occupied);
    case RIDER_JANGGI_ELEPHANT:
        return lame_leaper_attack(JanggiElephantDirections, s, occupied);
    default:
        return Bitboard(0);
    }
  }
#endif

}

/// safe_destination() returns the bitboard of target square for the given step
/// from the given square. If the step is off the board, returns empty bitboard.

inline Bitboard safe_destination(Square s, int step) {
    Square to = Square(s + step);
    return is_ok(to) && distance(s, to) <= 3 ? square_bb(to) : Bitboard(0);
}


/// Bitboards::pretty() returns an ASCII representation of a bitboard suitable
/// to be printed to standard output. Useful for debugging.

std::string Bitboards::pretty(Bitboard b) {

  std::string border;
  for (int f = 0; f < FILE_NB; ++f)
      border += "+---";
  border += "+\n";

  std::string s = border;

  for (Rank r = RANK_MAX; r >= RANK_1; --r)
  {
      for (File f = FILE_A; f <= FILE_MAX; ++f)
          s += (b & make_square(f, r)) ? "| X " : "|   ";

      std::string rankStr = std::to_string(int(r) + 1);
      if (rankStr.size() == 1)
          rankStr = " " + rankStr;
      s += "| " + rankStr + "\n";
      s += border;
  }

  s += "  ";
  for (File f = FILE_A; f <= FILE_MAX; ++f)
  {
      s += char('a' + f);
      if (f != FILE_MAX)
          s += "   ";
  }
  s += "\n";

  return s;
}

/// Bitboards::init_pieces() initializes piece move/attack bitboards and rider types

void Bitboards::init_pieces() {

  for (PieceType pt = PAWN; pt <= KING; ++pt)
  {
      const PieceInfo* pi = pieceMap.find(pt)->second;

      // Detect rider types
      for (auto modality : {MODALITY_QUIET, MODALITY_CAPTURE})
      {
          for (bool initial : {false, true})
          {
              // We do not support initial captures
              if (modality == MODALITY_CAPTURE && initial)
                  continue;
              auto& riderTypes = modality == MODALITY_CAPTURE ? AttackRiderTypes[pt] : MoveRiderTypes[initial][pt];
              riderTypes = NO_RIDER;
              for (auto const& [d, limit] : pi->steps[initial][modality])
              {
                  if (limit && LameDabbabaDirections.find(d) != LameDabbabaDirections.end())
                      riderTypes |= RIDER_LAME_DABBABA;
                  if (limit && HorseDirections.find(d) != HorseDirections.end())
                      riderTypes |= RIDER_HORSE;
                  if (limit && ElephantDirections.find(d) != ElephantDirections.end())
                      riderTypes |= RIDER_ELEPHANT;
                  if (limit && JanggiElephantDirections.find(d) != JanggiElephantDirections.end())
                      riderTypes |= RIDER_JANGGI_ELEPHANT;
              }
              for (auto const& [d, limit] : pi->slider[initial][modality])
              {
                  if (BishopDirections.find(d) != BishopDirections.end())
                      riderTypes |= RIDER_BISHOP;
                  if (RookDirectionsH.find(d) != RookDirectionsH.end())
                      riderTypes |= RIDER_ROOK_H;
                  if (RookDirectionsV.find(d) != RookDirectionsV.end())
                      riderTypes |= RIDER_ROOK_V;
                  if (HorseDirections.find(d) != HorseDirections.end())
                      riderTypes |= RIDER_NIGHTRIDER;
              }
              for (auto const& [d, limit] : pi->hopper[initial][modality])
              {
                  if (RookDirectionsH.find(d) != RookDirectionsH.end())
                      riderTypes |= limit == 1 ? RIDER_GRASSHOPPER_H : RIDER_CANNON_H;
                  if (RookDirectionsV.find(d) != RookDirectionsV.end())
                      riderTypes |= limit == 1 ? RIDER_GRASSHOPPER_V : RIDER_CANNON_V;
                  if (BishopDirections.find(d) != BishopDirections.end())
                      riderTypes |= limit == 1 ? RIDER_GRASSHOPPER_D : RIDER_CANNON_DIAG;
              }
          }
      }

      // Initialize move/attack bitboards
      for (Color c : { WHITE, BLACK })
      {
          for (Square s = SQ_A1; s <= SQ_MAX; ++s)
          {
              for (auto modality : {MODALITY_QUIET, MODALITY_CAPTURE})
              {
                  for (bool initial : {false, true})
                  {
                      // We do not support initial captures
                      if (modality == MODALITY_CAPTURE && initial)
                          continue;
                      auto& pseudo = modality == MODALITY_CAPTURE ? PseudoAttacks[c][pt][s] : PseudoMoves[initial][c][pt][s];
                      auto& leaper = modality == MODALITY_CAPTURE ? LeaperAttacks[c][pt][s] : LeaperMoves[initial][c][pt][s];
                      pseudo = 0;
                      leaper = 0;
                      for (auto const& [d, limit] : pi->steps[initial][modality])
                      {
                          pseudo |= safe_destination(s, c == WHITE ? d : -d);
                          if (!limit)
                              leaper |= safe_destination(s, c == WHITE ? d : -d);
                      }
                      pseudo |= sliding_attack<RIDER>(pi->slider[initial][modality], s, 0, c);
                      pseudo |= sliding_attack<HOPPER_RANGE>(pi->hopper[initial][modality], s, 0, c);
                  }
              }
          }
      }
  }
}


/// Bitboards::init() initializes various bitboard tables. It is called at
/// startup and relies on global objects to be already zero-initialized.

void Bitboards::init() {

  for (unsigned i = 0; i < (1 << 16); ++i)
      PopCnt16[i] = uint8_t(std::bitset<16>(i).count());

  for (Square s = SQ_A1; s <= SQ_MAX; ++s)
      SquareBB[s] = make_bitboard(s);

#if defined(VERY_LARGE_BOARDS)
  std::array<Bitboard, FILE_NB> fileMask{};
  std::array<Bitboard, RANK_NB> rankMask{};
  AllSquares = Bitboard(0);
  DarkSquares = Bitboard(0);

  for (File f = FILE_A; f <= FILE_MAX; ++f)
      for (Rank r = RANK_1; r <= RANK_MAX; ++r)
      {
          Square sq = make_square(f, r);
          Bitboard bb = SquareBB[sq];
          fileMask[f] |= bb;
          rankMask[r] |= bb;
          AllSquares |= bb;
          if (((int(f) ^ int(r)) & 1) == 1)
              DarkSquares |= bb;
      }

  std::array<Bitboard*, FILE_NB> fileGlobals = { &FileABB, &FileBBB, &FileCBB, &FileDBB,
                                                 &FileEBB, &FileFBB, &FileGBB, &FileHBB,
                                                 &FileIBB, &FileJBB, &FileKBB, &FileLBB,
                                                 &FileMBB, &FileNBB, &FileOBB, &FilePBB };
  for (int f = 0; f < FILE_NB; ++f)
      *fileGlobals[f] = fileMask[f];

  std::array<Bitboard*, RANK_NB> rankGlobals = { &Rank1BB, &Rank2BB, &Rank3BB, &Rank4BB,
                                                 &Rank5BB, &Rank6BB, &Rank7BB, &Rank8BB,
                                                 &Rank9BB, &Rank10BB, &Rank11BB, &Rank12BB,
                                                 &Rank13BB, &Rank14BB, &Rank15BB, &Rank16BB };
  for (int r = 0; r < RANK_NB; ++r)
      *rankGlobals[r] = rankMask[r];

  int halfFiles = FILE_NB / 2;
  QueenSide = Bitboard(0);
  KingSide = Bitboard(0);
  for (int f = 0; f < halfFiles; ++f)
      QueenSide |= fileMask[f];
  for (int f = halfFiles; f < FILE_NB; ++f)
      KingSide |= fileMask[f];

  int centerWidth = std::min(4, int(FILE_NB));
  int centerStart = std::max(0, halfFiles - centerWidth / 2);
  int centerEnd = std::min(int(FILE_NB) - 1, centerStart + centerWidth - 1);
  CenterFiles = Bitboard(0);
  for (int f = centerStart; f <= centerEnd; ++f)
      CenterFiles |= fileMask[f];

  int rankHalf = int(RANK_NB) / 2;
  int centerRankStart = std::max(0, rankHalf - 1);
  int centerRankEnd = std::min(int(RANK_NB) - 1, centerRankStart + 1);
  Center = Bitboard(0);
  for (int r = centerRankStart; r <= centerRankEnd; ++r)
      Center |= rankMask[r];
  Center &= CenterFiles;

  for (int f = 0; f < FILE_NB; ++f)
  {
      if (f < centerStart)
      {
          Bitboard flank = QueenSide;
          if (f == 0 && halfFiles > 0)
              flank ^= fileMask[halfFiles - 1];
          KingFlank[f] = flank;
      }
      else if (f <= centerEnd)
          KingFlank[f] = CenterFiles;
      else
      {
          Bitboard flank = KingSide;
          if (f == FILE_NB - 1 && halfFiles < FILE_NB)
              flank ^= fileMask[halfFiles];
          KingFlank[f] = flank;
      }
  }
#endif

  for (File f = FILE_A; f <= FILE_MAX; ++f)
      for (Rank r = RANK_1; r <= RANK_MAX; ++r)
          BoardSizeBB[f][r] = forward_file_bb(BLACK, make_square(f, r)) | SquareBB[make_square(f, r)] | (f > FILE_A ? BoardSizeBB[f - 1][r] : Bitboard(0));

  for (Square s1 = SQ_A1; s1 <= SQ_MAX; ++s1)
      for (Square s2 = SQ_A1; s2 <= SQ_MAX; ++s2)
              SquareDistance[s1][s2] = std::max(distance<File>(s1, s2), distance<Rank>(s1, s2));

#if !defined(VERY_LARGE_BOARDS)
#ifdef PRECOMPUTED_MAGICS
  init_magics<RIDER>(RookTableH, RookMagicsH, RookDirectionsH, RookMagicHInit);
  init_magics<RIDER>(RookTableV, RookMagicsV, RookDirectionsV, RookMagicVInit);
  init_magics<RIDER>(BishopTable, BishopMagics, BishopDirections, BishopMagicInit);
  init_magics<HOPPER>(CannonTableH, CannonMagicsH, RookDirectionsH, CannonMagicHInit);
  init_magics<HOPPER>(CannonTableV, CannonMagicsV, RookDirectionsV, CannonMagicVInit);
  init_magics<LAME_LEAPER>(LameDabbabaTable, LameDabbabaMagics, LameDabbabaDirections, LameDabbabaMagicInit);
  init_magics<LAME_LEAPER>(HorseTable, HorseMagics, HorseDirections, HorseMagicInit);
  init_magics<LAME_LEAPER>(ElephantTable, ElephantMagics, ElephantDirections, ElephantMagicInit);
  init_magics<LAME_LEAPER>(JanggiElephantTable, JanggiElephantMagics, JanggiElephantDirections, JanggiElephantMagicInit);
  init_magics<HOPPER>(CannonDiagTable, CannonDiagMagics, BishopDirections, CannonDiagMagicInit);
  init_magics<RIDER>(NightriderTable, NightriderMagics, HorseDirections, NightriderMagicInit);
  init_magics<HOPPER>(GrasshopperTableH, GrasshopperMagicsH, GrasshopperDirectionsH, GrasshopperMagicHInit);
  init_magics<HOPPER>(GrasshopperTableV, GrasshopperMagicsV, GrasshopperDirectionsV, GrasshopperMagicVInit);
  init_magics<HOPPER>(GrasshopperTableD, GrasshopperMagicsD, GrasshopperDirectionsD, GrasshopperMagicDInit);
#else
  init_magics<RIDER>(RookTableH, RookMagicsH, RookDirectionsH);
  init_magics<RIDER>(RookTableV, RookMagicsV, RookDirectionsV);
  init_magics<RIDER>(BishopTable, BishopMagics, BishopDirections);
  init_magics<HOPPER>(CannonTableH, CannonMagicsH, RookDirectionsH);
  init_magics<HOPPER>(CannonTableV, CannonMagicsV, RookDirectionsV);
  init_magics<LAME_LEAPER>(LameDabbabaTable, LameDabbabaMagics, LameDabbabaDirections);
  init_magics<LAME_LEAPER>(HorseTable, HorseMagics, HorseDirections);
  init_magics<LAME_LEAPER>(ElephantTable, ElephantMagics, ElephantDirections);
  init_magics<LAME_LEAPER>(JanggiElephantTable, JanggiElephantMagics, JanggiElephantDirections);
  init_magics<HOPPER>(CannonDiagTable, CannonDiagMagics, BishopDirections);
  init_magics<RIDER>(NightriderTable, NightriderMagics, HorseDirections);
  init_magics<HOPPER>(GrasshopperTableH, GrasshopperMagicsH, GrasshopperDirectionsH);
  init_magics<HOPPER>(GrasshopperTableV, GrasshopperMagicsV, GrasshopperDirectionsV);
  init_magics<HOPPER>(GrasshopperTableD, GrasshopperMagicsD, GrasshopperDirectionsD);
#endif
#endif

  init_pieces();

  for (Square s1 = SQ_A1; s1 <= SQ_MAX; ++s1)
  {
      for (PieceType pt : { BISHOP, ROOK })
          for (Square s2 = SQ_A1; s2 <= SQ_MAX; ++s2)
          {
              if (PseudoAttacks[WHITE][pt][s1] & s2)
              {
                  LineBB[s1][s2]    = (attacks_bb(WHITE, pt, s1, 0) & attacks_bb(WHITE, pt, s2, 0)) | s1 | s2;
                  BetweenBB[s1][s2] = (attacks_bb(WHITE, pt, s1, square_bb(s2)) & attacks_bb(WHITE, pt, s2, square_bb(s1)));
              }
              BetweenBB[s1][s2] |= s2;
          }
  }
}

namespace {

  // init_magics() computes all rook and bishop attacks at startup. Magic
  // bitboards are used to look up attacks of sliding pieces. As a reference see
  // www.chessprogramming.org/Magic_Bitboards. In particular, here we use the so
  // called "fancy" approach.

  template <MovementType MT>
#ifdef PRECOMPUTED_MAGICS
  void init_magics(std::vector<Bitboard>& table, Magic magics[], std::map<Direction, int> directions, const Bitboard magicsInit[]) {
#else
  void init_magics(std::vector<Bitboard>& table, Magic magics[], std::map<Direction, int> directions) {
#endif

    // Optimal PRNG seeds to pick the correct magics in the shortest time
#ifndef PRECOMPUTED_MAGICS
#ifdef LARGEBOARDS
    int seeds[][RANK_NB] = { { 734, 10316, 55013, 32803, 12281, 15100,  16645, 255, 346, 89123 },
                             { 734, 10316, 55013, 32803, 12281, 15100,  16645, 255, 346, 89123 } };
#else
    int seeds[][RANK_NB] = { { 8977, 44560, 54343, 38998,  5731, 95205, 104912, 17020 },
                             {  728, 10316, 55013, 32803, 12281, 15100,  16645,   255 } };
#endif
#endif

    std::vector<size_t> offsets(SQUARE_NB);
    size_t totalSize = 0;

    for (Square s = SQ_A1; s <= SQ_MAX; ++s)
    {
        Bitboard edges = ((Rank1BB | rank_bb(RANK_MAX)) & ~rank_bb(s)) | ((FileABB | file_bb(FILE_MAX)) & ~file_bb(s));
        Magic& m = magics[s];
        m.mask = (MT == LAME_LEAPER ? lame_leaper_path(directions, s)
                                    : sliding_attack<MT == HOPPER ? HOPPER_RANGE : MT>(directions, s, 0)) & ~edges;
#ifdef LARGEBOARDS
        m.shift = 128 - popcount(m.mask);
#else
        m.shift = (Is64Bit ? 64 : 32) - popcount(m.mask);
#endif
        offsets[s] = totalSize;
        totalSize += size_t(1) << popcount(m.mask);
    }

    table.assign(totalSize, Bitboard(0));

    const size_t maxSubsets = size_t(1) << (FILE_NB + RANK_NB - 4);
    std::vector<Bitboard> occupancy(maxSubsets);
    std::vector<Bitboard> reference(maxSubsets);
    std::vector<int> epoch(maxSubsets, 0);
    int cnt = 0;

    for (Square s = SQ_A1; s <= SQ_MAX; ++s)
    {
        Magic& m = magics[s];
        Bitboard mask = m.mask;
        m.attacks = table.data() + offsets[s];

        Bitboard b = 0;
        size_t index = 0;
        do {
            occupancy[index] = b;
            reference[index] = MT == LAME_LEAPER ? lame_leaper_attack(directions, s, b)
                                                 : sliding_attack<MT>(directions, s, b);
#if defined(USE_PEXT)
            if (HasPext)
                m.attacks[m.index(b)] = reference[index];
#endif
            ++index;
            b = (b - mask) & mask;
        } while (b);

#if defined(USE_PEXT)
        if (HasPext)
            continue;
#endif

        int subsetCountInt = int(index);

#ifndef PRECOMPUTED_MAGICS
        PRNG rng(seeds[Is64Bit][rank_of(s)]);
#endif

        for (int attempt = 0; attempt < subsetCountInt; )
        {
            for (m.magic = Bitboard(0); popcount((m.magic * mask) >> (SQUARE_NB - FILE_NB)) < FILE_NB - 2; )
            {
#ifdef PRECOMPUTED_MAGICS
                m.magic = magicsInit[s];
#elif defined(LARGEBOARDS)
                m.magic = (rng.sparse_rand<Bitboard>() << 64) ^ rng.sparse_rand<Bitboard>();
#else
                m.magic = rng.sparse_rand<Bitboard>();
#endif
            }

            for (++cnt, attempt = 0; attempt < subsetCountInt; ++attempt)
            {
                unsigned idx = m.index(occupancy[attempt]);

                if (epoch[idx] < cnt)
                {
                    epoch[idx] = cnt;
                    m.attacks[idx] = reference[attempt];
                }
                else if (m.attacks[idx] != reference[attempt])
                    break;
            }
        }
    }
  }
}

#if defined(VERY_LARGE_BOARDS)
Bitboard on_the_fly_rider_attacks(RiderType r, Square s, Bitboard occupied) {
  return on_the_fly_rider_attacks_impl(r, s, occupied);
}
#endif

} // namespace Stockfish
