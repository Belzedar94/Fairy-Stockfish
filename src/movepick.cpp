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
#include <cassert>

#include "movepick.h"

namespace Stockfish {

// Since continuation history grows quadratically with the number of piece types,
// we need to reserve a limited number of slots and map piece types to these slots
// in order to reduce memory consumption to a reasonable level.
int history_slot(Piece pc) {
    return pc == NO_PIECE ? 0 : (type_of(pc) == KING ? PIECE_SLOTS - 1 : type_of(pc) % (PIECE_SLOTS - 1)) + color_of(pc) * PIECE_SLOTS;
}

namespace {

  bool is_battle_kings(const Position& pos) {
    return  pos.gating()
        && pos.gating_piece_after(WHITE, PAWN)   == KNIGHT
        && pos.gating_piece_after(WHITE, KNIGHT) == BISHOP
        && pos.gating_piece_after(WHITE, BISHOP) == ROOK
        && pos.gating_piece_after(WHITE, ROOK)   == QUEEN
        && pos.gating_piece_after(WHITE, QUEEN)  == COMMONER
        && pos.gating_piece_after(BLACK, PAWN)   == KNIGHT
        && pos.gating_piece_after(BLACK, KNIGHT) == BISHOP
        && pos.gating_piece_after(BLACK, BISHOP) == ROOK
        && pos.gating_piece_after(BLACK, ROOK)   == QUEEN
        && pos.gating_piece_after(BLACK, QUEEN)  == COMMONER;
  }

  int battle_kings_stage_index(PieceType pt) {
    switch (pt)
    {
    case PAWN:     return 0;
    case KNIGHT:   return 1;
    case BISHOP:   return 2;
    case ROOK:     return 3;
    case QUEEN:    return 4;
    case COMMONER:
    case KING:     return 5;
    default:       return -1;
    }
  }

  int battle_kings_nearest_royal_distance(const Position& pos, Color c, Square sq) {
    Bitboard royals = pos.pieces(c, KING) | pos.pieces(c, COMMONER);
    if (!royals)
        return 0;

    const int boardSpan = std::max(int(pos.max_rank()) + 1, int(pos.max_file()) + 1);
    int best = boardSpan;

    while (royals)
    {
        Square ks = pop_lsb(royals);
        best = std::min(best, distance(ks, sq));
    }

    return best;
  }

  int battle_kings_adjustment(const Position& pos, Move m) {
    static constexpr int MoverBase[]   = { 720, 520, 280, -560, -1680, -5040 };
    static constexpr int SpawnBase[]   = { 640, 420, 260, -660, -2160, -5760 };
    static constexpr int CaptureBase[] = { 480, 760, 1100, 1560, 2120, 5600 };

    const Color us = pos.side_to_move();
    const Color them = ~us;
    const Square from = from_sq(m);
    const Square to = to_sq(m);

    int bonus = 0;

    const PieceType mover = type_of(pos.moved_piece(m));
    const int moverStage = battle_kings_stage_index(mover);
    if (moverStage >= 0)
        bonus += MoverBase[moverStage];

    const int boardHeight = int(pos.max_rank()) + 1;
    const int relFrom = int(relative_rank(us, from, pos.max_rank()));
    const int relTo   = int(relative_rank(us, to,   pos.max_rank()));
    const int forward = relTo - relFrom;
    const int fileCentrality = edge_distance(file_of(from), pos.max_file());

    if (moverStage >= 0)
    {
        if (moverStage <= 2)
            bonus += 85 * forward + 40 * fileCentrality;
        else
            bonus -= 75 * forward + 35 * fileCentrality;
    }

    const int lightInventory = pos.count<PAWN>(us) + pos.count<KNIGHT>(us) + pos.count<BISHOP>(us);
    const int heavyInventory = pos.count<ROOK>(us) + pos.count<QUEEN>(us);
    const int royalInventory = pos.count<COMMONER>(us) + pos.count<KING>(us);
    const int enemyRoyalInventory = pos.count<COMMONER>(them) + pos.count<KING>(them);

    if (PieceType spawn = gating_type(m); spawn != NO_PIECE_TYPE)
    {
        const int spawnStage = battle_kings_stage_index(spawn);
        bonus += SpawnBase[spawnStage];

        const int halfBoundary = (boardHeight - 1) / 2;
        const int homeDepth  = std::max(0, halfBoundary - relFrom + 1);
        const int enemyDepth = std::max(0, relFrom - halfBoundary);

        if (spawnStage >= 3)
            bonus += -170 * homeDepth + 120 * enemyDepth;
        else
            bonus += 70 * relFrom + 35 * fileCentrality;

        if (spawnStage <= 2 && heavyInventory > lightInventory)
            bonus += 60 * std::min(heavyInventory - lightInventory, 4);

        if (spawnStage >= 3 && heavyInventory >= lightInventory)
            bonus -= 80 * std::min(heavyInventory - lightInventory + 1, 5);

        const int balance = popcount(pos.attackers_to(from, us)) - popcount(pos.attackers_to(from, them));
        bonus += (spawnStage >= 3 ? 200 : 140) * balance;

        if (spawn == COMMONER)
        {
            bonus -= 1800 * (royalInventory - enemyRoyalInventory + 1);

            const int nearest = battle_kings_nearest_royal_distance(pos, us, from);
            if (nearest && nearest < 3)
                bonus -= 500 * (3 - nearest);
        }
        else
        {
            Bitboard mobility = pos.attacks_from(us, spawn, from) & ~pos.pieces(us);
            const int mobilityCount = popcount(mobility);
            bonus += 90 * mobilityCount;

            if (!mobilityCount)
                bonus -= spawnStage >= 3 ? 700 : 260;
        }

        if (moverStage >= 0)
        {
            if (spawnStage > moverStage)
                bonus += 80 * (spawnStage - moverStage);
            else
                bonus -= 60 * (moverStage - spawnStage);
        }
    }
    else if (mover == KING)
    {
        const int safety = popcount(pos.attackers_to(to, us)) - popcount(pos.attackers_to(to, them));
        bonus += 120 * safety;
    }

    if (pos.capture(m))
    {
        if (Piece captured = pos.piece_on(to); captured != NO_PIECE)
        {
            const PieceType victim = type_of(captured);
            const int victimStage = battle_kings_stage_index(victim);

            if (victimStage >= 0)
                bonus += CaptureBase[victimStage];

            if (victim == COMMONER || victim == KING)
                bonus += 5000;
            else if (victim == QUEEN)
                bonus += 2000;
        }
    }

    return bonus;
  }

  enum Stages {
    MAIN_TT, CAPTURE_INIT, GOOD_CAPTURE, REFUTATION, QUIET_INIT, QUIET, BAD_CAPTURE,
    EVASION_TT, EVASION_INIT, EVASION,
    PROBCUT_TT, PROBCUT_INIT, PROBCUT,
    QSEARCH_TT, QCAPTURE_INIT, QCAPTURE, QCHECK_INIT, QCHECK
  };

  // partial_insertion_sort() sorts moves in descending order up to and including
  // a given limit. The order of moves smaller than the limit is left unspecified.
  void partial_insertion_sort(ExtMove* begin, ExtMove* end, int limit) {

    for (ExtMove *sortedEnd = begin, *p = begin + 1; p < end; ++p)
        if (p->value >= limit)
        {
            ExtMove tmp = *p, *q;
            *p = *++sortedEnd;
            for (q = sortedEnd; q != begin && *(q - 1) < tmp; --q)
                *q = *(q - 1);
            *q = tmp;
        }
  }

} // namespace


/// Constructors of the MovePicker class. As arguments we pass information
/// to help it to return the (presumably) good moves first, to decide which
/// moves to return (in the quiescence search, for instance, we only want to
/// search captures, promotions, and some checks) and how important good move
/// ordering is at the current node.

/// MovePicker constructor for the main search
MovePicker::MovePicker(const Position& p, Move ttm, Depth d, const ButterflyHistory* mh, const GateHistory* dh, const LowPlyHistory* lp,
                       const CapturePieceToHistory* cph, const PieceToHistory** ch, Move cm, const Move* killers, int pl)
           : pos(p), mainHistory(mh), gateHistory(dh), lowPlyHistory(lp), captureHistory(cph), continuationHistory(ch),
             ttMove(ttm), refutations{{killers[0], 0}, {killers[1], 0}, {cm, 0}}, depth(d), ply(pl) {

  assert(d > 0);

  stage = (pos.checkers() ? EVASION_TT : MAIN_TT) +
          !(ttm && pos.pseudo_legal(ttm));
}

/// MovePicker constructor for quiescence search
MovePicker::MovePicker(const Position& p, Move ttm, Depth d, const ButterflyHistory* mh, const GateHistory* dh,
                       const CapturePieceToHistory* cph, const PieceToHistory** ch, Square rs)
           : pos(p), mainHistory(mh), gateHistory(dh), captureHistory(cph), continuationHistory(ch), ttMove(ttm), recaptureSquare(rs), depth(d) {

  assert(d <= 0);

  stage = (pos.checkers() ? EVASION_TT : QSEARCH_TT) +
          !(   ttm
            && (pos.checkers() || depth > DEPTH_QS_RECAPTURES || to_sq(ttm) == recaptureSquare)
            && pos.pseudo_legal(ttm));
}

/// MovePicker constructor for ProbCut: we generate captures with SEE greater
/// than or equal to the given threshold.
MovePicker::MovePicker(const Position& p, Move ttm, Value th, const GateHistory* dh, const CapturePieceToHistory* cph)
           : pos(p), gateHistory(dh), captureHistory(cph), ttMove(ttm), threshold(th) {

  assert(!pos.checkers());

  stage = PROBCUT_TT + !(ttm && pos.capture(ttm)
                             && pos.pseudo_legal(ttm)
                             && pos.see_ge(ttm, threshold));
}

/// MovePicker::score() assigns a numerical value to each move in a list, used
/// for sorting. Captures are ordered by Most Valuable Victim (MVV), preferring
/// captures with a good history. Quiets moves are ordered using the histories.
template<GenType Type>
void MovePicker::score() {

  static_assert(Type == CAPTURES || Type == QUIETS || Type == EVASIONS, "Wrong type");

  const bool battleKings = is_battle_kings(pos);

  for (auto& m : *this)
      if constexpr (Type == CAPTURES)
      {
          m.value =  int(PieceValue[MG][pos.piece_on(to_sq(m))]) * 6
                   + (*gateHistory)[pos.side_to_move()][gating_square(m)]
                   + (*captureHistory)[pos.moved_piece(m)][to_sq(m)][type_of(pos.piece_on(to_sq(m)))];

          if (battleKings)
              m.value += battle_kings_adjustment(pos, m);
      }

      else if constexpr (Type == QUIETS)
      {
          m.value =      (*mainHistory)[pos.side_to_move()][from_to(m)]
                   +     (*gateHistory)[pos.side_to_move()][gating_square(m)]
                   + 2 * (*continuationHistory[0])[history_slot(pos.moved_piece(m))][to_sq(m)]
                   +     (*continuationHistory[1])[history_slot(pos.moved_piece(m))][to_sq(m)]
                   +     (*continuationHistory[3])[history_slot(pos.moved_piece(m))][to_sq(m)]
                   +     (*continuationHistory[5])[history_slot(pos.moved_piece(m))][to_sq(m)]
                   + (ply < MAX_LPH ? std::min(4, depth / 3) * (*lowPlyHistory)[ply][from_to(m)] : 0);

          if (battleKings)
              m.value += battle_kings_adjustment(pos, m);
      }

      else // Type == EVASIONS
      {
          if (pos.capture(m))
              m.value =  PieceValue[MG][pos.piece_on(to_sq(m))]
                       - Value(type_of(pos.moved_piece(m)));
          else
              m.value =      (*mainHistory)[pos.side_to_move()][from_to(m)]
                       + 2 * (*continuationHistory[0])[history_slot(pos.moved_piece(m))][to_sq(m)]
                       - (1 << 28);
      }
}

/// MovePicker::select() returns the next move satisfying a predicate function.
/// It never returns the TT move.
template<MovePicker::PickType T, typename Pred>
Move MovePicker::select(Pred filter) {

  while (cur < endMoves)
  {
      if (T == Best)
          std::swap(*cur, *std::max_element(cur, endMoves));

      if (*cur != ttMove && filter())
          return *cur++;

      cur++;
  }
  return MOVE_NONE;
}

/// MovePicker::next_move() is the most important method of the MovePicker class. It
/// returns a new pseudo-legal move every time it is called until there are no more
/// moves left, picking the move with the highest score from a list of generated moves.
Move MovePicker::next_move(bool skipQuiets) {

top:
  switch (stage) {

  case MAIN_TT:
  case EVASION_TT:
  case QSEARCH_TT:
  case PROBCUT_TT:
      ++stage;
      assert(pos.legal(ttMove) == MoveList<LEGAL>(pos).contains(ttMove) || pos.virtual_drop(ttMove));
      return ttMove;

  case CAPTURE_INIT:
  case PROBCUT_INIT:
  case QCAPTURE_INIT:
      cur = endBadCaptures = moves;
      endMoves = generate<CAPTURES>(pos, cur);

      score<CAPTURES>();
      ++stage;
      goto top;

  case GOOD_CAPTURE:
      if (select<Best>([&](){
                       return pos.see_ge(*cur, Value(-69 * cur->value / 1024 - 500 * (pos.captures_to_hand() && pos.gives_check(*cur))))?
                              // Move losing capture to endBadCaptures to be tried later
                              true : (*endBadCaptures++ = *cur, false); }))
          return *(cur - 1);

      // Prepare the pointers to loop over the refutations array
      cur = std::begin(refutations);
      endMoves = std::end(refutations);

      // If the countermove is the same as a killer, skip it
      if (   refutations[0].move == refutations[2].move
          || refutations[1].move == refutations[2].move)
          --endMoves;

      ++stage;
      [[fallthrough]];

  case REFUTATION:
      if (select<Next>([&](){ return    *cur != MOVE_NONE
                                    && !pos.capture(*cur)
                                    &&  pos.pseudo_legal(*cur); }))
          return *(cur - 1);
      ++stage;
      [[fallthrough]];

  case QUIET_INIT:
      if (!skipQuiets && !(pos.must_capture() && pos.has_capture()))
      {
          cur = endBadCaptures;
          endMoves = generate<QUIETS>(pos, cur);

          score<QUIETS>();
          partial_insertion_sort(cur, endMoves, -3000 * depth);
      }

      ++stage;
      [[fallthrough]];

  case QUIET:
      if (   !skipQuiets
          && select<Next>([&](){return   *cur != refutations[0].move
                                      && *cur != refutations[1].move
                                      && *cur != refutations[2].move;}))
          return *(cur - 1);

      // Prepare the pointers to loop over the bad captures
      cur = moves;
      endMoves = endBadCaptures;

      ++stage;
      [[fallthrough]];

  case BAD_CAPTURE:
      return select<Next>([](){ return true; });

  case EVASION_INIT:
      cur = moves;
      endMoves = generate<EVASIONS>(pos, cur);

      score<EVASIONS>();
      ++stage;
      [[fallthrough]];

  case EVASION:
      return select<Best>([](){ return true; });

  case PROBCUT:
      return select<Best>([&](){ return pos.see_ge(*cur, threshold); });

  case QCAPTURE:
      if (select<Best>([&](){ return   depth > DEPTH_QS_RECAPTURES
                                    || to_sq(*cur) == recaptureSquare; }))
          return *(cur - 1);

      // If we did not find any move and we do not try checks, we have finished
      if (depth != DEPTH_QS_CHECKS)
          return MOVE_NONE;

      ++stage;
      [[fallthrough]];

  case QCHECK_INIT:
      cur = moves;
      endMoves = generate<QUIET_CHECKS>(pos, cur);

      ++stage;
      [[fallthrough]];

  case QCHECK:
      return select<Next>([](){ return true; });
  }

  assert(false);
  return MOVE_NONE; // Silence warning
}

} // namespace Stockfish
