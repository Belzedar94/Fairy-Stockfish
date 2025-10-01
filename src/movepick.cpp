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

  struct BattleKingsHeuristic {

    static constexpr std::array<PieceType, 6> CaptureOrder = {PAWN, KNIGHT, BISHOP, ROOK, QUEEN, COMMONER};

    bool enabled = false;
    Color us = WHITE;
    Color them = BLACK;
    std::array<int, CaptureOrder.size()> enemyRemaining{};
    int friendlyYoung = 0;
    int enemyYoung = 0;
    Rank maxRank = RANK_8;

    explicit BattleKingsHeuristic(const Position& pos) {
      enabled = is_battle_kings(pos);
      if (!enabled)
          return;

      us = pos.side_to_move();
      them = ~us;

      for (size_t i = 0; i < CaptureOrder.size(); ++i)
          enemyRemaining[i] = pos.count(them, CaptureOrder[i]);

      friendlyYoung = pos.count(us, PAWN) + pos.count(us, KNIGHT) + pos.count(us, BISHOP);
      enemyYoung = pos.count(them, PAWN) + pos.count(them, KNIGHT) + pos.count(them, BISHOP);
      maxRank = pos.max_rank();
    }

    int category_index(PieceType pt) const {
      for (size_t i = 0; i < CaptureOrder.size(); ++i)
          if (CaptureOrder[i] == pt)
              return int(i);
      return -1;
    }

    int next_target_index() const {
      for (size_t i = 0; i < CaptureOrder.size(); ++i)
          if (enemyRemaining[i])
              return int(i);
      return -1;
    }

    int young_remaining() const {
      return enemyRemaining[0] + enemyRemaining[1] + enemyRemaining[2];
    }

    int capture_score(PieceType victim) const {
      if (victim == NO_PIECE_TYPE)
          return 0;

      int idx = category_index(victim);
      if (idx == -1)
          return 0;

      if (victim == COMMONER)
          return 6000;

      int target = next_target_index();
      if (target == -1)
          return 0;

      int score = 700 + (int(CaptureOrder.size()) - 1 - idx) * 140;
      if (idx == target)
          score += 450;
      else if (idx > target)
          score -= 300 * (idx - target);
      if (enemyRemaining[idx] == 1)
          score += 180;
      return score;
    }

    int mover_score(PieceType mover, bool isCapture) const {
      int score = 0;

      switch (mover)
      {
      case PAWN:
          score += 460;
          if (!isCapture)
              score += 120;
          break;

      case KNIGHT:
          score += 320;
          break;

      case BISHOP:
          score += 200;
          break;

      case ROOK:
          score -= 220;
          if (isCapture)
              score += 140;
          break;

      case QUEEN:
          score -= 900;
          if (!isCapture)
              score -= 450;
          break;

      case COMMONER:
          score -= 480;
          break;

      default:
          break;
      }

      if (friendlyYoung <= enemyYoung && (mover == PAWN || mover == KNIGHT || mover == BISHOP))
          score += 180;

      if (friendlyYoung > enemyYoung && (mover == ROOK || mover == QUEEN))
          score -= 120;

      if (young_remaining())
      {
          if (mover == PAWN || mover == KNIGHT || mover == BISHOP)
              score += 120;
          else if (mover == ROOK || mover == QUEEN)
              score -= 140;
      }

      return score;
    }

    int gating_score(const Position& pos, Move m) const {
      PieceType gate = gating_type(m);
      if (gate == NO_PIECE_TYPE)
          return 0;

      int score = 0;
      switch (gate)
      {
      case KNIGHT:
          score += 1300;
          if (friendlyYoung <= enemyYoung)
              score += 160;
          if (young_remaining())
              score += 120;
          break;

      case BISHOP:
          score += 950;
          if (friendlyYoung <= enemyYoung)
              score += 120;
          if (young_remaining())
              score += 80;
          break;

      case ROOK:
          score -= 320;
          if (young_remaining())
              score -= 220;
          break;

      case QUEEN:
          score -= 1500;
          if (young_remaining())
              score -= 240;
          break;

      case COMMONER:
          score -= 3800;
          break;

      default:
          break;
      }

      Square sq = gating_square(m);
      if (sq != SQ_NONE)
      {
          int attackers = popcount(pos.attackers_to(sq, them));
          int defenders = popcount(pos.attackers_to(sq, us));

          if (gate == KNIGHT || gate == BISHOP)
          {
              score += defenders * 60;
              if (!attackers)
                  score += 160;
          }

          if (attackers > defenders)
              score -= 200 * (attackers - defenders);
      }

      return score;
    }

    int queen_penalty(const Position& pos, Move m, PieceType mover, bool isCapture, PieceType victim) const {
      if (mover != QUEEN)
          return 0;

      int score = -600;

      if (!isCapture)
          score -= 500;

      if (young_remaining())
          score -= 500;

      if (gating_type(m) == COMMONER)
          score -= 1600;

      if (victim == COMMONER)
          score += 6000;
      else if (isCapture)
          score += capture_score(victim) / 2;

      Square from = from_sq(m);
      if (from != SQ_NONE && pos.attackers_to(from, them))
          score -= 220;

      return score;
    }

    int score(const Position& pos, Move m) const {
      if (!enabled)
          return 0;

      Piece moved = pos.moved_piece(m);
      PieceType mover = type_of(moved);
      bool isCapture = pos.capture(m);
      Piece captured = isCapture ? pos.piece_on(to_sq(m)) : NO_PIECE;
      PieceType victim = captured != NO_PIECE ? type_of(captured) : NO_PIECE_TYPE;

      int total = mover_score(mover, isCapture);
      total += gating_score(pos, m);

      if (victim != NO_PIECE_TYPE)
          total += capture_score(victim);

      total += queen_penalty(pos, m, mover, isCapture, victim);

      Square from = from_sq(m);
      Square to = to_sq(m);
      if (from != SQ_NONE && to != SQ_NONE)
      {
          int fromRank = int(relative_rank(us, from, maxRank));
          int toRank = int(relative_rank(us, to, maxRank));
          int delta = toRank - fromRank;

          if (mover == PAWN)
              total += delta * 140;
          else if (mover == KNIGHT || mover == BISHOP)
              total += delta * 80;
          else if (mover == ROOK || mover == QUEEN)
              total += delta * 40;
      }

      return total;
    }
  };

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

  const BattleKingsHeuristic battle(pos);

  for (auto& m : *this)
      if constexpr (Type == CAPTURES)
      {
          m.value =  int(PieceValue[MG][pos.piece_on(to_sq(m))]) * 6
                   + (*gateHistory)[pos.side_to_move()][gating_square(m)]
                   + (*captureHistory)[pos.moved_piece(m)][to_sq(m)][type_of(pos.piece_on(to_sq(m)))];

          if (battle.enabled)
              m.value += battle.score(pos, m);
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

          if (battle.enabled)
              m.value += battle.score(pos, m);
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
