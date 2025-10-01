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

  constexpr PieceType BattleKingsOrder[] = {PAWN, KNIGHT, BISHOP, ROOK, QUEEN, COMMONER};
  constexpr int BattleKingsOrderCount = int(sizeof(BattleKingsOrder) / sizeof(BattleKingsOrder[0]));

  struct BattleKingsStats {
    PieceType target = NO_PIECE_TYPE;
    PieceType shortage = NO_PIECE_TYPE;
    int pressureDiff = 0;
    int kingBalance = 0;
    File maxFile = FILE_H;
    Rank maxRank = RANK_8;
  };

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

  constexpr int battle_kings_stage_index(PieceType pt) {
    for (int i = 0; i < BattleKingsOrderCount; ++i)
        if (BattleKingsOrder[i] == pt)
            return i;

    return BattleKingsOrderCount;
  }

  PieceType battle_kings_target_piece(const Position& pos, Color attacker) {
    Color defender = ~attacker;

    for (PieceType pt : BattleKingsOrder)
        if (pos.count(defender, pt))
            return pt;

    return NO_PIECE_TYPE;
  }

  PieceType battle_kings_shortage_piece(const Position& pos, Color us) {
    Color them = ~us;

    for (PieceType pt : BattleKingsOrder)
        if (pos.count(us, pt) < pos.count(them, pt))
            return pt;

    return NO_PIECE_TYPE;
  }

  BattleKingsStats battle_kings_collect(const Position& pos) {
    Color us = pos.side_to_move();
    BattleKingsStats stats;

    stats.target = battle_kings_target_piece(pos, us);
    stats.shortage = battle_kings_shortage_piece(pos, us);

    for (int i = 0; i < BattleKingsOrderCount; ++i)
    {
        PieceType pt = BattleKingsOrder[i];
        int weight = (i + 1) * 6;
        stats.pressureDiff += weight * (pos.count(us, pt) - pos.count(~us, pt));
    }

    stats.kingBalance = pos.count(us, COMMONER) - pos.count(~us, COMMONER);
    stats.maxFile = pos.max_file();
    stats.maxRank = pos.max_rank();

    return stats;
  }

  int battle_kings_mover_bonus(PieceType pt) {
    switch (pt)
    {
    case PAWN:     return 900;
    case KNIGHT:   return 600;
    case BISHOP:   return 400;
    case ROOK:     return -600;
    case QUEEN:    return -2000;
    case COMMONER: return -800;
    default:       return 0;
    }
  }

  int battle_kings_gate_bonus(PieceType pt) {
    switch (pt)
    {
    case KNIGHT:   return 1000;
    case BISHOP:   return 700;
    case ROOK:     return -500;
    case QUEEN:    return -1500;
    case COMMONER: return -4000;
    default:       return 0;
    }
  }

  int battle_kings_gate_alignment_bonus(PieceType shortage, PieceType gate) {
    if (shortage == NO_PIECE_TYPE || gate == NO_PIECE_TYPE)
        return 0;

    int gateStage = battle_kings_stage_index(gate);
    int shortageStage = battle_kings_stage_index(shortage);

    if (gateStage == shortageStage)
        return 1600;

    if (gateStage == shortageStage + 1)
        return 600;

    if (gateStage > shortageStage + 1)
        return -500 * (gateStage - shortageStage);

    return 300 * (shortageStage - gateStage);
  }

  int battle_kings_gate_pressure_bonus(int pressureDiff, PieceType gate) {
    if (gate == NO_PIECE_TYPE)
        return 0;

    int gateStage = battle_kings_stage_index(gate);
    if (gateStage >= BattleKingsOrderCount)
        return 0;

    if (pressureDiff > 0 && gateStage >= 3)
        return -std::min(2400, pressureDiff * (gateStage - 2) * 40);

    if (pressureDiff < 0 && gateStage <= 2)
        return std::min(1800, (-pressureDiff) * (3 - gateStage) * 30);

    return 0;
  }

  int battle_kings_centrality(File f, Rank r, File maxFile, Rank maxRank) {
    return edge_distance(f, maxFile) + edge_distance(r, maxRank);
  }

  int battle_kings_gate_safety(const Position& pos, Move m, PieceType gate, const BattleKingsStats& stats) {
    if (gate == NO_PIECE_TYPE)
        return 0;

    Square gateSq = gating_square(m);
    if (gateSq == SQ_NONE)
        return 0;

    Color us = pos.side_to_move();
    Square from = from_sq(m);
    Square to = to_sq(m);

    Bitboard occ = pos.pieces();

    if (from != SQ_NONE)
        occ ^= from;

    occ |= to;

    if (pos.capture(m))
    {
        Square csq = type_of(m) == EN_PASSANT ? pos.capture_square(to) : to;
        occ -= csq;
    }

    occ |= gateSq;

    Bitboard attackers = pos.attackers_to(gateSq, occ, ~us);
    Bitboard defenders = pos.attackers_to(gateSq, occ, us);

    int stage = battle_kings_stage_index(gate);
    int centrality = battle_kings_centrality(file_of(gateSq), rank_of(gateSq), stats.maxFile, stats.maxRank);
    int base = 140 + stage * 90 + centrality * 20;

    if (!attackers)
        return (stage <= 1 ? 120 : 60) * (1 + centrality);

    if (defenders)
        return -base / 2;

    int penalty = base;
    if (stage >= BattleKingsOrderCount - 1)
        penalty *= 3;
    else if (stage >= BattleKingsOrderCount - 2)
        penalty *= 2;

    return -penalty;
  }

  int battle_kings_capture_bonus(PieceType pt) {
    switch (pt)
    {
    case COMMONER: return 9000;
    case PAWN:     return 5000;
    case KNIGHT:   return 3200;
    case BISHOP:   return 2200;
    case ROOK:     return 1200;
    case QUEEN:    return 600;
    default:       return 0;
    }
  }

  int battle_kings_adjustment(const Position& pos, Move m, const BattleKingsStats& stats) {
    int bonus = 0;

    PieceType mover = type_of(pos.moved_piece(m));
    bonus += battle_kings_mover_bonus(mover);

    PieceType gate = gating_type(m);
    if (gate != NO_PIECE_TYPE)
    {
        bonus += battle_kings_gate_bonus(gate);
        bonus += battle_kings_gate_alignment_bonus(stats.shortage, gate);
        bonus += battle_kings_gate_pressure_bonus(stats.pressureDiff, gate);
        bonus += battle_kings_gate_safety(pos, m, gate, stats);

        if (gate == COMMONER)
        {
            if (stats.kingBalance >= 0)
                bonus -= 3000 + 300 * stats.kingBalance;
            else
                bonus += 400 * (-stats.kingBalance);
        }
    }

    if (stats.shortage != NO_PIECE_TYPE && mover == stats.shortage)
        bonus += 500;

    Piece captured = NO_PIECE;
    PieceType victim = NO_PIECE_TYPE;

    if (pos.capture(m))
    {
        Square capSq = type_of(m) == EN_PASSANT ? pos.capture_square(to_sq(m)) : to_sq(m);
        captured = pos.piece_on(capSq);
        if (captured != NO_PIECE)
            victim = type_of(captured);
    }

    if (victim != NO_PIECE_TYPE)
    {
        int victimStage = battle_kings_stage_index(victim);

        bonus += battle_kings_capture_bonus(victim);

        if (stats.target == victim)
            bonus += 2800;
        else if (stats.target != NO_PIECE_TYPE && victimStage > battle_kings_stage_index(stats.target))
            bonus += 400;

        if (victimStage >= BattleKingsOrderCount - 2)
        {
            int needKings = std::max(0, -stats.kingBalance);
            bonus += 1800 + 250 * needKings;
        }

        if (stats.shortage == victim)
            bonus += 900;

        if (mover == QUEEN && victimStage <= battle_kings_stage_index(ROOK))
            bonus -= 2600;
    }

    if (stats.pressureDiff > 0)
    {
        int moverStage = battle_kings_stage_index(mover);
        if (moverStage >= 3)
        {
            Square from = from_sq(m);
            Square to = to_sq(m);
            Rank relFrom = from != SQ_NONE ? relative_rank(pos.side_to_move(), from, stats.maxRank) : RANK_1;
            Rank relTo = relative_rank(pos.side_to_move(), to, stats.maxRank);

            if (relTo <= relFrom)
                bonus += 300;
            else
                bonus -= 200;
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
  BattleKingsStats battleStats;

  if (battleKings)
      battleStats = battle_kings_collect(pos);

  for (auto& m : *this)
      if constexpr (Type == CAPTURES)
      {
          m.value =  int(PieceValue[MG][pos.piece_on(to_sq(m))]) * 6
                   + (*gateHistory)[pos.side_to_move()][gating_square(m)]
                   + (*captureHistory)[pos.moved_piece(m)][to_sq(m)][type_of(pos.piece_on(to_sq(m)))];

          if (battleKings)
              m.value += battle_kings_adjustment(pos, m, battleStats);
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
              m.value += battle_kings_adjustment(pos, m, battleStats);
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
