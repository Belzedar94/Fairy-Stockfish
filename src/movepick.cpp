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

#include <cassert>

#include "bitboard.h"
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

  constexpr PieceType BattleKingsOrder[] = {PAWN, KNIGHT, BISHOP, ROOK, QUEEN, COMMONER};
  constexpr int BattleKingsOrderCount = int(sizeof(BattleKingsOrder) / sizeof(BattleKingsOrder[0]));

  constexpr int StageBaseMover[BattleKingsOrderCount]   = {900, 620, 280, -420, -1300, -1800};
  constexpr int StageBaseGate[BattleKingsOrderCount]    = {740, 540, 260, -560, -1700, -3400};
  constexpr int StageScarcityScale[BattleKingsOrderCount] = {140, 160, 190, 250, 380, 560};
  constexpr int StageAgingScale[BattleKingsOrderCount]    = {70, 90, 130, 200, 300, 420};
  constexpr int StageSafetyScale[BattleKingsOrderCount]   = {160, 190, 220, 320, 460, 680};
  constexpr int StageCaptureBase[BattleKingsOrderCount]   = {1800, 1600, 1300, 1100, 900, 700};
  constexpr int StageCaptureSwing[BattleKingsOrderCount]  = {140, 180, 220, 300, 380, 480};
  constexpr int StageCentralWeight[BattleKingsOrderCount] = {140, 170, 180, 110, 70, 40};
  constexpr int StageForwardWeight[BattleKingsOrderCount] = {220, 180, 130, 50, 20, 0};

  int battle_kings_stage_index(PieceType pt) {
    for (int i = 0; i < BattleKingsOrderCount; ++i)
        if (BattleKingsOrder[i] == pt)
            return i;

    return BattleKingsOrderCount;
  }

  struct BattleKingsContext {
    Color us;
    Color them;
    int boardArea;
    int occupancy;
    int ourCount[BattleKingsOrderCount];
    int theirCount[BattleKingsOrderCount];

    explicit BattleKingsContext(const Position& pos)
      : us(pos.side_to_move()),
        them(~pos.side_to_move()),
        boardArea((int(pos.max_file()) + 1) * (int(pos.max_rank()) + 1)),
        occupancy(popcount(pos.pieces())) {

      for (int i = 0; i < BattleKingsOrderCount; ++i)
      {
          ourCount[i]   = pos.count(us,   BattleKingsOrder[i]);
          theirCount[i] = pos.count(them, BattleKingsOrder[i]);
      }
    }
  };

  int battle_kings_mover_flow(const BattleKingsContext& ctx, PieceType mover) {
    int stage = battle_kings_stage_index(mover);
    if (stage >= BattleKingsOrderCount)
        return 0;

    int bonus = StageBaseMover[stage];
    int imbalance = ctx.ourCount[stage] - ctx.theirCount[stage];

    bonus -= imbalance * StageScarcityScale[stage];

    int oldLead = 0;
    for (int i = stage + 1; i < BattleKingsOrderCount; ++i)
        oldLead += std::max(0, ctx.ourCount[i] - ctx.theirCount[i]);

    if (stage <= 2)
        bonus += oldLead * StageAgingScale[stage];
    else
        bonus -= oldLead * StageAgingScale[stage];

    return bonus;
  }

  int battle_kings_spawn_preference(const BattleKingsContext& ctx, PieceType gate) {
    int stage = battle_kings_stage_index(gate);
    if (stage >= BattleKingsOrderCount)
        return 0;

    int bonus = StageBaseGate[stage];
    int scarcity = ctx.theirCount[stage] - ctx.ourCount[stage];

    bonus += scarcity * (StageScarcityScale[stage] + StageScarcityScale[stage] / 2);

    int excessAged = 0;
    for (int i = stage + 1; i < BattleKingsOrderCount; ++i)
        excessAged += std::max(0, ctx.ourCount[i] - ctx.theirCount[i]);

    bonus -= excessAged * StageAgingScale[stage];

    if (stage <= 2)
        bonus += (ctx.boardArea - ctx.occupancy) * 8;
    else
        bonus -= ctx.occupancy * (stage - 2) * 6;

    if (ctx.ourCount[stage] > ctx.theirCount[stage] + 1)
        bonus -= (ctx.ourCount[stage] - ctx.theirCount[stage] - 1) * StageScarcityScale[stage] * (stage + 1);

    return bonus;
  }

  int battle_kings_spawn_safety(const Position& pos, const BattleKingsContext& ctx, Move m, PieceType gate, Bitboard occAfter) {
    Square gateSq = gating_square(m);
    if (gateSq == SQ_NONE)
        gateSq = from_sq(m);

    if (gateSq == SQ_NONE)
        return 0;

    int stage = battle_kings_stage_index(gate);
    if (stage >= BattleKingsOrderCount)
        return 0;

    Bitboard enemyAttackers = pos.attackers_to(gateSq, occAfter, ctx.them);
    Bitboard ourAttackers   = pos.attackers_to(gateSq, occAfter, ctx.us);

    int enemyCount = popcount(enemyAttackers);
    int ourCount   = popcount(ourAttackers);

    int deficit = enemyCount - ourCount;
    int scale = StageSafetyScale[stage];

    if (deficit <= 0)
        return (ourCount + 1) * (scale / 2 + (stage <= 2 ? 70 : 30));

    int territory = int(relative_rank(ctx.us, gateSq, pos.max_rank()));
    return -(deficit * scale * (stage + 1) + territory * (stage + 1) * 20);
  }

  int battle_kings_flow_bonus(const Position& pos, const BattleKingsContext& ctx, Move m, PieceType mover, bool isCapture) {
    Square from = from_sq(m);
    Square to = to_sq(m);
    if (from == SQ_NONE || to == SQ_NONE)
        return 0;

    int stage = battle_kings_stage_index(mover);
    if (stage >= BattleKingsOrderCount)
        stage = BattleKingsOrderCount - 1;

    int forwardGain = int(relative_rank(ctx.us, to, pos.max_rank())) - int(relative_rank(ctx.us, from, pos.max_rank()));
    int centralGain =   edge_distance(file_of(to), pos.max_file())
                      + edge_distance(rank_of(to), pos.max_rank())
                      - edge_distance(file_of(from), pos.max_file())
                      - edge_distance(rank_of(from), pos.max_rank());

    int bonus = forwardGain * StageForwardWeight[stage] + centralGain * StageCentralWeight[stage];

    if (ctx.occupancy > ctx.boardArea * 3 / 4 && stage <= 2)
        bonus += (ctx.occupancy - ctx.boardArea * 3 / 4) * 6;

    if (forwardGain < 0 && stage <= 2)
        bonus += forwardGain * StageForwardWeight[stage] / 2;

    if (isCapture)
        bonus /= 2;

    return bonus;
  }

  int battle_kings_capture_swing(const BattleKingsContext& ctx, PieceType mover, PieceType victim) {
    int victimStage = battle_kings_stage_index(victim);
    if (victimStage >= BattleKingsOrderCount)
        return 0;

    int moverStage = battle_kings_stage_index(mover);
    if (moverStage >= BattleKingsOrderCount)
        moverStage = BattleKingsOrderCount - 1;

    int bonus = StageCaptureBase[victimStage];

    int scarcity = ctx.theirCount[victimStage] - ctx.ourCount[victimStage];
    bonus += scarcity * StageCaptureSwing[victimStage];

    if (victim == COMMONER)
        bonus += 3000 + 400 * ctx.theirCount[victimStage];

    bonus += (victimStage - moverStage) * StageCaptureSwing[victimStage];

    if (victimStage >= 3 && ctx.occupancy > ctx.boardArea / 2)
        bonus += (ctx.occupancy - ctx.boardArea / 2) * 6;

    return bonus;
  }

  int battle_kings_adjustment(const Position& pos, Move m) {
    BattleKingsContext ctx(pos);

    Square from = from_sq(m);
    Square to = to_sq(m);
    PieceType mover = type_of(pos.moved_piece(m));
    PieceType gate = gating_type(m);
    Square gateSq = gating_square(m);
    if (gateSq == SQ_NONE)
        gateSq = from;

    Square captureSq = SQ_NONE;
    Piece captured = NO_PIECE;
    if (pos.capture(m))
    {
        captureSq = type_of(m) == EN_PASSANT ? pos.capture_square(to) : to;
        captured = pos.piece_on(captureSq);
    }

    Bitboard occAfter = pos.pieces();
    if (is_ok(from))
        occAfter &= ~square_bb(from);
    if (is_ok(captureSq))
        occAfter &= ~square_bb(captureSq);
    if (is_ok(to))
        occAfter |= square_bb(to);
    if (gate != NO_PIECE_TYPE && is_ok(gateSq))
        occAfter |= square_bb(gateSq);

    int bonus = battle_kings_mover_flow(ctx, mover);

    if (gate != NO_PIECE_TYPE)
    {
        bonus += battle_kings_spawn_preference(ctx, gate);
        bonus += battle_kings_spawn_safety(pos, ctx, m, gate, occAfter);
    }
    else if (mover == COMMONER || mover == KING)
    {
        int idx = battle_kings_stage_index(COMMONER);
        if (idx < BattleKingsOrderCount)
            bonus -= (ctx.ourCount[idx] - ctx.theirCount[idx]) * 400;
    }

    if (captured != NO_PIECE)
        bonus += battle_kings_capture_swing(ctx, mover, type_of(captured));

    bonus += battle_kings_flow_bonus(pos, ctx, m, mover, captured != NO_PIECE);

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
