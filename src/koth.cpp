/*
  KOTH-Stockfish, a King of the Hill specialization of Stockfish
  Copyright (C) 2026 The KOTH-Stockfish developers

  KOTH-Stockfish is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
*/

#include "koth.h"

#include <sstream>
#include <utility>

#include "movegen.h"
#include "position.h"

namespace Stockfish::Koth {

namespace {

void add_predicate(GameStatus& status, TerminalPredicate predicate) {
    status.predicates |= static_cast<std::uint16_t>(predicate);
}

struct PositionSnapshot {
    std::string       fen;
    StateInfo         state;
    std::vector<Move> moves;
};

bool same_transition(const TransitionEvent& left, const TransitionEvent& right) {
    return left.move == right.move && left.kind == right.kind
        && left.kingEnteredHill == right.kingEnteredHill;
}

bool same_state(const StateInfo& left, const StateInfo& right) {
    return left.materialKey == right.materialKey && left.pawnKey == right.pawnKey
        && left.minorPieceKey == right.minorPieceKey
        && left.nonPawnKey[WHITE] == right.nonPawnKey[WHITE]
        && left.nonPawnKey[BLACK] == right.nonPawnKey[BLACK]
        && left.nonPawnMaterial[WHITE] == right.nonPawnMaterial[WHITE]
        && left.nonPawnMaterial[BLACK] == right.nonPawnMaterial[BLACK]
        && left.castlingRights == right.castlingRights && left.rule50 == right.rule50
        && left.pliesFromNull == right.pliesFromNull && left.epSquare == right.epSquare
        && left.fenEpSquare == right.fenEpSquare && left.key == right.key
        && left.checkersBB == right.checkersBB && left.previous == right.previous
        && left.blockersForKing[WHITE] == right.blockersForKing[WHITE]
        && left.blockersForKing[BLACK] == right.blockersForKing[BLACK]
        && left.pinners[WHITE] == right.pinners[WHITE]
        && left.pinners[BLACK] == right.pinners[BLACK]
        && left.checkSquares[PAWN] == right.checkSquares[PAWN]
        && left.checkSquares[KNIGHT] == right.checkSquares[KNIGHT]
        && left.checkSquares[BISHOP] == right.checkSquares[BISHOP]
        && left.checkSquares[ROOK] == right.checkSquares[ROOK]
        && left.checkSquares[QUEEN] == right.checkSquares[QUEEN]
        && left.checkSquares[KING] == right.checkSquares[KING]
        && left.capturedPiece == right.capturedPiece && left.repetition == right.repetition
        && left.repetitionCount == right.repetitionCount
        && same_transition(left.transition, right.transition);
}

PositionSnapshot snapshot(const Position& position) {
    return {position.fen(), *position.state(), physical_moves(position)};
}

bool same_snapshot(const PositionSnapshot& left, const PositionSnapshot& right) {
    return left.fen == right.fen && same_state(left.state, right.state)
        && left.moves == right.moves;
}

std::optional<std::string> test_position_round_trip(std::string_view fen, int& checks) {
    StateInfo rootState;
    Position  position;
    if (auto error = position.set(std::string(fen), false, &rootState))
        return "root parse failed: " + std::string(error->what());

    const auto root  = snapshot(position);
    const auto moves = physical_moves(position);

    for (Move move : moves)
    {
        const Piece  moved        = position.moved_piece(move);
        const Square from         = move.from_sq();
        const Square encodedTo    = move.to_sq();
        const bool   expectedHill = type_of(moved) == KING && move.type_of() != CASTLING
                               && is_hill(encodedTo) && !is_hill(from);

        StateInfo childState;
        position.do_move(move, childState);
        ++checks;

        if (position.transition().kind != TransitionKind::ACCEPTED_MOVE
            || position.transition().move != move
            || position.transition().kingEnteredHill != expectedHill)
            return "real-move transition mismatch move=" + std::to_string(move.raw());

        const auto child = snapshot(position);
        for (Move reply : physical_moves(position))
        {
            StateInfo grandchildState;
            position.do_move(reply, grandchildState);
            position.undo_move(reply);
            ++checks;
            if (!same_snapshot(snapshot(position), child))
                return "nested make/undo mismatch move=" + std::to_string(move.raw())
                     + " reply=" + std::to_string(reply.raw());
        }

        position.undo_move(move);
        ++checks;
        if (!same_snapshot(snapshot(position), root))
            return "make/undo mismatch move=" + std::to_string(move.raw());
    }

    if (!position.checkers())
    {
        StateInfo nullState;
        position.do_null_move(nullState);
        ++checks;
        if (position.transition().kind != TransitionKind::SEARCH_NULL
            || position.transition().move != Move::none() || position.transition().kingEnteredHill
            || position.repetition_count() != 0 || position.ep_square() != SQ_NONE
            || position.fen_ep_square() != SQ_NONE)
            return "null-move state was not cleared";
        position.undo_null_move();
        ++checks;
        if (!same_snapshot(snapshot(position), root))
            return "null make/undo mismatch";
    }

    return std::nullopt;
}

}  // namespace

bool any_king_on_hill(const Position& position) {
    return bool(position.pieces(KING) & HillSquares);
}

bool is_canonical_uci_move(std::string_view token) {
    if (token.size() != 4 && token.size() != 5)
        return false;

    const auto file = [](char c) { return c >= 'a' && c <= 'h'; };
    const auto rank = [](char c) { return c >= '1' && c <= '8'; };

    if (!file(token[0]) || !rank(token[1]) || !file(token[2]) || !rank(token[3]))
        return false;

    return token.size() == 4 || token[4] == 'q' || token[4] == 'r' || token[4] == 'b'
        || token[4] == 'n';
}

GameStatus classify(const Position& position, AdjudicationContext context) {
    GameStatus status;

    if (context == AdjudicationContext::SEARCH_NULL)
        return status;

    const bool hasLegalMove = MoveList<LEGAL>(position).size() != 0;
    const bool inCheck      = bool(position.checkers());

    if (!hasLegalMove && inCheck)
        add_predicate(status, TerminalPredicate::CHECKMATE);

    const auto& transition = position.transition();
    if (context != AdjudicationContext::RAW_ROOT_VALIDATION
        && transition.kind == TransitionKind::ACCEPTED_MOVE && transition.kingEnteredHill)
        add_predicate(status, TerminalPredicate::HILL);

    if (!hasLegalMove && !inCheck)
        add_predicate(status, TerminalPredicate::STALEMATE);

    if (position.repetition_count() >= 3)
        add_predicate(status, TerminalPredicate::REPETITION3_AUTO);

    if (position.repetition_count() >= 5)
        add_predicate(status, TerminalPredicate::REPETITION5_DIAGNOSTIC);

    if (position.rule50_count() >= 100)
        add_predicate(status, TerminalPredicate::RULE50_AUTO);

    if (status.has(TerminalPredicate::CHECKMATE))
        status.primary = PrimaryReason::CHECKMATE;
    else if (status.has(TerminalPredicate::HILL))
        status.primary = PrimaryReason::HILL;
    else if (status.has(TerminalPredicate::STALEMATE))
        status.primary = PrimaryReason::STALEMATE;
    else if (status.has(TerminalPredicate::REPETITION3_AUTO)
             || status.has(TerminalPredicate::RULE50_AUTO))
        status.primary = PrimaryReason::AUTOMATIC_DRAW;

    if (status.primary == PrimaryReason::CHECKMATE || status.primary == PrimaryReason::HILL)
        status.winner = ~position.side_to_move();

    return status;
}

std::vector<Move> physical_moves(const Position& position) {
    std::vector<Move> moves;
    for (Move move : MoveList<LEGAL>(position))
        moves.push_back(move);
    return moves;
}

std::vector<Move> game_moves(const Position& position) {
    if (classify(position).terminal())
        return {};
    return physical_moves(position);
}

bool is_immediate_hill_move(const Position& position, Move move) {
    return move.is_ok() && move.type_of() != CASTLING && type_of(position.moved_piece(move)) == KING
        && !is_hill(move.from_sq()) && is_hill(move.to_sq());
}

std::vector<Move> immediate_hill_moves(const Position& position) {
    std::vector<Move> moves;

    const auto context = position.transition().kind == TransitionKind::SEARCH_NULL
                         ? AdjudicationContext::SEARCH_NULL
                         : AdjudicationContext::ACCEPTED_TRAJECTORY;
    if (classify(position, context).terminal())
        return moves;

    for (Move move : MoveList<LEGAL>(position))
        if (is_immediate_hill_move(position, move))
            moves.push_back(move);

    return moves;
}

std::string primary_reason_name(PrimaryReason reason) {
    switch (reason)
    {
    case PrimaryReason::NONE :
        return "NONE";
    case PrimaryReason::CHECKMATE :
        return "CHECKMATE";
    case PrimaryReason::HILL :
        return "HILL";
    case PrimaryReason::STALEMATE :
        return "STALEMATE";
    case PrimaryReason::AUTOMATIC_DRAW :
        return "AUTOMATIC_DRAW";
    }
    return "UNKNOWN";
}

std::string serialize(const GameStatus& status) {
    std::ostringstream out;
    out << "profile=KOTH_LICHESS_V1 terminal=" << (status.terminal() ? "true" : "false")
        << " predicates=" << status.predicates << " primary=" << primary_reason_name(status.primary)
        << " winner=";

    if (!status.winner)
        out << "none";
    else
        out << (*status.winner == WHITE ? "white" : "black");

    return out.str();
}

std::string state_selftest() {
    static constexpr std::string_view TestFens[] = {
      "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1",
      "r3k2r/8/8/8/8/8/8/R3K2R w KQkq - 0 1",
      "7k/P7/8/8/8/8/8/7K w - - 0 1",
      "7k/8/8/8/3pP3/8/8/K7 b - e3 0 1",
      "7k/8/8/8/8/2K5/8/8 w - - 99 1",
      "7k/8/8/8/2R5/8/8/K7 w - - 0 1",
    };

    int checks = 0;
    for (std::string_view fen : TestFens)
        if (auto failure = test_position_round_trip(fen, checks))
            return "FAIL checks=" + std::to_string(checks) + " detail=" + *failure;

    return "PASS cases=" + std::to_string(std::size(TestFens))
         + " checks=" + std::to_string(checks);
}

}  // namespace Stockfish::Koth
