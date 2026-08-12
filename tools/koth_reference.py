#!/usr/bin/env python3
"""Independent executable reference for KOTH_LICHESS_V1.

This module deliberately depends on python-chess only for standard physical
chess state and legal moves. KOTH trajectory and terminal semantics live here,
not in python-chess's variant result helpers.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from dataclasses import asdict, dataclass
from typing import Any, Final, Iterable

import chess


INPUT_SCHEMA: Final = "koth-reference-input-v1"
OUTPUT_SCHEMA: Final = "koth-reference-output-v1"
PROFILE: Final = "KOTH_LICHESS_V1"
REFERENCE_ID: Final = "koth-python-reference-v1"
GOALS: Final = frozenset((chess.D4, chess.E4, chess.D5, chess.E5))
PREDICATE_ORDER: Final = (
    "CHECKMATE",
    "HILL",
    "STALEMATE",
    "THREEFOLD",
    "FIVEFOLD",
    "FIFTY_MOVE",
)


class ReferenceError(Exception):
    """Stable fail-closed input or trajectory error."""

    def __init__(self, code: str, detail: str) -> None:
        super().__init__(detail)
        self.code = code
        self.detail = detail


@dataclass(frozen=True)
class TerminalStatus:
    terminal: bool
    predicates: tuple[str, ...]
    primary_reason: str
    winner: str | None
    result: str


def canonical_fen(board: chess.Board) -> str:
    return board.fen(en_passant="fen")


def sorted_legal_uci(board: chess.Board) -> list[str]:
    return sorted(move.uci() for move in board.legal_moves)


def legal_moves_sha256(moves: Iterable[str]) -> str:
    payload = ("\n".join(moves) + "\n").encode("ascii")
    return hashlib.sha256(payload).hexdigest().upper()


def goal_occupants(board: chess.Board) -> tuple[str, ...]:
    occupants: list[str] = []
    for color, label in ((chess.WHITE, "white"), (chess.BLACK, "black")):
        square = board.king(color)
        if square in GOALS:
            occupants.append(label)
    return tuple(occupants)


def terminal_status(
    board: chess.Board,
    *,
    hill_transition: bool,
    mover: chess.Color | None,
) -> TerminalStatus:
    present: set[str] = set()
    if board.is_checkmate():
        present.add("CHECKMATE")
    if hill_transition:
        present.add("HILL")
    if board.is_stalemate():
        present.add("STALEMATE")
    if board.is_repetition(3):
        present.add("THREEFOLD")
    if board.is_repetition(5):
        present.add("FIVEFOLD")
    if board.halfmove_clock >= 100:
        present.add("FIFTY_MOVE")

    predicates = tuple(name for name in PREDICATE_ORDER if name in present)
    if "CHECKMATE" in present:
        primary = "CHECKMATE"
        winner_color = not board.turn
    elif "HILL" in present:
        primary = "HILL"
        winner_color = mover
    elif "STALEMATE" in present:
        primary = "STALEMATE"
        winner_color = None
    elif present.intersection(("THREEFOLD", "FIVEFOLD", "FIFTY_MOVE")):
        primary = "AUTOMATIC_DRAW"
        winner_color = None
    else:
        primary = "NONE"
        winner_color = None

    if winner_color is chess.WHITE:
        winner, result = "white", "1-0"
    elif winner_color is chess.BLACK:
        winner, result = "black", "0-1"
    else:
        winner = None
        result = "1/2-1/2" if primary in ("STALEMATE", "AUTOMATIC_DRAW") else "*"

    return TerminalStatus(
        terminal=primary != "NONE",
        predicates=predicates,
        primary_reason=primary,
        winner=winner,
        result=result,
    )


def checked_root(root_fen: str) -> chess.Board:
    if not isinstance(root_fen, str) or len(root_fen.split()) != 6:
        raise ReferenceError("KREF_E_FEN_FIELDS", "root_fen must contain six fields")
    try:
        board = chess.Board(root_fen, chess960=False)
    except ValueError as exc:
        raise ReferenceError("KREF_E_FEN_PARSE", str(exc)) from exc

    if board.status() != chess.STATUS_VALID:
        raise ReferenceError("KREF_E_FEN_INVALID", f"status={board.status()}")
    if canonical_fen(board) != root_fen:
        raise ReferenceError("KREF_E_FEN_NONCANONICAL", "root_fen does not round-trip exactly")

    occupants = goal_occupants(board)
    if occupants:
        raise ReferenceError(
            "KREF_E_AMBIGUOUS_GOAL_ROOT",
            "goal occupancy lacks an authenticated predecessor transition",
        )

    status = terminal_status(board, hill_transition=False, mover=None)
    if status.terminal:
        raise ReferenceError(
            "KREF_E_TERMINAL_ROOT",
            f"already-terminal root: {status.primary_reason}",
        )
    return board


def koth_san(standard_san: str, status: TerminalStatus) -> str:
    if "CHECKMATE" in status.predicates or "HILL" in status.predicates:
        if standard_san.endswith(("+", "#")):
            standard_san = standard_san[:-1]
        return standard_san + "#"
    return standard_san


def replay(request: dict[str, Any]) -> dict[str, Any]:
    if request.get("schema") != INPUT_SCHEMA:
        raise ReferenceError("KREF_E_SCHEMA", f"expected {INPUT_SCHEMA}")
    if request.get("rule_profile") != PROFILE:
        raise ReferenceError("KREF_E_PROFILE", f"expected {PROFILE}")

    moves = request.get("moves")
    if not isinstance(moves, list) or any(not isinstance(move, str) for move in moves):
        raise ReferenceError("KREF_E_MOVES", "moves must be a JSON array of strings")

    root_fen = request.get("root_fen")
    board = checked_root(root_fen)
    records: list[dict[str, Any]] = []
    last_status = terminal_status(board, hill_transition=False, mover=None)

    for index, token in enumerate(moves, start=1):
        if last_status.terminal:
            raise ReferenceError(
                "KREF_E_MOVE_AFTER_TERMINAL",
                f"move token follows terminal ply {index - 1}",
            )
        try:
            move = chess.Move.from_uci(token)
        except ValueError as exc:
            raise ReferenceError("KREF_E_MOVE_PARSE", f"ply={index}") from exc
        if move not in board.legal_moves:
            raise ReferenceError("KREF_E_MOVE_ILLEGAL", f"ply={index} move={token}")

        pre_fen = canonical_fen(board)
        moving_piece = board.piece_at(move.from_square)
        mover = board.turn
        standard_san = board.san(move)
        board.push(move)

        hill_transition = bool(
            moving_piece
            and moving_piece.piece_type == chess.KING
            and moving_piece.color == mover
            and move.to_square in GOALS
            and board.king(mover) == move.to_square
        )
        if len(goal_occupants(board)) > 1:
            raise ReferenceError("KREF_E_INVALID_TRAJECTORY", f"ply={index}")

        last_status = terminal_status(
            board,
            hill_transition=hill_transition,
            mover=mover,
        )
        physical_moves = sorted_legal_uci(board)
        records.append(
            {
                "ply": index,
                "move": token,
                "pre_fen": pre_fen,
                "post_fen": canonical_fen(board),
                "san": koth_san(standard_san, last_status),
                "physical_legal_moves_sha256": legal_moves_sha256(physical_moves),
                "physical_legal_move_count": len(physical_moves),
                "game_legal_move_count": 0 if last_status.terminal else len(physical_moves),
                "terminal": asdict(last_status),
            }
        )

    return {
        "schema": OUTPUT_SCHEMA,
        "reference_id": REFERENCE_ID,
        "rule_profile": PROFILE,
        "accepted": True,
        "root_fen": root_fen,
        "moves": moves,
        "plies": records,
        "final_fen": canonical_fen(board),
        "terminal": asdict(last_status),
    }


def failure(exc: ReferenceError) -> dict[str, Any]:
    return {
        "schema": OUTPUT_SCHEMA,
        "reference_id": REFERENCE_ID,
        "rule_profile": PROFILE,
        "accepted": False,
        "error": {"code": exc.code, "detail": exc.detail},
    }


def identity() -> dict[str, Any]:
    return {
        "schema": "koth-reference-identity-v1",
        "reference_id": REFERENCE_ID,
        "rule_profile": PROFILE,
        "python_version": platform.python_version(),
        "python_implementation": platform.python_implementation(),
        "chess_version": chess.__version__,
        "goal_squares": ["d4", "e4", "d5", "e5"],
    }


def read_request() -> dict[str, Any]:
    try:
        value = json.loads(sys.stdin.buffer.read().decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ReferenceError("KREF_E_JSON", "invalid UTF-8 JSON input") from exc
    if not isinstance(value, dict):
        raise ReferenceError("KREF_E_JSON_TYPE", "input must be a JSON object")
    return value


def emit(value: dict[str, Any]) -> None:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    sys.stdout.buffer.write(payload + b"\n")
    sys.stdout.buffer.flush()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("identity", "replay"))
    args = parser.parse_args(argv)

    if args.command == "identity":
        emit(identity())
        return 0

    try:
        emit(replay(read_request()))
        return 0
    except ReferenceError as exc:
        emit(failure(exc))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
