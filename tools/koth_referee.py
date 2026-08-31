#!/usr/bin/env python3
"""Independent, trajectory-authenticated KOTH_LICHESS_V1 referee.

The referee deliberately uses python-chess only for orthodox physical move
legality and notation. KOTH terminal timing, predicate precedence, automatic
draw policy, clock ordering, persistence, and fail-closed root handling are
implemented here from the project authority contract.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import chess
import chess.pgn
import chess.variant


PROFILE = "KOTH_LICHESS_V1"
CLOCK_PROFILE = "KOTH_CLOCK_V1"
SCHEMA = "KOTH_REFEREE_RECORD_V1"
REFEREE_VERSION = 1
PYTHON_CHESS_VERSION = "1.11.2"
PYTHON_CHESS_COMMIT = "9c24454dcea4f8a30259d811a2f10b26e911deb4"
PYTHON_CHESS_HASHES = {
    "__init__.py": "1286FF19809DF5AD3B649493C8B2AFB516AAFB1D9C7784B6E22E6C3D57E06E1A",
    "variant.py": "B7A9A5DD62CC5CCD810D4F64B70E4959CD4B1ACF927AA5DDC1E460DFD180B57B",
}

GOAL_SQUARES = frozenset(
    (chess.D4, chess.E4, chess.D5, chess.E5)
)

PREDICATE_BITS = {
    "CHECKMATE": 1 << 0,
    "HILL": 1 << 1,
    "STALEMATE": 1 << 2,
    "REPETITION3_AUTO": 1 << 3,
    "REPETITION5_DIAGNOSTIC": 1 << 4,
    "RULE50_AUTO": 1 << 5,
}


class RefereeError(RuntimeError):
    """A fail-closed referee error with a stable machine code."""

    def __init__(self, code: str, detail: str = "") -> None:
        self.code = code
        self.detail = detail
        super().__init__(f"{code}{': ' + detail if detail else ''}")


@dataclass(frozen=True)
class ClockConfig:
    initial_ms: int
    increment_ms: int = 0

    def __post_init__(self) -> None:
        if self.initial_ms <= 0:
            raise RefereeError("INVALID_CLOCK_INITIAL", str(self.initial_ms))
        if self.increment_ms < 0:
            raise RefereeError("INVALID_CLOCK_INCREMENT", str(self.increment_ms))


@dataclass(frozen=True)
class MoveTiming:
    elapsed_ms: int
    lag_compensation_ms: int = 0

    def __post_init__(self) -> None:
        if self.elapsed_ms < 0:
            raise RefereeError("INVALID_ELAPSED_TIME", str(self.elapsed_ms))
        if self.lag_compensation_ms < 0 or self.lag_compensation_ms > self.elapsed_ms:
            raise RefereeError(
                "INVALID_LAG_COMPENSATION",
                f"elapsed={self.elapsed_ms} lag={self.lag_compensation_ms}",
            )

    @property
    def charge_ms(self) -> int:
        return self.elapsed_ms - self.lag_compensation_ms


@dataclass(frozen=True)
class BoardStatus:
    predicates: tuple[str, ...]
    bitmask: int
    primary: str
    winner: str | None

    @property
    def terminal(self) -> bool:
        return self.primary != "NONE"

    def as_dict(self) -> dict[str, Any]:
        return {
            "terminal": self.terminal,
            "predicates": list(self.predicates),
            "bitmask": self.bitmask,
            "primary": self.primary,
            "winner": self.winner,
        }


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def dependency_identity() -> dict[str, Any]:
    package_root = Path(chess.__file__).resolve().parent
    return {
        "name": "python-chess",
        "version": chess.__version__,
        "commit": PYTHON_CHESS_COMMIT,
        "license": "GPL-3.0+",
        "files": {
            "__init__.py": _sha256(package_root / "__init__.py"),
            "variant.py": _sha256(package_root / "variant.py"),
        },
    }


def verify_dependency_identity() -> dict[str, Any]:
    identity = dependency_identity()
    if identity["version"] != PYTHON_CHESS_VERSION:
        raise RefereeError(
            "REFERENCE_VERSION_MISMATCH",
            f"expected={PYTHON_CHESS_VERSION} actual={identity['version']}",
        )
    for name, expected in PYTHON_CHESS_HASHES.items():
        actual = identity["files"][name]
        if actual != expected:
            raise RefereeError(
                "REFERENCE_SOURCE_MISMATCH",
                f"file={name} expected={expected} actual={actual}",
            )
    return identity


def _color_name(color: chess.Color) -> str:
    return "white" if color == chess.WHITE else "black"


def _result(winner: str | None, terminal: bool) -> str:
    if not terminal:
        return "*"
    if winner == "white":
        return "1-0"
    if winner == "black":
        return "0-1"
    return "1/2-1/2"


def _canonical_json(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")


def _canonical_fen(fen: str) -> tuple[str, chess.Board]:
    fields = fen.split()
    if len(fields) != 6:
        raise RefereeError(
            "INVALID_FEN_FIELD_COUNT", f"expected=6 actual={len(fields)}"
        )
    canonical = " ".join(fields)
    try:
        board = chess.Board(canonical)
    except ValueError as error:
        raise RefereeError("INVALID_FEN", str(error)) from error

    if board.status() != chess.STATUS_VALID:
        raise RefereeError("INVALID_FEN_STATUS", str(board.status()))
    emitted = board.fen(en_passant="fen")
    if emitted != canonical:
        raise RefereeError(
            "NONCANONICAL_OR_INCONSISTENT_FEN",
            f"input={canonical} emitted={emitted}",
        )
    return canonical, board


def _raw_fen(board: chess.Board) -> str:
    return board.fen(en_passant="fen")


def _epd(board: chess.Board) -> str:
    return " ".join(_raw_fen(board).split()[:4])


def _occurrence_key(board: chess.Board) -> str:
    return " ".join(board.fen(en_passant="legal").split()[:4])


def _king_on_goal(board: chess.Board) -> bool:
    return any(board.piece_at(square) and board.piece_at(square).piece_type == chess.KING
               for square in GOAL_SQUARES)


def _parse_uci(token: str) -> chess.Move:
    if len(token) not in (4, 5) or token != token.lower():
        raise RefereeError("NONCANONICAL_UCI_MOVE", token)
    try:
        move = chess.Move.from_uci(token)
    except ValueError as error:
        raise RefereeError("NONCANONICAL_UCI_MOVE", token) from error
    if not move or move.uci() != token:
        raise RefereeError("NONCANONICAL_UCI_MOVE", token)
    return move


class Referee:
    """Transactional KOTH referee for one root-plus-trajectory game."""

    def __init__(
        self,
        root_fen: str = chess.STARTING_FEN,
        clock: ClockConfig | None = None,
        *,
        verify_reference: bool = True,
    ) -> None:
        self.reference = (
            verify_dependency_identity() if verify_reference else dependency_identity()
        )
        self.root_fen, self.board = _canonical_fen(root_fen)
        self.notation_board = chess.variant.KingOfTheHillBoard(self.root_fen)
        self.clock = clock
        self.clocks = (
            {"white": clock.initial_ms, "black": clock.initial_ms}
            if clock
            else None
        )
        self.accepted_moves: list[dict[str, Any]] = []
        self.forfeit_event: dict[str, Any] | None = None
        self.occurrences = {_occurrence_key(self.board): 1}
        self.status = self._classify(hill_transition=False)
        self.referee_termination = "BOARD" if self.status.terminal else "NONE"
        self.referee_winner = self.status.winner

        if _king_on_goal(self.board):
            raise RefereeError("AMBIGUOUS_GOAL_ROOT")
        if self.status.terminal:
            raise RefereeError("TERMINAL_RAW_ROOT", self.status.primary)

    @property
    def terminal(self) -> bool:
        return self.referee_termination != "NONE"

    @property
    def winner(self) -> str | None:
        return self.referee_winner

    @property
    def result(self) -> str:
        return _result(self.winner, self.terminal)

    def physical_moves(self) -> tuple[str, ...]:
        return tuple(sorted(move.uci() for move in self.board.legal_moves))

    def game_moves(self) -> tuple[str, ...]:
        return () if self.terminal else self.physical_moves()

    def _classify(self, *, hill_transition: bool) -> BoardStatus:
        predicates: list[str] = []
        if self.board.is_checkmate():
            predicates.append("CHECKMATE")
        if hill_transition:
            predicates.append("HILL")
        if self.board.is_stalemate():
            predicates.append("STALEMATE")

        count = self.occurrences.get(_occurrence_key(self.board), 1)
        if count >= 3:
            predicates.append("REPETITION3_AUTO")
        if count >= 5:
            predicates.append("REPETITION5_DIAGNOSTIC")
        if self.board.halfmove_clock >= 100:
            predicates.append("RULE50_AUTO")

        predicate_set = set(predicates)
        if "CHECKMATE" in predicate_set:
            primary = "CHECKMATE"
        elif "HILL" in predicate_set:
            primary = "HILL"
        elif "STALEMATE" in predicate_set:
            primary = "STALEMATE"
        elif "REPETITION3_AUTO" in predicate_set or "RULE50_AUTO" in predicate_set:
            primary = "AUTOMATIC_DRAW"
        else:
            primary = "NONE"

        winner = (
            _color_name(not self.board.turn)
            if primary in ("CHECKMATE", "HILL")
            else None
        )
        bitmask = sum(PREDICATE_BITS[name] for name in predicates)
        return BoardStatus(tuple(predicates), bitmask, primary, winner)

    def submit(self, uci: str, timing: MoveTiming | None = None) -> BoardStatus:
        if self.terminal:
            raise RefereeError("POSTTERMINAL_MOVE", uci)
        if (self.clock is None) != (timing is None):
            raise RefereeError(
                "CLOCK_TIMING_CONTRACT",
                "timing is required exactly when a clock is configured",
            )

        move = _parse_uci(uci)
        if not self.board.is_legal(move):
            raise RefereeError("ILLEGAL_MOVE", uci)

        mover = self.board.turn
        mover_name = _color_name(mover)
        pre_fen = _raw_fen(self.board)
        pre_epd = _epd(self.board)
        clocks_before = dict(self.clocks) if self.clocks is not None else None

        if self.clocks is not None and timing is not None:
            remaining = self.clocks[mover_name] - timing.charge_ms
            if remaining <= 0:
                self.clocks[mover_name] = max(0, remaining)
                self.referee_termination = "TIME_FORFEIT"
                self.referee_winner = _color_name(not mover)
                self.forfeit_event = {
                    "attempted_uci": uci,
                    "mover": mover_name,
                    "elapsed_ms": timing.elapsed_ms,
                    "lag_compensation_ms": timing.lag_compensation_ms,
                    "charged_ms": timing.charge_ms,
                    "clocks_before_ms": clocks_before,
                    "clocks_after_ms": dict(self.clocks),
                    "board_progress_saved": False,
                    "pre_fen": pre_fen,
                }
                return self.status

        candidate = self.board.copy(stack=True)
        notation_candidate = self.notation_board.copy(stack=True)
        moved_piece = candidate.piece_at(move.from_square)
        hill_transition = bool(
            moved_piece
            and moved_piece.piece_type == chess.KING
            and move.to_square in GOAL_SQUARES
            and move.from_square not in GOAL_SQUARES
        )
        san = notation_candidate.san(move)
        candidate.push(move)
        notation_candidate.push(move)
        if _raw_fen(candidate) != _raw_fen(notation_candidate):
            raise RefereeError("NOTATION_BOARD_DIVERGENCE")

        occurrence_key = _occurrence_key(candidate)
        occurrence_count = self.occurrences.get(occurrence_key, 0) + 1
        candidate_occurrences = dict(self.occurrences)
        candidate_occurrences[occurrence_key] = occurrence_count

        old_board = self.board
        old_occurrences = self.occurrences
        self.board = candidate
        self.occurrences = candidate_occurrences
        candidate_status = self._classify(hill_transition=hill_transition)
        self.board = old_board
        self.occurrences = old_occurrences

        if self.clocks is not None and timing is not None and self.clock is not None:
            self.clocks[mover_name] = (
                self.clocks[mover_name] - timing.charge_ms + self.clock.increment_ms
            )

        clocks_after = dict(self.clocks) if self.clocks is not None else None
        move_record = {
            "ply": len(self.accepted_moves) + 1,
            "mover": mover_name,
            "uci": uci,
            "san": san,
            "pre_fen": pre_fen,
            "pre_epd": pre_epd,
            "post_fen": _raw_fen(candidate),
            "post_epd": _epd(candidate),
            "position_occurrence_count": occurrence_count,
            "hill_transition": hill_transition,
            "board_status": candidate_status.as_dict(),
            "timing": (
                {
                    "elapsed_ms": timing.elapsed_ms,
                    "lag_compensation_ms": timing.lag_compensation_ms,
                    "charged_ms": timing.charge_ms,
                }
                if timing is not None
                else None
            ),
            "clocks_before_ms": clocks_before,
            "clocks_after_ms": clocks_after,
        }

        self.board = candidate
        self.notation_board = notation_candidate
        self.occurrences = candidate_occurrences
        self.accepted_moves.append(move_record)
        self.status = candidate_status
        if self.status.terminal:
            self.referee_termination = "BOARD"
            self.referee_winner = self.status.winner
        return self.status

    def trajectory_id(self) -> str:
        identity = {
            "profile": PROFILE,
            "root_fen": self.root_fen,
            "accepted_uci": [move["uci"] for move in self.accepted_moves],
        }
        return hashlib.sha256(_canonical_json(identity)).hexdigest().upper()

    def record(self) -> dict[str, Any]:
        record: dict[str, Any] = {
            "schema": SCHEMA,
            "schema_version": REFEREE_VERSION,
            "profile": PROFILE,
            "clock_profile": CLOCK_PROFILE,
            "goal_squares": ["d4", "e4", "d5", "e5"],
            "reference": self.reference,
            "root_fen": self.root_fen,
            "accepted_moves": self.accepted_moves,
            "accepted_uci": [move["uci"] for move in self.accepted_moves],
            "final_fen": _raw_fen(self.board),
            "final_epd": _epd(self.board),
            "physical_moves": list(self.physical_moves()),
            "game_moves": list(self.game_moves()),
            "board_status": self.status.as_dict(),
            "referee_termination": self.referee_termination,
            "winner": self.winner,
            "result": self.result,
            "clock": (
                {
                    "initial_ms": self.clock.initial_ms,
                    "increment_ms": self.clock.increment_ms,
                    "remaining_ms": self.clocks,
                }
                if self.clock is not None
                else None
            ),
            "forfeit_event": self.forfeit_event,
            "trajectory_id": self.trajectory_id(),
        }
        record["record_sha256"] = hashlib.sha256(_canonical_json(record)).hexdigest().upper()
        return record

    def pgn(self, headers: Mapping[str, str] | None = None) -> str:
        game = chess.pgn.Game()
        game.setup(chess.variant.KingOfTheHillBoard(self.root_fen))
        game.headers["Event"] = "KOTH-Stockfish certified referee"
        game.headers["Site"] = "local"
        game.headers["Date"] = "????.??.??"
        game.headers["Round"] = "-"
        game.headers["White"] = "White"
        game.headers["Black"] = "Black"
        game.headers["Result"] = self.result
        game.headers["Termination"] = (
            "Time forfeit" if self.referee_termination == "TIME_FORFEIT" else "Normal"
            if self.terminal else "Unterminated"
        )
        game.headers["KOTHProfile"] = PROFILE
        game.headers["KOTHClock"] = CLOCK_PROFILE
        game.headers["KOTHPredicates"] = ",".join(self.status.predicates) or "NONE"
        game.headers["KOTHPrimary"] = self.status.primary
        game.headers["KOTHTrajectory"] = self.trajectory_id()
        if headers:
            for name, value in headers.items():
                if not isinstance(name, str) or not isinstance(value, str):
                    raise RefereeError("INVALID_PGN_HEADER_TYPE")
                game.headers[name] = value

        node: chess.pgn.GameNode = game
        for move_record in self.accepted_moves:
            node = node.add_variation(chess.Move.from_uci(move_record["uci"]))

        exporter = chess.pgn.StringExporter(
            headers=True, variations=False, comments=False
        )
        return game.accept(exporter) + "\n"


def replay(
    root_fen: str,
    moves: Iterable[Mapping[str, Any]],
    clock: ClockConfig | None = None,
) -> Referee:
    referee = Referee(root_fen, clock)
    for item in moves:
        if "uci" not in item:
            raise RefereeError("MISSING_MOVE_UCI")
        timing = None
        if clock is not None:
            timing = MoveTiming(
                int(item.get("elapsed_ms", -1)),
                int(item.get("lag_compensation_ms", 0)),
            )
        referee.submit(str(item["uci"]), timing)
    return referee


def _clock_from_payload(payload: Mapping[str, Any]) -> ClockConfig | None:
    raw = payload.get("clock")
    if raw is None:
        return None
    if not isinstance(raw, Mapping):
        raise RefereeError("INVALID_CLOCK_PAYLOAD")
    return ClockConfig(int(raw["initial_ms"]), int(raw.get("increment_ms", 0)))


def run_payload(payload: Mapping[str, Any]) -> dict[str, Any]:
    root_fen = str(payload.get("root_fen", chess.STARTING_FEN))
    raw_moves = payload.get("moves", [])
    if not isinstance(raw_moves, list) or not all(isinstance(item, Mapping) for item in raw_moves):
        raise RefereeError("INVALID_MOVES_PAYLOAD")
    referee = replay(root_fen, raw_moves, _clock_from_payload(payload))
    record = referee.record()
    pgn = referee.pgn(payload.get("pgn_headers"))
    return {
        "record": record,
        "pgn": pgn,
        "pgn_sha256": hashlib.sha256(pgn.encode("utf-8")).hexdigest().upper(),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        help="JSON request file; stdin is used when omitted",
    )
    args = parser.parse_args(argv)

    try:
        text = args.input.read_text(encoding="utf-8") if args.input else sys.stdin.read()
        payload = json.loads(text)
        if not isinstance(payload, Mapping):
            raise RefereeError("INVALID_TOP_LEVEL_PAYLOAD")
        print(json.dumps(run_payload(payload), ensure_ascii=False, sort_keys=True, indent=2))
        return 0
    except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError, RefereeError) as error:
        if isinstance(error, RefereeError):
            code, detail = error.code, error.detail
        else:
            code, detail = "INVALID_REQUEST", str(error)
        print(
            json.dumps({"error": {"code": code, "detail": detail}}, sort_keys=True),
            file=sys.stderr,
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
