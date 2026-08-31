#!/usr/bin/env python3
"""Independent referee, persistence, clock, and engine differential tests."""

from __future__ import annotations

import argparse
import io
import json
import random
import subprocess
import sys
from pathlib import Path
from typing import Callable

import chess
import chess.pgn


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import koth_referee as kr  # noqa: E402


class TestFailure(RuntimeError):
    pass


class EngineClient:
    def __init__(self, executable: Path) -> None:
        self.process = subprocess.Popen(
            [str(executable)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )
        if self.process.stdin is None or self.process.stdout is None:
            raise TestFailure("could not open engine pipes")
        self.stdin = self.process.stdin
        self.stdout = self.process.stdout
        self.transact("uci", terminal="uciok")

    def transact(self, *commands: str, terminal: str = "readyok") -> list[str]:
        for command in commands:
            self.stdin.write(command + "\n")
        if terminal == "readyok":
            self.stdin.write("isready\n")
        self.stdin.flush()
        lines: list[str] = []
        while True:
            line = self.stdout.readline()
            if not line:
                raise TestFailure(f"engine exited before {terminal}: {lines!r}")
            line = line.rstrip("\r\n")
            lines.append(line)
            if line == terminal:
                return lines

    def snapshot(self, root_fen: str, moves: list[str]) -> dict[str, object]:
        command = f"position fen {root_fen}"
        if moves:
            command += " moves " + " ".join(moves)
        lines = self.transact(command, "d", "kothstatus", "kothmoves")
        if any("error command=position" in line for line in lines):
            raise TestFailure("engine rejected referee trajectory: " + "\n".join(lines))

        fen_lines = [line[5:] for line in lines if line.startswith("Fen: ")]
        status_lines = [
            line.removeprefix("info string koth ")
            for line in lines
            if line.startswith("info string koth profile=")
        ]
        move_lines = [
            line.removeprefix("info string koth ")
            for line in lines
            if line.startswith("info string koth physical_count=")
        ]
        if len(fen_lines) != 1 or len(status_lines) != 1 or len(move_lines) != 1:
            raise TestFailure(f"incomplete engine snapshot: {lines!r}")

        status = dict(part.split("=", 1) for part in status_lines[0].split())
        physical_field = move_lines[0].split(" physical=", 1)[1].split(" game=", 1)[0]
        game_field = move_lines[0].split(" game=", 1)[1]
        return {
            "fen": fen_lines[0],
            "bitmask": int(status["predicates"]),
            "primary": status["primary"],
            "winner": status["winner"],
            "terminal": status["terminal"] == "true",
            "physical": tuple(sorted(move for move in physical_field.split(",") if move)),
            "game": tuple(sorted(move for move in game_field.split(",") if move)),
        }

    def close(self) -> None:
        if self.process.poll() is None:
            self.stdin.write("quit\n")
            self.stdin.flush()
        try:
            code = self.process.wait(timeout=10)
        except subprocess.TimeoutExpired as error:
            self.process.kill()
            raise TestFailure("engine did not stop") from error
        if code != 0:
            raise TestFailure(f"engine exit code {code}")


class Suite:
    def __init__(self, engine: EngineClient) -> None:
        self.engine = engine
        self.assertions = 0
        self.differential_positions = 0

    def check(self, condition: bool, message: str) -> None:
        self.assertions += 1
        if not condition:
            raise TestFailure(message)

    def expect_error(self, code: str, action: Callable[[], object]) -> None:
        try:
            action()
        except kr.RefereeError as error:
            self.check(error.code == code, f"expected {code}, got {error}")
        else:
            raise TestFailure(f"expected referee error {code}")

    def compare_engine(self, referee: kr.Referee) -> None:
        moves = [item["uci"] for item in referee.accepted_moves]
        observed = self.engine.snapshot(referee.root_fen, moves)
        self.differential_positions += 1
        self.check(observed["fen"] == referee.record()["final_fen"], str(observed))
        self.check(observed["bitmask"] == referee.status.bitmask, str(observed))
        self.check(observed["primary"] == referee.status.primary, str(observed))
        self.check(observed["winner"] == (referee.status.winner or "none"), str(observed))
        self.check(observed["terminal"] == referee.status.terminal, str(observed))
        self.check(observed["physical"] == referee.physical_moves(), str(observed))
        self.check(observed["game"] == referee.game_moves(), str(observed))

    def test_reference_identity(self) -> None:
        identity = kr.verify_dependency_identity()
        self.check(identity["version"] == "1.11.2", str(identity))
        self.check(identity["commit"] == kr.PYTHON_CHESS_COMMIT, str(identity))
        original = kr.PYTHON_CHESS_HASHES["variant.py"]
        try:
            kr.PYTHON_CHESS_HASHES["variant.py"] = "0" * 64
            self.expect_error("REFERENCE_SOURCE_MISMATCH", kr.verify_dependency_identity)
        finally:
            kr.PYTHON_CHESS_HASHES["variant.py"] = original

    def test_goal_and_precedence(self) -> None:
        legal_entries = [
            ("7k/8/8/8/8/2K5/8/8 w - - 0 1", "c3d4", "white"),
            ("k7/8/8/8/8/5K2/8/8 w - - 0 1", "f3e4", "white"),
            ("8/8/2K5/8/8/8/8/7k w - - 0 1", "c6d5", "white"),
            ("8/8/5K2/8/8/8/8/k7 w - - 0 1", "f6e5", "white"),
            ("8/8/8/2k5/8/8/8/7K b - - 0 1", "c5d4", "black"),
            ("8/8/8/5k2/8/8/8/K7 b - - 0 1", "f5e4", "black"),
            ("7K/8/8/8/2k5/8/8/8 b - - 0 1", "c4d5", "black"),
            ("K7/8/8/8/5k2/8/8/8 b - - 0 1", "f4e5", "black"),
        ]
        for fen, move, winner in legal_entries:
            referee = kr.Referee(fen)
            status = referee.submit(move)
            self.check(status.predicates == ("HILL",), str(status))
            self.check(status.winner == winner, str(status))
            self.check(referee.game_moves() == (), str(referee.game_moves()))
            self.check(bool(referee.physical_moves()), "physical domain disappeared")
            self.compare_engine(referee)

        cases = [
            (
                "1rkr4/1p1p4/8/8/8/2K5/8/2R5 w - - 0 1",
                "c3d4",
                ("CHECKMATE", "HILL"),
                3,
                "CHECKMATE",
            ),
            (
                "k7/2Q5/8/8/8/2K5/8/8 w - - 0 1",
                "c3d4",
                ("HILL", "STALEMATE"),
                6,
                "HILL",
            ),
            (
                "7k/8/8/8/8/2K5/8/8 w - - 99 1",
                "c3d4",
                ("HILL", "RULE50_AUTO"),
                34,
                "HILL",
            ),
            (
                "7k/8/5KQ1/8/8/8/8/8 w - - 0 1",
                "g6g7",
                ("CHECKMATE",),
                1,
                "CHECKMATE",
            ),
            (
                "k7/8/2K5/8/1Q6/8/8/8 w - - 0 1",
                "b4b6",
                ("STALEMATE",),
                4,
                "STALEMATE",
            ),
        ]
        for fen, move, predicates, bitmask, primary in cases:
            referee = kr.Referee(fen)
            status = referee.submit(move)
            self.check(status.predicates == predicates, str(status))
            self.check(status.bitmask == bitmask, str(status))
            self.check(status.primary == primary, str(status))
            self.compare_engine(referee)

    def test_rejections_and_rollback(self) -> None:
        for fen in (
            "7k/8/8/8/3K4/8/8/8 b - - 0 1",
            "8/8/8/4k3/8/8/8/K7 w - - 0 1",
        ):
            self.expect_error("AMBIGUOUS_GOAL_ROOT", lambda fen=fen: kr.Referee(fen))

        for fen in (
            "7k/6Q1/6K1/8/8/8/8/8 b - - 0 1",
            "k7/2Q5/8/8/8/8/8/K7 b - - 0 1",
            "7k/8/8/8/8/2K5/8/8 w - - 100 1",
        ):
            self.expect_error("TERMINAL_RAW_ROOT", lambda fen=fen: kr.Referee(fen))

        attacked = kr.Referee("3r3k/8/8/8/8/2K5/8/8 w - - 0 1")
        before = attacked.record()
        self.expect_error("ILLEGAL_MOVE", lambda: attacked.submit("c3d4"))
        self.check(attacked.record() == before, "illegal move mutated referee state")

        attacked_entries = [
            ("3r3k/8/8/8/8/2K5/8/8 w - - 0 1", "c3d4"),
            ("k3r3/8/8/8/8/5K2/8/8 w - - 0 1", "f3e4"),
            ("8/8/2K5/8/8/8/8/3r3k w - - 0 1", "c6d5"),
            ("8/8/5K2/8/8/8/8/k3r3 w - - 0 1", "f6e5"),
            ("3R3K/8/8/8/8/2k5/8/8 b - - 0 1", "c3d4"),
            ("K3R3/8/8/8/8/5k2/8/8 b - - 0 1", "f3e4"),
            ("8/8/2k5/8/8/8/8/3R3K b - - 0 1", "c6d5"),
            ("8/8/5k2/8/8/8/8/K3R3 b - - 0 1", "f6e5"),
        ]
        for fen, move in attacked_entries:
            referee = kr.Referee(fen)
            snapshot = referee.record()
            self.expect_error("ILLEGAL_MOVE", lambda referee=referee, move=move: referee.submit(move))
            self.check(referee.record() == snapshot, f"attacked entry mutated state: {fen} {move}")

        for move in ("C3D4", "c3d9", "0000"):
            referee = kr.Referee("7k/8/8/8/8/2K5/8/8 w - - 0 1")
            self.expect_error("NONCANONICAL_UCI_MOVE", lambda move=move: referee.submit(move))

        promoted_king = kr.Referee("7k/8/8/8/8/2K5/8/8 w - - 0 1")
        self.expect_error("ILLEGAL_MOVE", lambda: promoted_king.submit("c3d4q"))

        terminal = kr.Referee("7k/8/8/8/8/2K5/8/8 w - - 0 1")
        terminal.submit("c3d4")
        self.expect_error("POSTTERMINAL_MOVE", lambda: terminal.submit("h8h7"))

        self.expect_error(
            "INVALID_FEN_STATUS",
            lambda: kr.Referee(
                "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq e3 0 1"
            ),
        )
        self.expect_error(
            "NONCANONICAL_OR_INCONSISTENT_FEN",
            lambda: kr.Referee(
                "rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w QKkq - 0 1"
            ),
        )

    def test_draws_and_persistence(self) -> None:
        referee = kr.Referee()
        moves = ["g1f3", "g8f6", "f3g1", "f6g8"] * 2
        for move in moves:
            referee.submit(move)
        self.check(referee.status.predicates == ("REPETITION3_AUTO",), str(referee.status))
        self.check(referee.status.primary == "AUTOMATIC_DRAW", str(referee.status))
        self.check(referee.result == "1/2-1/2", referee.result)
        self.compare_engine(referee)

        source = kr.Referee("7k/8/8/8/8/2K5/8/8 w - - 0 1")
        source.submit("c3d4")
        record = source.record()
        clone = kr.replay(record["root_fen"], [{"uci": move} for move in record["accepted_uci"]])
        self.check(clone.record() == record, "canonical record replay changed bytes")
        self.check(record["final_epd"] == " ".join(record["final_fen"].split()[:4]), str(record))
        without_hash = dict(record)
        observed_hash = without_hash.pop("record_sha256")
        expected_hash = kr.hashlib.sha256(kr._canonical_json(without_hash)).hexdigest().upper()
        self.check(observed_hash == expected_hash, "record hash does not authenticate payload")

        pgn = source.pgn()
        self.check('[Variant "King of the Hill"]' in pgn, pgn)
        self.check('[Termination "Normal"]' in pgn, pgn)
        self.check('[KOTHPredicates "HILL"]' in pgn, pgn)
        self.check("1. Kd4# 1-0" in pgn, pgn)
        parsed = chess.pgn.read_game(io.StringIO(pgn))
        self.check(parsed is not None, pgn)
        assert parsed is not None
        self.check(parsed.headers["Result"] == "1-0", str(parsed.headers))
        self.check([move.uci() for move in parsed.mainline_moves()] == ["c3d4"], pgn)

        output = kr.run_payload(
            {
                "root_fen": source.root_fen,
                "moves": [{"uci": "c3d4"}],
                "pgn_headers": {"White": "Candidate", "Black": "Comparator"},
            }
        )
        self.check(output["record"] == record, json.dumps(output, sort_keys=True))
        self.check(output["pgn_sha256"] == kr.hashlib.sha256(output["pgn"].encode()).hexdigest().upper(), str(output))

    def test_clock_ordering(self) -> None:
        root = "7k/8/8/8/8/2K5/8/8 w - - 0 1"

        boundary = kr.Referee(root, kr.ClockConfig(100, 0))
        before = boundary.record()["final_fen"]
        boundary.submit("c3d4", kr.MoveTiming(100, 0))
        self.check(boundary.referee_termination == "TIME_FORFEIT", str(boundary.record()))
        self.check(boundary.result == "0-1", boundary.result)
        self.check(boundary.record()["final_fen"] == before, str(boundary.record()))
        self.check(boundary.accepted_moves == [], str(boundary.accepted_moves))
        self.check(boundary.forfeit_event is not None, str(boundary.record()))
        assert boundary.forfeit_event is not None
        self.check(not boundary.forfeit_event["board_progress_saved"], str(boundary.forfeit_event))
        self.check('[Termination "Time forfeit"]' in boundary.pgn(), boundary.pgn())

        compensated = kr.Referee(root, kr.ClockConfig(100, 0))
        compensated.submit("c3d4", kr.MoveTiming(110, 20))
        self.check(compensated.referee_termination == "BOARD", str(compensated.record()))
        self.check(compensated.result == "1-0", compensated.result)
        self.check(compensated.clocks == {"white": 10, "black": 100}, str(compensated.clocks))

        late = kr.Referee(root, kr.ClockConfig(100, 0))
        late.submit("c3d4", kr.MoveTiming(120, 10))
        self.check(late.referee_termination == "TIME_FORFEIT", str(late.record()))

        increment = kr.Referee(clock=kr.ClockConfig(100, 10))
        increment.submit("e2e4", kr.MoveTiming(30, 0))
        self.check(increment.clocks == {"white": 80, "black": 100}, str(increment.clocks))

        no_clock = kr.Referee()
        self.expect_error(
            "CLOCK_TIMING_CONTRACT",
            lambda: no_clock.submit("e2e4", kr.MoveTiming(1, 0)),
        )
        with_clock = kr.Referee(clock=kr.ClockConfig(100))
        self.expect_error("CLOCK_TIMING_CONTRACT", lambda: with_clock.submit("e2e4"))
        self.expect_error(
            "INVALID_LAG_COMPENSATION", lambda: kr.MoveTiming(10, 11)
        )

    def test_generated_differential(self) -> None:
        rng = random.Random(0x4B4F5448)
        for _ in range(8):
            referee = kr.Referee()
            self.compare_engine(referee)
            for _ply in range(24):
                if referee.terminal:
                    break
                moves = list(referee.board.legal_moves)
                self.check(bool(moves), "nonterminal state had no physical move")
                referee.submit(rng.choice(moves).uci())
                self.compare_engine(referee)

    def run(self) -> None:
        self.test_reference_identity()
        self.test_goal_and_precedence()
        self.test_rejections_and_rollback()
        self.test_draws_and_persistence()
        self.test_clock_ordering()
        self.test_generated_differential()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("engine", type=Path)
    args = parser.parse_args()
    executable = args.engine.resolve()
    if not executable.is_file():
        raise TestFailure(f"engine not found: {executable}")

    engine = EngineClient(executable)
    suite = Suite(engine)
    try:
        suite.run()
    finally:
        engine.close()
    print(
        "PASS koth_referee "
        f"assertions={suite.assertions} differential_positions={suite.differential_positions}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except TestFailure as error:
        print(f"FAIL {error}", file=sys.stderr)
        raise SystemExit(1)
