#!/usr/bin/env python3
"""Protocol-level conformance checks for the KOTH rule core.

This suite intentionally does not call the search evaluator. It checks the
authoritative root-plus-trajectory ingress, terminal classification, move-domain
split, repetition, FEN/EP handling, and fail-closed protocol behavior.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


class ConformanceFailure(RuntimeError):
    pass


class EngineSession:
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
            raise ConformanceFailure("failed to open engine pipes")

        self.stdin = self.process.stdin
        self.stdout = self.process.stdout
        self.uci_output = self.transact("uci", terminal="uciok")

    def transact(self, *commands: str, terminal: str = "readyok") -> list[str]:
        for command in commands:
            self.stdin.write(command + "\n")
        if terminal == "readyok":
            self.stdin.write("isready\n")
        self.stdin.flush()

        output: list[str] = []
        while True:
            line = self.stdout.readline()
            if line == "":
                raise ConformanceFailure(
                    f"engine exited before {terminal}; output={output!r}"
                )
            line = line.rstrip("\r\n")
            output.append(line)
            if line == terminal:
                return output

    def close(self) -> None:
        if self.process.poll() is None:
            try:
                self.stdin.write("quit\n")
                self.stdin.flush()
            except (BrokenPipeError, OSError):
                pass
        try:
            exit_code = self.process.wait(timeout=10)
        except subprocess.TimeoutExpired as error:
            self.process.kill()
            raise ConformanceFailure("engine did not exit after quit") from error
        if exit_code != 0:
            raise ConformanceFailure(f"engine exit code {exit_code}")


@dataclass(frozen=True)
class Snapshot:
    fen: str
    status: str
    moves: str


class Suite:
    def __init__(self, engine: EngineSession) -> None:
        self.engine = engine
        self.count = 0

    def check(self, condition: bool, message: str) -> None:
        self.count += 1
        if not condition:
            raise ConformanceFailure(message)

    @staticmethod
    def joined(lines: list[str]) -> str:
        return "\n".join(lines)

    def snapshot(self) -> Snapshot:
        output = self.engine.transact("d", "kothstatus", "kothmoves")
        fen_lines = [line for line in output if line.startswith("Fen: ")]
        status_lines = [
            line
            for line in output
            if line.startswith("info string koth profile=KOTH_LICHESS_V1")
        ]
        move_lines = [
            line
            for line in output
            if line.startswith("info string koth physical_count=")
        ]
        self.check(len(fen_lines) == 1, f"missing/duplicate FEN: {output!r}")
        self.check(len(status_lines) == 1, f"missing/duplicate status: {output!r}")
        self.check(len(move_lines) == 1, f"missing/duplicate move domains: {output!r}")
        return Snapshot(fen_lines[0][5:], status_lines[0], move_lines[0])

    def set_position(self, command: str) -> list[str]:
        return self.engine.transact(command)

    def assert_hill(self, fen: str, move: str, winner: str) -> Snapshot:
        output = self.set_position(f"position fen {fen} moves {move}")
        self.check("error command=position" not in self.joined(output), self.joined(output))
        state = self.snapshot()
        self.check("terminal=true" in state.status, state.status)
        self.check("predicates=2" in state.status, state.status)
        self.check("primary=HILL" in state.status, state.status)
        self.check(f"winner={winner}" in state.status, state.status)
        self.check("transition=ACCEPTED_MOVE" in state.status, state.status)
        self.check("king_entered_hill=true" in state.status, state.status)
        self.check("game_count=0" in state.moves, state.moves)
        return state

    def run(self) -> None:
        self.check(
            any(line.startswith("id name KOTH-Stockfish ") for line in self.engine.uci_output),
            self.joined(self.engine.uci_output),
        )
        self.check(
            any(
                line
                == "id author the KOTH-Stockfish developers; based on Stockfish (see AUTHORS file)"
                for line in self.engine.uci_output
            ),
            self.joined(self.engine.uci_output),
        )
        option_lines = [line for line in self.engine.uci_output if line.startswith("option name ")]
        exact_variant = (
            "option name UCI_Variant type combo default kingofthehill var kingofthehill"
        )
        self.check(exact_variant in option_lines, self.joined(option_lines))
        for forbidden in ("UCI_Chess960", "UCI_ShowWDL", "Syzygy", "EvalFile"):
            self.check(
                not any(forbidden in line for line in option_lines),
                f"forbidden option exposed: {forbidden}",
            )

        invalid_option = self.engine.transact(
            "setoption name UCI_Variant value chess"
        )
        self.check(
            any("code=INVALID_OPTION_VALUE" in line for line in invalid_option),
            self.joined(invalid_option),
        )
        valid_option = self.engine.transact(
            "setoption name UCI_Variant value kingofthehill"
        )
        self.check(
            not any("error" in line.lower() for line in valid_option),
            self.joined(valid_option),
        )

        white_entries = [
            ("7k/8/8/8/8/2K5/8/8 w - - 0 1", "c3d4"),
            ("k7/8/8/8/8/5K2/8/8 w - - 0 1", "f3e4"),
            ("8/8/2K5/8/8/8/8/7k w - - 0 1", "c6d5"),
            ("8/8/5K2/8/8/8/8/k7 w - - 0 1", "f6e5"),
        ]
        for fen, move in white_entries:
            self.assert_hill(fen, move, "white")

        black_entries = [
            ("8/8/8/2k5/8/8/8/7K b - - 0 1", "c5d4"),
            ("8/8/8/5k2/8/8/8/K7 b - - 0 1", "f5e4"),
            ("7K/8/8/8/2k5/8/8/8 b - - 0 1", "c4d5"),
            ("K7/8/8/8/5k2/8/8/8 b - - 0 1", "f4e5"),
        ]
        for fen, move in black_entries:
            self.assert_hill(fen, move, "black")

        attacked_entries = [
            ("3r3k/8/8/8/8/2K5/8/8 w - - 0 1", "c3d4"),
            ("4r2k/8/8/8/8/5K2/8/8 w - - 0 1", "f3e4"),
            ("8/8/2K5/8/8/8/8/3r3k w - - 0 1", "c6d5"),
            ("8/8/5K2/8/8/8/8/4r2k w - - 0 1", "f6e5"),
        ]
        for fen, move in attacked_entries:
            self.set_position("position startpos")
            before = self.snapshot()
            output = self.set_position(f"position fen {fen} moves {move}")
            self.check(
                any("code=ILLEGAL_TRAJECTORY_MOVE" in line for line in output),
                self.joined(output),
            )
            self.check(self.snapshot() == before, "illegal hill entry mutated live state")

        black_attacked_entries = [
            ("3R3K/8/8/8/8/2k5/8/8 b - - 0 1", "c3d4"),
            ("K3R3/8/8/8/8/5k2/8/8 b - - 0 1", "f3e4"),
            ("8/8/2k5/8/8/8/8/3R3K b - - 0 1", "c6d5"),
            ("8/8/5k2/8/8/8/8/K3R3 b - - 0 1", "f6e5"),
        ]
        for fen, move in black_attacked_entries:
            self.set_position("position startpos")
            before = self.snapshot()
            output = self.set_position(f"position fen {fen} moves {move}")
            self.check(
                any("code=ILLEGAL_TRAJECTORY_MOVE" in line for line in output),
                self.joined(output),
            )
            self.check(self.snapshot() == before, "illegal black hill entry mutated live state")

        output = self.set_position(
            "position fen 7k/8/8/8/2R5/8/8/K7 w - - 0 1 moves c4d4"
        )
        self.check("error command=position" not in self.joined(output), self.joined(output))
        nonking = self.snapshot()
        self.check("terminal=false" in nonking.status, nonking.status)
        self.check("king_entered_hill=false" in nonking.status, nonking.status)

        for goal_root in (
            "7k/8/8/8/3K4/8/8/8 b - - 0 1",
            "8/8/8/4k3/8/8/8/K7 w - - 0 1",
        ):
            self.set_position("position startpos")
            before = self.snapshot()
            output = self.set_position(f"position fen {goal_root}")
            self.check(
                any("code=AMBIGUOUS_GOAL_ROOT" in line for line in output),
                self.joined(output),
            )
            self.check(self.snapshot() == before, "goal root mutated live state")

        self.set_position("position startpos")
        before = self.snapshot()
        output = self.set_position(
            "position fen 7k/8/8/8/8/2K5/8/8 w - - 0 1 moves c3d4 h8h7"
        )
        self.check(
            any("code=POSTTERMINAL_TRAJECTORY_TAIL" in line for line in output),
            self.joined(output),
        )
        self.check(self.snapshot() == before, "postterminal tail mutated live state")

        output = self.set_position(
            "position fen 1rkr4/1p1p4/8/8/8/2K5/8/2R5 w - - 0 1 moves c3d4"
        )
        self.check("error command=position" not in self.joined(output), self.joined(output))
        mate_hill = self.snapshot()
        self.check("predicates=3" in mate_hill.status, mate_hill.status)
        self.check("primary=CHECKMATE" in mate_hill.status, mate_hill.status)

        output = self.set_position(
            "position fen k7/2Q5/8/8/8/2K5/8/8 w - - 0 1 moves c3d4"
        )
        self.check("error command=position" not in self.joined(output), self.joined(output))
        stale_hill = self.snapshot()
        self.check("predicates=6" in stale_hill.status, stale_hill.status)
        self.check("primary=HILL" in stale_hill.status, stale_hill.status)

        output = self.set_position(
            "position fen 7k/8/8/8/8/2K5/8/8 w - - 99 1 moves c3d4"
        )
        self.check("error command=position" not in self.joined(output), self.joined(output))
        rule50_hill = self.snapshot()
        self.check("predicates=34" in rule50_hill.status, rule50_hill.status)
        self.check("primary=HILL" in rule50_hill.status, rule50_hill.status)

        output = self.set_position(
            "position startpos moves g1f3 g8f6 f3g1 f6g8 "
            "g1f3 g8f6 f3g1 f6g8"
        )
        self.check("error command=position" not in self.joined(output), self.joined(output))
        repetition = self.snapshot()
        self.check("predicates=8" in repetition.status, repetition.status)
        self.check("repetition_count=3" in repetition.status, repetition.status)
        self.check("primary=AUTOMATIC_DRAW" in repetition.status, repetition.status)

        self.set_position("position startpos moves e2e4")
        raw_ep = self.snapshot()
        self.check(raw_ep.fen.endswith(" b KQkq e3 0 1"), raw_ep.fen)
        self.check("fen_ep=e3 legal_ep=-" in raw_ep.status, raw_ep.status)
        output = self.set_position(f"position fen {raw_ep.fen}")
        self.check("error command=position" not in self.joined(output), self.joined(output))
        self.check(self.snapshot().fen == raw_ep.fen, "raw EP FEN did not round-trip")

        self.set_position(
            "position fen 7k/8/8/8/3p4/8/4P3/K7 w - - 0 1 moves e2e4"
        )
        legal_ep = self.snapshot()
        self.check("fen_ep=e3 legal_ep=e3" in legal_ep.status, legal_ep.status)

        parser_rejections = [
            "position fen 7k/8/8/8/8/2K5/8/8 w - - 0 0",
            "position fen rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w QKkq - 0 1",
            "position fen rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KK - 0 1",
            "position fen 7k/8/8/8/8/2K5/8/8 w K - 0 1",
            "position fen rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w H - 0 1",
            "position fen rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq e3 0 1",
            "position startpos moves E2E4",
            "position startpos moves",
            "position startpos trailing",
        ]
        for command in parser_rejections:
            self.set_position("position startpos")
            before = self.snapshot()
            output = self.set_position(command)
            self.check(
                any("error" in line.lower() for line in output),
                f"parser accepted {command!r}: {output!r}",
            )
            self.check(self.snapshot() == before, f"parser failure mutated state: {command}")

        for terminal_root in (
            "7k/6Q1/6K1/8/8/8/8/8 b - - 0 1",
            "k7/2Q5/8/8/8/8/8/K7 b - - 0 1",
            "7k/8/8/8/8/2K5/8/8 w - - 100 1",
        ):
            self.set_position("position startpos")
            before = self.snapshot()
            output = self.set_position(f"position fen {terminal_root}")
            self.check(
                any("code=TERMINAL_RAW_ROOT" in line for line in output),
                self.joined(output),
            )
            self.check(self.snapshot() == before, "terminal raw root mutated live state")

        terminal_state = self.assert_hill(
            "7k/8/8/8/8/2K5/8/8 w - - 0 1", "c3d4", "white"
        )
        self.check("physical_count=3" in terminal_state.moves, terminal_state.moves)
        terminal_perft = self.engine.transact("go perft 1")
        self.check(
            any("Nodes searched: 0" in line for line in terminal_perft),
            self.joined(terminal_perft),
        )

        self.set_position("position startpos")
        start_perft = self.engine.transact("go perft 2")
        self.check(
            any("Nodes searched: 400" in line for line in start_perft),
            self.joined(start_perft),
        )

        blocked_search = self.engine.transact("go depth 1")
        self.check(
            any("code=KOTH_EVALUATOR_NOT_AUTHENTICATED" in line for line in blocked_search),
            self.joined(blocked_search),
        )
        self.check("bestmove (none)" in blocked_search, self.joined(blocked_search))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("engine", type=Path)
    args = parser.parse_args()
    executable = args.engine.resolve()
    if not executable.is_file():
        raise ConformanceFailure(f"engine not found: {executable}")

    session = EngineSession(executable)
    suite = Suite(session)
    try:
        suite.run()
    finally:
        session.close()

    print(f"PASS koth_rule_core assertions={suite.count}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ConformanceFailure as error:
        print(f"FAIL {error}", file=sys.stderr)
        raise SystemExit(1)
