#!/usr/bin/env python3
"""End-to-end correctness checks for the KOTH T0 search contract.

The suite treats UCI output as the public boundary. It exercises authoritative
terminal classification without an evaluator, legal and attacked hill entries,
terminal precedence in PVs, qsearch quiet wins, root filtering, MultiPV, and the
T0 bans on transposition-table and Syzygy authority.
"""

from __future__ import annotations

import argparse
import queue
import re
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path
from typing import Callable


class SearchFailure(RuntimeError):
    pass


class EngineSession:
    def __init__(
        self, executable: Path, network: Path | None, cwd: Path | None = None
    ) -> None:
        self.process = subprocess.Popen(
            [str(executable)],
            cwd=cwd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            encoding="utf-8",
            errors="replace",
            bufsize=1,
        )
        if self.process.stdin is None or self.process.stdout is None:
            raise SearchFailure("failed to open engine pipes")

        self.stdin = self.process.stdin
        self.stdout = self.process.stdout
        self.lines: queue.Queue[str | None] = queue.Queue()
        self.reader = threading.Thread(target=self._read_output, daemon=True)
        self.reader.start()

        self.send("uci")
        self.uci_output = self.read_until(lambda line: line == "uciok")
        if network is not None:
            self.send(f"setoption name EvalFile value {network}")
            loaded = self.ready()
            if not any("network loaded=true" in line for line in loaded):
                raise SearchFailure(f"legacy network did not load: {loaded!r}")

    def _read_output(self) -> None:
        try:
            for line in self.stdout:
                self.lines.put(line.rstrip("\r\n"))
        finally:
            self.lines.put(None)

    def send(self, command: str) -> None:
        if self.process.poll() is not None:
            raise SearchFailure(
                f"engine exited with code {self.process.returncode} before {command!r}"
            )
        self.stdin.write(command + "\n")
        self.stdin.flush()

    def read_until(
        self, predicate: Callable[[str], bool], timeout: float = 30.0
    ) -> list[str]:
        deadline = time.monotonic() + timeout
        output: list[str] = []
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise SearchFailure(f"engine output timeout; output={output!r}")
            try:
                line = self.lines.get(timeout=remaining)
            except queue.Empty as error:
                raise SearchFailure(f"engine output timeout; output={output!r}") from error
            if line is None:
                raise SearchFailure(
                    f"engine exited with code {self.process.poll()}; output={output!r}"
                )
            output.append(line)
            if predicate(line):
                return output

    def ready(self) -> list[str]:
        self.send("isready")
        return self.read_until(lambda line: line == "readyok")

    def setoption(self, name: str, value: str | int) -> list[str]:
        self.send(f"setoption name {name} value {value}")
        return self.ready()

    def search(
        self,
        position: str,
        go: str,
        *,
        synchronize_position: bool = True,
    ) -> list[str]:
        self.send(position)
        if synchronize_position:
            positioned = self.ready()
            if any("error command=position" in line for line in positioned):
                raise SearchFailure("\n".join(positioned))
        self.send(go)
        return self.read_until(lambda line: line.startswith("bestmove "))

    def close(self) -> None:
        if self.process.poll() is None:
            try:
                self.send("quit")
            except (BrokenPipeError, OSError, SearchFailure):
                pass
        try:
            exit_code = self.process.wait(timeout=10)
        except subprocess.TimeoutExpired as error:
            self.process.kill()
            self.process.wait(timeout=10)
            raise SearchFailure("engine did not exit after quit") from error
        if exit_code != 0:
            raise SearchFailure(f"engine exit code {exit_code}")


class Suite:
    LEGAL_ENTRIES = (
        ("7k/8/8/8/8/2K5/8/8 w - - 0 1", "c3d4", "white"),
        ("k7/8/8/8/8/5K2/8/8 w - - 0 1", "f3e4", "white"),
        ("8/8/2K5/8/8/8/8/7k w - - 0 1", "c6d5", "white"),
        ("8/8/5K2/8/8/8/8/k7 w - - 0 1", "f6e5", "white"),
        ("8/8/8/2k5/8/8/8/7K b - - 0 1", "c5d4", "black"),
        ("8/8/8/5k2/8/8/8/K7 b - - 0 1", "f5e4", "black"),
        ("7K/8/8/8/2k5/8/8/8 b - - 0 1", "c4d5", "black"),
        ("K7/8/8/8/5k2/8/8/8 b - - 0 1", "f4e5", "black"),
    )

    ATTACKED_ENTRIES = (
        ("3r3k/8/8/8/8/2K5/8/8 w - - 0 1", "c3d4"),
        ("k3r3/8/8/8/8/5K2/8/8 w - - 0 1", "f3e4"),
        ("8/8/2K5/8/8/8/8/3r3k w - - 0 1", "c6d5"),
        ("8/8/5K2/8/8/8/8/k3r3 w - - 0 1", "f6e5"),
        ("3R3K/8/8/8/8/2k5/8/8 b - - 0 1", "c3d4"),
        ("K3R3/8/8/8/8/5k2/8/8 b - - 0 1", "f3e4"),
        ("8/8/2k5/8/8/8/8/3R3K b - - 0 1", "c6d5"),
        ("8/8/5k2/8/8/8/8/K3R3 b - - 0 1", "f6e5"),
    )

    def __init__(self, executable: Path, network: Path) -> None:
        self.executable = executable
        self.network = network
        self.assertions = 0
        self.searches = 0

    def check(self, condition: bool, message: str) -> None:
        self.assertions += 1
        if not condition:
            raise SearchFailure(message)

    @staticmethod
    def joined(lines: list[str]) -> str:
        return "\n".join(lines)

    def checked_search(
        self, engine: EngineSession, position: str, go: str
    ) -> list[str]:
        output = engine.search(position, go)
        self.searches += 1
        self.check(
            not any("code=KOTH_INVALID_PV" in line for line in output),
            self.joined(output),
        )
        for line in output:
            if line.startswith("info depth ") and " nodes " in f" {line} ":
                self.check(" hashfull 0 " in f" {line} ", line)
                self.check(" tbhits 0 " in f" {line} ", line)
        return output

    def assert_terminal_pv(
        self,
        output: list[str],
        *,
        predicates: int,
        primary: str,
        winner: str,
        move: str,
        mate: int = 1,
    ) -> None:
        joined = self.joined(output)
        expected_status = (
            "info string koth pv_terminal profile=KOTH_LICHESS_V1 "
            f"terminal=true predicates={predicates} primary={primary} winner={winner}"
        )
        self.check(expected_status in output, joined)
        pv_lines = [line for line in output if line.startswith("info depth ")]
        self.check(bool(pv_lines), joined)
        self.check(any(f" score mate {mate} " in f" {line} " for line in pv_lines), joined)
        self.check(any(re.search(rf"\bpv {re.escape(move)}$", line) for line in pv_lines), joined)
        self.check(output[-1] == f"bestmove {move}", joined)
        self.check(" ponder " not in output[-1], output[-1])

    def test_protocol_surface(self, engine: EngineSession) -> None:
        joined = self.joined(engine.uci_output)
        self.check("option name SyzygyPath" not in joined, joined)
        self.check("option name EvalFile type string" in joined, joined)
        self.check(
            "option name UCI_Variant type combo default kingofthehill var kingofthehill"
            in joined,
            joined,
        )

        output = engine.search("position startpos", "go mate 1")
        self.searches += 1
        self.check(
            "info string error code=UNSUPPORTED_KOTH_GO_LIMIT limit=mate" in output,
            self.joined(output),
        )
        self.check(output[-1] == "bestmove (none)", self.joined(output))

        for command in ("bench", "speedtest"):
            engine.send(command)
            engine.send("isready")
            output = engine.read_until(lambda line: line == "readyok")
            self.check(
                f"info string error code=UNSUPPORTED_KOTH_COMMAND command={command}"
                in output,
                self.joined(output),
            )

    def test_terminal_before_evaluator(self) -> None:
        cases = (
            (
                "position fen 7k/8/8/8/8/2K5/8/8 w - - 0 1 moves c3d4",
                2,
                "HILL",
                "white",
            ),
            (
                "position fen 1rkr4/1p1p4/8/8/8/2K5/8/2R5 w - - 0 1 moves c3d4",
                3,
                "CHECKMATE",
                "white",
            ),
            (
                "position fen k7/2Q5/8/8/8/2K5/8/8 w - - 0 1 moves c3d4",
                6,
                "HILL",
                "white",
            ),
            (
                "position fen 7k/8/8/8/8/2K5/8/8 w - - 99 1 moves c3d4",
                34,
                "HILL",
                "white",
            ),
            (
                "position startpos moves g1f3 g8f6 f3g1 f6g8 "
                "g1f3 g8f6 f3g1 f6g8",
                8,
                "AUTOMATIC_DRAW",
                "none",
            ),
            (
                "position fen 7k/8/8/8/8/8/8/K7 w - - 99 1 moves a1a2",
                32,
                "AUTOMATIC_DRAW",
                "none",
            ),
        )

        with tempfile.TemporaryDirectory(prefix="koth-terminal-no-net-") as directory:
            engine = EngineSession(self.executable, None, Path(directory))
            try:
                for position, predicates, primary, winner in cases:
                    output = engine.search(
                        position, "go depth 64", synchronize_position=False
                    )
                    self.searches += 1
                    joined = self.joined(output)
                    self.check("error command=position" not in joined, joined)
                    self.check("network loaded=true" not in joined, joined)
                    self.check(not any(line.startswith("info depth ") for line in output), joined)
                    self.check(
                        (
                            "info string koth profile=KOTH_LICHESS_V1 terminal=true "
                            f"predicates={predicates} primary={primary} winner={winner}"
                        )
                        in output,
                        joined,
                    )
                    self.check(output[-1] == "bestmove (none)", joined)
                    self.check(" ponder " not in output[-1], output[-1])
            finally:
                engine.close()

    def test_all_goal_entries(self, engine: EngineSession) -> None:
        for fen, move, winner in self.LEGAL_ENTRIES:
            output = self.checked_search(
                engine,
                f"position fen {fen}",
                f"go depth 1 searchmoves {move}",
            )
            self.assert_terminal_pv(
                output,
                predicates=2,
                primary="HILL",
                winner=winner,
                move=move,
            )

        limited = self.checked_search(
            engine,
            "position fen 7k/8/8/8/8/2K5/8/8 w - - 0 1",
            "go nodes 1 searchmoves c3d4",
        )
        self.check(limited[-1] == "bestmove c3d4", self.joined(limited))
        self.check(
            any("primary=HILL winner=white" in line for line in limited),
            self.joined(limited),
        )

    def test_terminal_precedence_in_pv(self, engine: EngineSession) -> None:
        cases = (
            (
                "1rkr4/1p1p4/8/8/8/2K5/8/2R5 w - - 0 1",
                3,
                "CHECKMATE",
            ),
            ("k7/2Q5/8/8/8/2K5/8/8 w - - 0 1", 6, "HILL"),
            ("7k/8/8/8/8/2K5/8/8 w - - 99 1", 34, "HILL"),
        )
        for fen, predicates, primary in cases:
            output = self.checked_search(
                engine, f"position fen {fen}", "go depth 1 searchmoves c3d4"
            )
            self.assert_terminal_pv(
                output,
                predicates=predicates,
                primary=primary,
                winner="white",
                move="c3d4",
            )

    def test_attacked_entries_and_root_filter(self, engine: EngineSession) -> None:
        for fen, attacked_move in self.ATTACKED_ENTRIES:
            unrestricted = self.checked_search(
                engine, f"position fen {fen}", "go depth 1"
            )
            self.check(
                unrestricted[-1] != f"bestmove {attacked_move}",
                self.joined(unrestricted),
            )

            restricted = self.checked_search(
                engine,
                f"position fen {fen}",
                f"go depth 1 searchmoves {attacked_move}",
            )
            self.check(restricted[-1] == "bestmove (none)", self.joined(restricted))
            self.check(
                not any(f" pv {attacked_move}" in line for line in restricted),
                self.joined(restricted),
            )

    def test_multipv_and_qsearch(self, engine: EngineSession) -> None:
        engine.setoption("MultiPV", 2)
        output = self.checked_search(
            engine,
            "position fen 7k/8/8/8/8/3K4/8/8 w - - 0 1",
            "go depth 1 searchmoves d3d4 d3e4",
        )
        pv_lines = [line for line in output if line.startswith("info depth 1 ")]
        self.check(len(pv_lines) == 2, self.joined(output))
        self.check(
            {re.search(r"\bpv (\S+)$", line).group(1) for line in pv_lines}
            == {"d3d4", "d3e4"},
            self.joined(output),
        )
        self.check(all(" score mate 1 " in f" {line} " for line in pv_lines), self.joined(output))
        self.check(
            sum(line.startswith("info string koth pv_terminal ") for line in output) == 2,
            self.joined(output),
        )
        self.check(output[-1] in {"bestmove d3d4", "bestmove d3e4"}, self.joined(output))

        engine.setoption("MultiPV", 1)
        restricted = self.checked_search(
            engine,
            "position fen 7k/8/8/8/8/3K4/8/8 w - - 0 1",
            "go depth 1 searchmoves d3e4",
        )
        self.check(restricted[-1] == "bestmove d3e4", self.joined(restricted))

        qsearch = self.checked_search(
            engine,
            "position fen 7k/8/8/8/8/2K5/8/6R1 w - - 0 1",
            "go depth 1 searchmoves g1g8",
        )
        joined = self.joined(qsearch)
        self.check(" score mate 2 " in f" {joined} ", joined)
        self.check(
            any(
                re.search(r"\bpv g1g8 \S+ c3d4$", line)
                for line in qsearch
                if line.startswith("info depth ")
            ),
            joined,
        )
        self.check(
            any("primary=HILL winner=white" in line for line in qsearch), joined
        )
        self.check(qsearch[-1].startswith("bestmove g1g8"), joined)

    @staticmethod
    def deterministic_signature(output: list[str]) -> tuple[str, str, str, str]:
        final_info = [line for line in output if line.startswith("info depth 3 ")][-1]
        score = re.search(r"\bscore (cp -?\d+|mate -?\d+)", final_info)
        nodes = re.search(r"\bnodes (\d+)", final_info)
        pv = re.search(r"\bpv (.*)$", final_info)
        if score is None or nodes is None or pv is None:
            raise SearchFailure(f"incomplete deterministic signature: {final_info}")
        return output[-1], score.group(1), nodes.group(1), pv.group(1)

    def test_t0_determinism(self, engine: EngineSession) -> None:
        signatures: list[tuple[str, str, str, str]] = []
        for _ in range(3):
            engine.send("ucinewgame")
            engine.ready()
            output = self.checked_search(engine, "position startpos", "go depth 3")
            signatures.append(self.deterministic_signature(output))
        self.check(signatures[0] == signatures[1] == signatures[2], str(signatures))

        engine.setoption("Threads", 2)
        two_thread = self.checked_search(
            engine,
            "position fen 7k/8/8/8/8/2K5/8/8 w - - 0 1",
            "go depth 8",
        )
        self.check(two_thread[-1] == "bestmove c3d4", self.joined(two_thread))
        self.check(
            all(" hashfull 0 " in f" {line} " for line in two_thread if line.startswith("info depth ")),
            self.joined(two_thread),
        )
        engine.setoption("Threads", 1)

    def run(self) -> None:
        self.test_terminal_before_evaluator()
        engine = EngineSession(self.executable, self.network)
        try:
            self.test_protocol_surface(engine)
            self.test_all_goal_entries(engine)
            self.test_terminal_precedence_in_pv(engine)
            self.test_attacked_entries_and_root_filter(engine)
            self.test_multipv_and_qsearch(engine)
            self.test_t0_determinism(engine)
        finally:
            engine.close()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("engine", type=Path)
    parser.add_argument("network", type=Path)
    args = parser.parse_args()
    executable = args.engine.resolve()
    network = args.network.resolve()
    if not executable.is_file():
        raise SearchFailure(f"engine not found: {executable}")
    if not network.is_file():
        raise SearchFailure(f"network not found: {network}")

    suite = Suite(executable, network)
    suite.run()
    print(
        "PASS koth_search_t0 "
        f"assertions={suite.assertions} searches={suite.searches} "
        "goal_entries=8 attacked_entries=8 tt=disabled syzygy=disabled"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, subprocess.SubprocessError, SearchFailure) as error:
        print(f"FAIL {error}", file=sys.stderr)
        raise SystemExit(1)
