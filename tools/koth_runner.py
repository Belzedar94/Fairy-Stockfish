#!/usr/bin/env python3
"""Strict two-engine UCI runner for KOTH correctness and local panels."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import queue
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Callable

import koth_referee as referee_module
import koth_book


class RunnerFailure(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


class UciEngine:
    def __init__(self, executable: Path, network: Path, role: str) -> None:
        self.executable = executable
        self.network = network
        self.role = role
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
            raise RunnerFailure(f"{role}: failed to open engine pipes")
        self.stdin = self.process.stdin
        self.stdout = self.process.stdout
        self.lines: queue.Queue[str | None] = queue.Queue()
        self.reader = threading.Thread(target=self._read_output, daemon=True)
        self.reader.start()

        self.send("uci")
        self.uci_output = self.read_until(lambda line: line == "uciok")
        joined = "\n".join(self.uci_output)
        required_variant = (
            "option name UCI_Variant type combo default kingofthehill var kingofthehill"
        )
        if required_variant not in joined:
            raise RunnerFailure(f"{role}: exact KOTH UCI option missing")
        for forbidden in ("SyzygyPath", "UCI_Chess960", "UCI_ShowWDL"):
            if f"option name {forbidden}" in joined:
                raise RunnerFailure(f"{role}: forbidden option exposed: {forbidden}")

        for name, value in (
            ("EvalFile", str(network)),
            ("UCI_Variant", "kingofthehill"),
            ("Threads", "1"),
            ("Hash", "16"),
            ("Ponder", "false"),
            ("MultiPV", "1"),
        ):
            self.send(f"setoption name {name} value {value}")
        loaded = self.ready()
        if not any("network loaded=true" in line for line in loaded):
            raise RunnerFailure(f"{role}: external network did not load: {loaded!r}")
        self.configuration_output = loaded

    def _read_output(self) -> None:
        try:
            for line in self.stdout:
                self.lines.put(line.rstrip("\r\n"))
        finally:
            self.lines.put(None)

    def send(self, command: str) -> None:
        if self.process.poll() is not None:
            raise RunnerFailure(f"{self.role}: engine exited before {command!r}")
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
                raise RunnerFailure(f"{self.role}: output timeout: {output!r}")
            try:
                line = self.lines.get(timeout=remaining)
            except queue.Empty as error:
                raise RunnerFailure(f"{self.role}: output timeout: {output!r}") from error
            if line is None:
                raise RunnerFailure(
                    f"{self.role}: engine exited with code {self.process.poll()}: {output!r}"
                )
            output.append(line)
            if predicate(line):
                return output

    def ready(self) -> list[str]:
        self.send("isready")
        return self.read_until(lambda line: line == "readyok")

    def choose_move(
        self,
        root_fen: str,
        moves: list[str],
        clocks: dict[str, int],
        increment_ms: int,
    ) -> tuple[str, int, list[str]]:
        position = f"position fen {root_fen}"
        if moves:
            position += " moves " + " ".join(moves)
        self.send(position)
        positioned = self.ready()
        if any("error command=position" in line for line in positioned):
            raise RunnerFailure(f"{self.role}: position rejected: {positioned!r}")

        command = (
            f"go wtime {clocks['white']} btime {clocks['black']} "
            f"winc {increment_ms} binc {increment_ms}"
        )
        started = time.monotonic_ns()
        self.send(command)
        output = self.read_until(lambda line: line.startswith("bestmove "))
        elapsed_ms = max(1, math.ceil((time.monotonic_ns() - started) / 1_000_000))
        fields = output[-1].split()
        if len(fields) < 2 or fields[1] == "(none)":
            raise RunnerFailure(f"{self.role}: no move returned: {output!r}")
        if any("code=KOTH_INVALID_PV" in line for line in output):
            raise RunnerFailure(f"{self.role}: invalid PV: {output!r}")
        return fields[1], elapsed_ms, positioned + output

    def identity(self) -> dict[str, object]:
        names = [line for line in self.uci_output if line.startswith("id name ")]
        authors = [line for line in self.uci_output if line.startswith("id author ")]
        return {
            "role": self.role,
            "id_name": names[0][8:] if len(names) == 1 else None,
            "id_author": authors[0][10:] if len(authors) == 1 else None,
            "engine_bytes": self.executable.stat().st_size,
            "engine_sha256": sha256_file(self.executable),
            "network_basename": self.network.name,
            "network_bytes": self.network.stat().st_size,
            "network_sha256": sha256_file(self.network),
            "options": {
                "UCI_Variant": "kingofthehill",
                "Threads": 1,
                "Hash": 16,
                "Ponder": False,
                "MultiPV": 1,
                "EvalFile": self.network.name,
            },
        }

    def close(self) -> None:
        if self.process.poll() is None:
            try:
                self.send("quit")
            except (BrokenPipeError, OSError, RunnerFailure):
                pass
        try:
            code = self.process.wait(timeout=10)
        except subprocess.TimeoutExpired as error:
            self.process.kill()
            self.process.wait(timeout=10)
            raise RunnerFailure(f"{self.role}: engine did not exit") from error
        if code != 0:
            raise RunnerFailure(f"{self.role}: engine exit code {code}")


def run_game(
    white: UciEngine,
    black: UciEngine,
    root_fen: str,
    initial_ms: int,
    increment_ms: int,
    max_plies: int,
    require_terminal: bool,
    book: dict[str, object] | None,
) -> dict[str, object]:
    referee = referee_module.Referee(
        root_fen, referee_module.ClockConfig(initial_ms, increment_ms)
    )
    events: list[dict[str, object]] = []

    while not referee.terminal and len(referee.accepted_moves) < max_plies:
        mover = referee.board.turn
        color = "white" if mover else "black"
        engine = white if mover else black
        before = dict(referee.clocks)
        move, elapsed_ms, output = engine.choose_move(
            root_fen,
            [record["uci"] for record in referee.accepted_moves],
            before,
            increment_ms,
        )
        status = referee.submit(move, referee_module.MoveTiming(elapsed_ms, 0))
        events.append(
            {
                "ply": len(referee.accepted_moves),
                "engine_role": engine.role,
                "mover": color,
                "uci": move,
                "elapsed_ms": elapsed_ms,
                "clocks_before_ms": before,
                "clocks_after_ms": dict(referee.clocks),
                "terminal": status.terminal,
                "primary": status.primary,
                "winner": status.winner,
                "uci_output": output,
            }
        )

    if require_terminal and not referee.terminal:
        raise RunnerFailure(f"game did not terminate within {max_plies} plies")

    pgn = referee.pgn(
        {
            "White": white.role,
            "Black": black.role,
            "KOTHRunner": "KOTH_UCI_RUNNER_V1",
        }
    )
    return {
        "schema": "KOTH_UCI_RUNNER_RECORD_V1",
        "schema_version": 1,
        "evidence_class": "E1_ENGINEERING",
        "strength_claim": False,
        "profile": referee_module.PROFILE,
        "clock_profile": referee_module.CLOCK_PROFILE,
        "root_fen": root_fen,
        "initial_ms": initial_ms,
        "increment_ms": increment_ms,
        "max_plies": max_plies,
        "book": book,
        "engines": {"white": white.identity(), "black": black.identity()},
        "events": events,
        "referee": referee.record(),
        "pgn": pgn,
        "pgn_sha256": hashlib.sha256(pgn.encode("utf-8")).hexdigest().upper(),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--white-engine", type=Path, required=True)
    parser.add_argument("--black-engine", type=Path, required=True)
    parser.add_argument("--white-network", type=Path, required=True)
    parser.add_argument("--black-network", type=Path, required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--root-fen")
    source.add_argument("--book", type=Path)
    parser.add_argument("--book-index", type=int, default=0)
    parser.add_argument("--initial-ms", type=int, default=5_000)
    parser.add_argument("--increment-ms", type=int, default=0)
    parser.add_argument("--max-plies", type=int, default=256)
    parser.add_argument("--require-terminal", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    paths = [
        args.white_engine.resolve(),
        args.black_engine.resolve(),
        args.white_network.resolve(),
        args.black_network.resolve(),
    ]
    if not all(path.is_file() for path in paths):
        raise RunnerFailure(f"missing engine or network input: {paths}")
    if args.initial_ms <= 0 or args.increment_ms < 0 or args.max_plies <= 0:
        raise RunnerFailure("invalid clock or ply limit")

    book_context = None
    if args.book is None:
        root_fen = args.root_fen
        assert root_fen is not None
    else:
        book_path = args.book.resolve()
        if not book_path.is_file():
            raise RunnerFailure(f"book not found: {book_path}")
        book_manifest = koth_book.load_book(book_path)
        records = book_manifest["records"]
        assert isinstance(records, list)
        if args.book_index < 0 or args.book_index >= len(records):
            raise RunnerFailure(f"book index out of range: {args.book_index}")
        entry = records[args.book_index]
        assert isinstance(entry, dict)
        root_fen = str(entry["fen"])
        book_context = {
            "basename": book_path.name,
            "bytes": book_manifest["bytes"],
            "sha256": book_manifest["sha256"],
            "role": book_manifest["role"],
            "strength_book": book_manifest["strength_book"],
            "pairing_contract": book_manifest["pairing_contract"],
            "index": args.book_index,
            "record_id": entry["id"],
        }

    white = UciEngine(paths[0], paths[2], "candidate-white")
    try:
        black = UciEngine(paths[1], paths[3], "comparator-black")
        try:
            result = run_game(
                white,
                black,
                root_fen,
                args.initial_ms,
                args.increment_ms,
                args.max_plies,
                args.require_terminal,
                book_context,
            )
        finally:
            black.close()
    finally:
        white.close()

    rendered = (json.dumps(result, indent=2, sort_keys=True) + "\n").encode("utf-8")
    if args.output is None:
        sys.stdout.buffer.write(rendered)
    else:
        output = args.output.resolve()
        if output.exists():
            raise RunnerFailure(f"output already exists: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(rendered)
        print(
            "PASS koth_runner "
            f"plies={len(result['events'])} result={result['referee']['result']} "
            f"output_sha256={hashlib.sha256(rendered).hexdigest().upper()}"
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (
        OSError,
        subprocess.SubprocessError,
        referee_module.RefereeError,
        RunnerFailure,
        ValueError,
    ) as error:
        print(f"FAIL {error}", file=sys.stderr)
        raise SystemExit(1)
