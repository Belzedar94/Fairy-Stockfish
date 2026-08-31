#!/usr/bin/env python3
"""Deterministic fixed-work bench for the KOTH search contract.

The digest excludes elapsed time and NPS. It authenticates the executable and
external legacy network separately, then hashes normalized semantic search
records produced with one thread, a fixed hash size, TT authority disabled, and
Syzygy unavailable.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import queue
import re
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable


NETWORK_BYTES = 47_721_371
NETWORK_SHA256 = "978B86D0E6A45E05F9F1375DCED129CEA0ACEA13041EA65960691632EDC47AF7"
NETWORK_NAMES = {"kingofthehill-978b86d0e6a4.nnue", "KOTH_v1.nnue"}


class BenchFailure(RuntimeError):
    pass


@dataclass(frozen=True)
class Case:
    case_id: str
    position: str
    go: str
    expected_bestmove: str | None = None
    expected_primary: str | None = None
    expected_winner: str | None = None


CASES = (
    Case("startpos-depth4", "position startpos", "go depth 4"),
    Case(
        "middlegame-nodes2048",
        "position fen r3k2r/p1ppqpb1/bn2pnp1/2pP4/1p2P3/2N2N2/PPQ1BPPP/R1B1K2R w KQkq - 0 10",
        "go nodes 2048",
    ),
    Case(
        "endgame-nodes2048",
        "position fen 8/2p5/3p4/1P1Pp2k/4Pp2/5P1K/8/8 w - - 0 40",
        "go nodes 2048",
    ),
    Case(
        "qsearch-quiet-hill",
        "position fen 7k/8/8/8/8/2K5/8/6R1 w - - 0 1",
        "go depth 1 searchmoves g1g8",
        "g1g8",
        "HILL",
        "white",
    ),
    Case(
        "goal-d4-white",
        "position fen 7k/8/8/8/8/2K5/8/8 w - - 0 1",
        "go depth 1 searchmoves c3d4",
        "c3d4",
        "HILL",
        "white",
    ),
    Case(
        "goal-e4-white",
        "position fen k7/8/8/8/8/5K2/8/8 w - - 0 1",
        "go depth 1 searchmoves f3e4",
        "f3e4",
        "HILL",
        "white",
    ),
    Case(
        "goal-d5-white",
        "position fen 8/8/2K5/8/8/8/8/7k w - - 0 1",
        "go depth 1 searchmoves c6d5",
        "c6d5",
        "HILL",
        "white",
    ),
    Case(
        "goal-e5-white",
        "position fen 8/8/5K2/8/8/8/8/k7 w - - 0 1",
        "go depth 1 searchmoves f6e5",
        "f6e5",
        "HILL",
        "white",
    ),
    Case(
        "goal-d4-black",
        "position fen 8/8/8/2k5/8/8/8/7K b - - 0 1",
        "go depth 1 searchmoves c5d4",
        "c5d4",
        "HILL",
        "black",
    ),
    Case(
        "goal-e4-black",
        "position fen 8/8/8/5k2/8/8/8/K7 b - - 0 1",
        "go depth 1 searchmoves f5e4",
        "f5e4",
        "HILL",
        "black",
    ),
    Case(
        "goal-d5-black",
        "position fen 7K/8/8/8/2k5/8/8/8 b - - 0 1",
        "go depth 1 searchmoves c4d5",
        "c4d5",
        "HILL",
        "black",
    ),
    Case(
        "goal-e5-black",
        "position fen K7/8/8/8/5k2/8/8/8 b - - 0 1",
        "go depth 1 searchmoves f4e5",
        "f4e5",
        "HILL",
        "black",
    ),
    Case(
        "precedence-mate-hill",
        "position fen 1rkr4/1p1p4/8/8/8/2K5/8/2R5 w - - 0 1",
        "go depth 1 searchmoves c3d4",
        "c3d4",
        "CHECKMATE",
        "white",
    ),
    Case(
        "precedence-stalemate-hill",
        "position fen k7/2Q5/8/8/8/2K5/8/8 w - - 0 1",
        "go depth 1 searchmoves c3d4",
        "c3d4",
        "HILL",
        "white",
    ),
    Case(
        "precedence-rule50-hill",
        "position fen 7k/8/8/8/8/2K5/8/8 w - - 99 1",
        "go depth 1 searchmoves c3d4",
        "c3d4",
        "HILL",
        "white",
    ),
)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def canonical_bytes(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


class Engine:
    def __init__(self, executable: Path, network: Path) -> None:
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
            raise BenchFailure("failed to open engine pipes")
        self.stdin = self.process.stdin
        self.stdout = self.process.stdout
        self.lines: queue.Queue[str | None] = queue.Queue()
        self.reader = threading.Thread(target=self._read_output, daemon=True)
        self.reader.start()

        self.send("uci")
        self.uci_output = self.read_until(lambda line: line == "uciok")
        for name, value in (
            ("EvalFile", str(network)),
            ("Threads", "1"),
            ("Hash", "16"),
            ("MultiPV", "1"),
        ):
            self.send(f"setoption name {name} value {value}")
        ready = self.ready()
        if not any("network loaded=true" in line for line in ready):
            raise BenchFailure(f"network failed to load: {ready!r}")

    def _read_output(self) -> None:
        try:
            for line in self.stdout:
                self.lines.put(line.rstrip("\r\n"))
        finally:
            self.lines.put(None)

    def send(self, command: str) -> None:
        if self.process.poll() is not None:
            raise BenchFailure(f"engine exited before command {command!r}")
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
                raise BenchFailure(f"engine output timeout: {output!r}")
            try:
                line = self.lines.get(timeout=remaining)
            except queue.Empty as error:
                raise BenchFailure(f"engine output timeout: {output!r}") from error
            if line is None:
                raise BenchFailure(
                    f"engine exited with code {self.process.poll()}: {output!r}"
                )
            output.append(line)
            if predicate(line):
                return output

    def ready(self) -> list[str]:
        self.send("isready")
        return self.read_until(lambda line: line == "readyok")

    def run_case(self, case: Case) -> list[str]:
        self.send("ucinewgame")
        self.ready()
        self.send(case.position)
        positioned = self.ready()
        if any("error command=position" in line for line in positioned):
            raise BenchFailure(f"{case.case_id}: {positioned!r}")
        self.send(case.go)
        return self.read_until(lambda line: line.startswith("bestmove "))

    def close(self) -> None:
        if self.process.poll() is None:
            self.send("quit")
        try:
            code = self.process.wait(timeout=10)
        except subprocess.TimeoutExpired as error:
            self.process.kill()
            self.process.wait(timeout=10)
            raise BenchFailure("engine did not exit after quit") from error
        if code != 0:
            raise BenchFailure(f"engine exit code {code}")


def field(pattern: str, line: str, name: str) -> str:
    match = re.search(pattern, line)
    if match is None:
        raise BenchFailure(f"missing {name} in {line!r}")
    return match.group(1)


def normalize(case: Case, output: list[str]) -> dict[str, object]:
    if any("code=KOTH_INVALID_PV" in line for line in output):
        raise BenchFailure(f"{case.case_id}: invalid PV: {output!r}")
    info_lines = [
        line
        for line in output
        if line.startswith("info depth ") and " nodes " in f" {line} "
    ]
    if not info_lines:
        raise BenchFailure(f"{case.case_id}: no complete info line: {output!r}")
    final = info_lines[-1]
    bestmove = field(r"^bestmove (\S+)", output[-1], "bestmove")
    ponder_match = re.search(r"\bponder (\S+)", output[-1])
    terminal_lines = [
        line.removeprefix("info string koth pv_terminal ")
        for line in output
        if line.startswith("info string koth pv_terminal ")
    ]
    terminal = terminal_lines[-1] if terminal_lines else None

    if case.expected_bestmove is not None and bestmove != case.expected_bestmove:
        raise BenchFailure(
            f"{case.case_id}: bestmove {bestmove} != {case.expected_bestmove}"
        )
    if case.expected_primary is not None:
        if terminal is None or f"primary={case.expected_primary}" not in terminal:
            raise BenchFailure(f"{case.case_id}: wrong terminal status {terminal!r}")
    if case.expected_winner is not None:
        if terminal is None or f"winner={case.expected_winner}" not in terminal:
            raise BenchFailure(f"{case.case_id}: wrong terminal winner {terminal!r}")

    hashfull = int(field(r"\bhashfull (\d+)", final, "hashfull"))
    tbhits = int(field(r"\btbhits (\d+)", final, "tbhits"))
    if hashfull != 0 or tbhits != 0:
        raise BenchFailure(
            f"{case.case_id}: forbidden authority hashfull={hashfull} tbhits={tbhits}"
        )

    return {
        "id": case.case_id,
        "bestmove": bestmove,
        "ponder": ponder_match.group(1) if ponder_match else None,
        "depth": int(field(r"\bdepth (\d+)", final, "depth")),
        "seldepth": int(field(r"\bseldepth (\d+)", final, "seldepth")),
        "score": field(r"\bscore ((?:cp|mate) -?\d+)", final, "score"),
        "nodes": int(field(r"\bnodes (\d+)", final, "nodes")),
        "hashfull": hashfull,
        "tbhits": tbhits,
        "pv": field(r"\bpv (.*)$", final, "pv"),
        "terminal": terminal,
    }


def execute_run(executable: Path, network: Path) -> list[dict[str, object]]:
    engine = Engine(executable, network)
    try:
        return [normalize(case, engine.run_case(case)) for case in CASES]
    finally:
        engine.close()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("engine", type=Path)
    parser.add_argument("network", type=Path)
    parser.add_argument("--runs", type=int, default=3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    executable = args.engine.resolve()
    network = args.network.resolve()
    if not executable.is_file():
        raise BenchFailure(f"engine not found: {executable}")
    if not network.is_file():
        raise BenchFailure(f"network not found: {network}")
    if args.runs < 2:
        raise BenchFailure("--runs must be at least 2")
    if network.name not in NETWORK_NAMES:
        raise BenchFailure(f"network basename is not allowed: {network.name}")
    if network.stat().st_size != NETWORK_BYTES:
        raise BenchFailure(f"network size mismatch: {network.stat().st_size}")
    network_sha = file_sha256(network)
    if network_sha != NETWORK_SHA256:
        raise BenchFailure(f"network SHA-256 mismatch: {network_sha}")

    runs = [execute_run(executable, network) for _ in range(args.runs)]
    digests = [hashlib.sha256(canonical_bytes(run)).hexdigest().upper() for run in runs]
    if len(set(digests)) != 1 or any(run != runs[0] for run in runs[1:]):
        raise BenchFailure(f"nondeterministic semantic results: {digests}")

    manifest = {
        "schema_version": 1,
        "contract": "KOTH_SEARCH_T0_BENCH_V1",
        "engine": {
            "bytes": executable.stat().st_size,
            "sha256": file_sha256(executable),
        },
        "network": {
            "basename": network.name,
            "bytes": NETWORK_BYTES,
            "sha256": network_sha,
            "external": True,
        },
        "settings": {
            "threads": 1,
            "hash_mib": 16,
            "multipv": 1,
            "tt_authority": "disabled",
            "syzygy": "disabled",
            "timing_fields_in_digest": False,
        },
        "case_count": len(CASES),
        "runs": args.runs,
        "run_digests": digests,
        "all_runs_equal": True,
        "semantic_digest": digests[0],
        "cases": runs[0],
    }
    rendered = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8")
    if args.output is not None:
        output = args.output.resolve()
        if output.exists():
            raise BenchFailure(f"output already exists: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(rendered)
        print(
            "PASS koth_bench "
            f"cases={len(CASES)} runs={args.runs} digest={digests[0]} "
            f"output_sha256={hashlib.sha256(rendered).hexdigest().upper()}"
        )
    else:
        sys.stdout.buffer.write(rendered)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, subprocess.SubprocessError, BenchFailure, ValueError) as error:
        print(f"FAIL {error}", file=sys.stderr)
        raise SystemExit(1)
