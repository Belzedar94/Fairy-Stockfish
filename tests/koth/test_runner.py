#!/usr/bin/env python3
"""One-move white and black canaries for the strict KOTH UCI runner."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path


class RunnerTestFailure(RuntimeError):
    pass


class Suite:
    CASES = (
        ("7k/8/8/8/8/2K5/8/8 w - - 0 1", "white", "1-0"),
        ("8/8/8/2k5/8/8/8/7K b - - 0 1", "black", "0-1"),
    )

    def __init__(self, runner: Path, engine: Path, network: Path) -> None:
        self.runner = runner
        self.engine = engine
        self.network = network
        self.assertions = 0

    def check(self, condition: bool, message: str) -> None:
        self.assertions += 1
        if not condition:
            raise RunnerTestFailure(message)

    def run(self) -> None:
        with tempfile.TemporaryDirectory(prefix="koth-runner-canary-") as directory:
            for index, (fen, winner, result) in enumerate(self.CASES):
                output = Path(directory) / f"runner-{index}.json"
                completed = subprocess.run(
                    [
                        sys.executable,
                        str(self.runner),
                        "--white-engine",
                        str(self.engine),
                        "--black-engine",
                        str(self.engine),
                        "--white-network",
                        str(self.network),
                        "--black-network",
                        str(self.network),
                        "--root-fen",
                        fen,
                        "--initial-ms",
                        "5000",
                        "--max-plies",
                        "1",
                        "--require-terminal",
                        "--output",
                        str(output),
                    ],
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    timeout=30,
                    check=False,
                )
                self.check(completed.returncode == 0, completed.stderr)
                self.check(completed.stdout.startswith("PASS koth_runner "), completed.stdout)
                record = json.loads(output.read_text(encoding="utf-8"))
                self.check(record["schema"] == "KOTH_UCI_RUNNER_RECORD_V1", str(record))
                self.check(record["evidence_class"] == "E1_ENGINEERING", str(record))
                self.check(record["strength_claim"] is False, str(record))
                self.check(record["profile"] == "KOTH_LICHESS_V1", str(record))
                self.check(len(record["events"]) == 1, str(record["events"]))
                event = record["events"][0]
                self.check(event["terminal"] is True, str(event))
                self.check(event["primary"] == "HILL", str(event))
                self.check(event["winner"] == winner, str(event))
                self.check(event["elapsed_ms"] > 0, str(event))
                self.check(
                    event["clocks_after_ms"][event["mover"]]
                    < event["clocks_before_ms"][event["mover"]],
                    str(event),
                )
                self.check(record["referee"]["result"] == result, str(record["referee"]))
                self.check(
                    record["referee"]["board_status"]["winner"] == winner,
                    str(record["referee"]),
                )
                self.check("[KOTHProfile \"KOTH_LICHESS_V1\"]" in record["pgn"], record["pgn"])
                for identity in record["engines"].values():
                    self.check(identity["options"]["UCI_Variant"] == "kingofthehill", str(identity))
                    self.check(identity["options"]["Threads"] == 1, str(identity))
                    self.check(identity["options"]["Hash"] == 16, str(identity))
                    self.check(identity["options"]["Ponder"] is False, str(identity))
                    self.check(identity["network_sha256"] == self._network_sha(), str(identity))

    def _network_sha(self) -> str:
        import hashlib

        digest = hashlib.sha256()
        with self.network.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        return digest.hexdigest().upper()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("runner", type=Path)
    parser.add_argument("engine", type=Path)
    parser.add_argument("network", type=Path)
    args = parser.parse_args()
    paths = [args.runner.resolve(), args.engine.resolve(), args.network.resolve()]
    if not all(path.is_file() for path in paths):
        raise RunnerTestFailure(f"missing input: {paths}")
    suite = Suite(*paths)
    suite.run()
    print(f"PASS koth_runner_canary assertions={suite.assertions} games=2")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, subprocess.SubprocessError, RunnerTestFailure, ValueError) as error:
        print(f"FAIL {error}", file=sys.stderr)
        raise SystemExit(1)
