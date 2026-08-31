#!/usr/bin/env python3
"""Independent smoke and integrity checks for a packaged KOTH runtime."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

from test_search import EngineSession, SearchFailure


NETWORK_SHA256 = "978B86D0E6A45E05F9F1375DCED129CEA0ACEA13041EA65960691632EDC47AF7"
EXPECTED_FILES = {
    "AUTHORS",
    "books/koth-runner-canary-v1.epd",
    "Copying.txt",
    "KOTH-RUNTIME.md",
    "KOTH-Stockfish.exe",
    "manifest.json",
    "SHA256SUMS",
    "tools/koth_bench.py",
    "tools/koth_book.py",
    "tools/koth_referee.py",
    "tools/koth-referee-requirements.txt",
    "tools/koth_runner.py",
}


class PackageTestFailure(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


class Suite:
    def __init__(self, package: Path, network: Path) -> None:
        self.package = package
        self.network = network
        self.assertions = 0

    def check(self, condition: bool, message: str) -> None:
        self.assertions += 1
        if not condition:
            raise PackageTestFailure(message)

    def verify_archive(self, destination: Path) -> tuple[Path, dict[str, object]]:
        with zipfile.ZipFile(self.package) as archive:
            infos = archive.infolist()
            names = [info.filename for info in infos]
            self.check(len(names) == len(set(names)), f"duplicate ZIP entries: {names}")
            self.check(
                all(name.startswith("KOTH-Stockfish/") for name in names),
                f"wrong package root: {names}",
            )
            for name in names:
                relative = Path(name).relative_to("KOTH-Stockfish")
                self.check(
                    not relative.is_absolute() and ".." not in relative.parts,
                    f"unsafe ZIP entry: {name}",
                )
            archive.extractall(destination)

        root = destination / "KOTH-Stockfish"
        files = {
            path.relative_to(root).as_posix()
            for path in root.rglob("*")
            if path.is_file()
        }
        self.check(files == EXPECTED_FILES, f"package file set mismatch: {files}")
        self.check(not any(name.lower().endswith(".nnue") for name in files), str(files))

        manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
        self.check(manifest["format"] == "KOTH_RUNTIME_V1", str(manifest))
        self.check(manifest["target"] == "windows-x86-64-sse2", str(manifest))
        self.check(manifest["network"]["included"] is False, str(manifest))
        self.check(manifest["network"]["external"] is True, str(manifest))
        self.check(manifest["network"]["sha256"] == NETWORK_SHA256, str(manifest))
        self.check(manifest["network"]["alias_is_default"] is False, str(manifest))

        checksum_lines = (root / "SHA256SUMS").read_text(encoding="ascii").splitlines()
        observed: dict[str, str] = {}
        for line in checksum_lines:
            match = re.fullmatch(r"([0-9A-F]{64}) \*(.+)", line)
            self.check(match is not None, f"invalid checksum line: {line!r}")
            assert match is not None
            observed[match.group(2)] = match.group(1)
        self.check(observed.keys() == files - {"SHA256SUMS"}, str(observed.keys()))
        for name, expected in observed.items():
            self.check(sha256_file(root / name) == expected, name)
        self.check(
            sha256_file(root / "KOTH-Stockfish.exe") == manifest["engine"]["sha256"],
            str(manifest["engine"]),
        )
        return root, manifest

    def verify_runtime(self, root: Path) -> None:
        executable = root / "KOTH-Stockfish.exe"

        missing = subprocess.run(
            [str(executable)],
            input="uci\nisready\n",
            cwd=root,
            text=True,
            encoding="utf-8",
            errors="replace",
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=20,
            check=False,
        )
        self.check(missing.returncode != 0, missing.stdout)
        self.check("code=KOTH_NET_MISSING" in missing.stdout, missing.stdout)
        self.check("readyok" not in missing.stdout, missing.stdout)

        engine = EngineSession(executable, None, root)
        try:
            uci = "\n".join(engine.uci_output)
            self.check(
                "option name UCI_Variant type combo default kingofthehill var kingofthehill"
                in uci,
                uci,
            )
            self.check("option name SyzygyPath" not in uci, uci)

            terminal = engine.search(
                "position fen 7k/8/8/8/8/2K5/8/8 w - - 0 1 moves c3d4",
                "go depth 64",
                synchronize_position=False,
            )
            self.check(
                any("primary=HILL winner=white" in line for line in terminal),
                "\n".join(terminal),
            )
            self.check(terminal[-1] == "bestmove (none)", "\n".join(terminal))

            loaded = engine.setoption("EvalFile", self.network)
            self.check(any("network loaded=true" in line for line in loaded), str(loaded))
            searched = engine.search(
                "position fen 7k/8/8/8/8/2K5/8/8 w - - 0 1",
                "go depth 1 searchmoves c3d4",
            )
            self.check(searched[-1] == "bestmove c3d4", "\n".join(searched))
            self.check(
                any("primary=HILL winner=white" in line for line in searched),
                "\n".join(searched),
            )
        finally:
            engine.close()

        bench = subprocess.run(
            [
                sys.executable,
                str(root / "tools" / "koth_bench.py"),
                str(executable),
                str(self.network),
                "--runs",
                "2",
            ],
            cwd=root,
            text=True,
            encoding="utf-8",
            errors="replace",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=60,
            check=False,
        )
        self.check(bench.returncode == 0, bench.stderr)
        bench_record = json.loads(bench.stdout)
        self.check(bench_record["all_runs_equal"] is True, bench.stdout)
        self.check(
            bench_record["semantic_digest"]
            == "9A66D1BAE80FF3B81B3254D18A1E82BC0202BDA49BEE080307136D1A976280CC",
            bench.stdout,
        )

        runner_output = root / "runner-smoke.json"
        runner = subprocess.run(
            [
                sys.executable,
                str(root / "tools" / "koth_runner.py"),
                "--white-engine",
                str(executable),
                "--black-engine",
                str(executable),
                "--white-network",
                str(self.network),
                "--black-network",
                str(self.network),
                "--book",
                str(root / "books" / "koth-runner-canary-v1.epd"),
                "--book-index",
                "0",
                "--initial-ms",
                "5000",
                "--max-plies",
                "1",
                "--require-terminal",
                "--output",
                str(runner_output),
            ],
            cwd=root,
            text=True,
            encoding="utf-8",
            errors="replace",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=30,
            check=False,
        )
        self.check(runner.returncode == 0, runner.stderr)
        runner_record = json.loads(runner_output.read_text(encoding="utf-8"))
        self.check(runner_record["strength_claim"] is False, str(runner_record))
        self.check(runner_record["book"]["strength_book"] is False, str(runner_record))
        self.check(
            runner_record["book"]["record_id"] == "KOTH-CANARY-GOAL-W-D4",
            str(runner_record),
        )
        self.check(runner_record["referee"]["result"] == "1-0", str(runner_record))
        self.check(
            runner_record["referee"]["board_status"]["primary"] == "HILL",
            str(runner_record),
        )

        payload = {
            "root_fen": "7k/8/8/8/8/2K5/8/8 w - - 0 1",
            "moves": [{"uci": "c3d4"}],
        }
        referee = subprocess.run(
            [sys.executable, str(root / "tools" / "koth_referee.py")],
            input=json.dumps(payload),
            text=True,
            encoding="utf-8",
            errors="replace",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=20,
            check=False,
        )
        self.check(referee.returncode == 0, referee.stderr)
        record = json.loads(referee.stdout)["record"]
        self.check(record["profile"] == "KOTH_LICHESS_V1", str(record))
        self.check(record["board_status"]["primary"] == "HILL", str(record))
        self.check(record["board_status"]["winner"] == "white", str(record))
        self.check(record["game_moves"] == [], str(record))
        self.check(bool(record["physical_moves"]), str(record))

    def run(self) -> None:
        with tempfile.TemporaryDirectory(prefix="koth-package-smoke-") as directory:
            root, _manifest = self.verify_archive(Path(directory))
            self.verify_runtime(root)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("package", type=Path)
    parser.add_argument("network", type=Path)
    args = parser.parse_args()
    package = args.package.resolve()
    network = args.network.resolve()
    if not package.is_file():
        raise PackageTestFailure(f"package not found: {package}")
    if not network.is_file():
        raise PackageTestFailure(f"network not found: {network}")
    if sha256_file(network) != NETWORK_SHA256:
        raise PackageTestFailure("external network SHA-256 mismatch")

    suite = Suite(package, network)
    suite.run()
    print(f"PASS koth_package assertions={suite.assertions} network_included=false")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (
        OSError,
        subprocess.SubprocessError,
        SearchFailure,
        PackageTestFailure,
        KeyError,
        ValueError,
        zipfile.BadZipFile,
    ) as error:
        print(f"FAIL {error}", file=sys.stderr)
        raise SystemExit(1)
