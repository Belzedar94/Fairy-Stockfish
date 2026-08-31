#!/usr/bin/env python3
"""Exact legacy-network loader and scalar-evaluator conformance tests."""

from __future__ import annotations

import argparse
import hashlib
import os
import random
import re
import shutil
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

import chess

from reference_legacy_nnue import (
    FIRST_LAYER_OFFSET,
    LAYER_BYTES,
    LegacyNetwork,
    SHA256,
)


CANONICAL_NAME = "kingofthehill-978b86d0e6a4.nnue"
ALIAS_NAME = "KOTH_v1.nnue"
GOALS = {chess.D4, chess.E4, chess.D5, chess.E5}
RAW_PATTERN = re.compile(
    r"terminal=false bucket=(?P<bucket>\d+) "
    r"lazy=(?:true|false) "
    r"psqt_raw=(?P<psqt_raw>-?\d+) positional_raw=(?P<positional_raw>-?\d+) "
    r"psqt=(?P<psqt>-?\d+) positional=(?P<positional>-?\d+) "
    r"total=(?P<total>-?\d+)"
)


class TestFailure(RuntimeError):
    pass


def check(condition: bool, message: str) -> None:
    if not condition:
        raise TestFailure(message)


def link_or_copy(source: Path, destination: Path) -> None:
    try:
        os.link(source, destination)
    except OSError:
        shutil.copyfile(source, destination)


def run_engine(
    executable: Path, commands: list[str], cwd: Path | None = None
) -> tuple[int, str]:
    completed = subprocess.run(
        [str(executable)],
        input="\n".join(commands) + "\n",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        encoding="utf-8",
        errors="replace",
        cwd=cwd,
        timeout=30,
        check=False,
    )
    return completed.returncode, completed.stdout


class EngineSession:
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
            raise TestFailure("could not open engine pipes")
        self.stdin = self.process.stdin
        self.stdout = self.process.stdout
        self.transact(["uci"], "uciok")
        loaded = self.transact(
            [f"setoption name EvalFile value {network}", "isready"], "readyok"
        )
        check("network loaded=true" in "\n".join(loaded), "network did not load")

    def transact(self, commands: list[str], terminal: str) -> list[str]:
        for command in commands:
            self.stdin.write(command + "\n")
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

    def evaluate(self, fen: str) -> dict[str, int]:
        lines = self.transact(
            [f"position fen {fen}", "kothneteval", "isready"], "readyok"
        )
        check(
            not any("error command=position" in line for line in lines),
            f"engine rejected reference FEN {fen}: {lines!r}",
        )
        matches = [RAW_PATTERN.search(line) for line in lines]
        matches = [match for match in matches if match is not None]
        check(len(matches) == 1, f"missing/duplicate raw evaluation: {lines!r}")
        return {key: int(value) for key, value in matches[0].groupdict().items()}

    def close(self) -> None:
        if self.process.poll() is None:
            self.stdin.write("quit\n")
            self.stdin.flush()
        try:
            code = self.process.wait(timeout=10)
        except subprocess.TimeoutExpired as error:
            self.process.kill()
            raise TestFailure("engine did not stop") from error
        check(code == 0, f"positive engine exited with code {code}")


class HistoricalEngineSession:
    """Diagnostic lane built from official Stockfish e8d64af1."""

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
            raise TestFailure("could not open historical engine pipes")
        self.stdin = self.process.stdin
        self.stdout = self.process.stdout
        self.transact(["uci"], "uciok")
        self.transact(
            [f"setoption name EvalFile value {network}", "isready"], "readyok"
        )

    def transact(self, commands: list[str], terminal: str) -> list[str]:
        for command in commands:
            self.stdin.write(command + "\n")
        self.stdin.flush()
        lines: list[str] = []
        while True:
            line = self.stdout.readline()
            if not line:
                raise TestFailure(
                    f"historical engine exited before {terminal}: {lines!r}"
                )
            line = line.rstrip("\r\n")
            lines.append(line)
            if line == terminal:
                return lines

    def evaluate(self, fen: str) -> int:
        lines = self.transact(
            [f"position fen {fen}", "nnueraw", "isready"], "readyok"
        )
        raw = [line for line in lines if line.startswith("nnueraw ")]
        check(len(raw) == 1, f"missing historical raw value for {fen}: {lines!r}")
        return int(raw[0].split()[1])

    def close(self) -> None:
        if self.process.poll() is None:
            self.stdin.write("quit\n")
            self.stdin.flush()
        try:
            code = self.process.wait(timeout=10)
        except subprocess.TimeoutExpired as error:
            self.process.kill()
            raise TestFailure("historical engine did not stop") from error
        check(code == 0, f"historical engine exited with code {code}")


def synthetic_fen(rng: random.Random, piece_count: int) -> str:
    eligible = [square for square in chess.SQUARES if square not in {chess.A1, chess.H8}]
    material: list[chess.Piece] = []
    for color in chess.COLORS:
        material.extend([chess.Piece(chess.PAWN, color)] * 8)
        material.extend([chess.Piece(chess.KNIGHT, color)] * 2)
        material.extend([chess.Piece(chess.BISHOP, color)] * 2)
        material.extend([chess.Piece(chess.ROOK, color)] * 2)
        material.append(chess.Piece(chess.QUEEN, color))
    for _ in range(20_000):
        board = chess.Board(None)
        board.set_piece_at(chess.A1, chess.Piece(chess.KING, chess.WHITE))
        board.set_piece_at(chess.H8, chess.Piece(chess.KING, chess.BLACK))
        shuffled_material = list(material)
        rng.shuffle(shuffled_material)
        squares = rng.sample(eligible, piece_count - 2)
        valid_placement = True
        for square, piece in zip(squares, shuffled_material):
            if piece.piece_type == chess.PAWN and chess.square_rank(square) in {0, 7}:
                valid_placement = False
                break
            board.set_piece_at(square, piece)
        if not valid_placement:
            continue
        board.turn = rng.choice(chess.COLORS)
        board.castling_rights = chess.BB_EMPTY
        board.ep_square = None
        board.halfmove_clock = 0
        board.fullmove_number = 1
        if board.is_valid() and not board.is_checkmate() and not board.is_stalemate():
            return board.fen(en_passant="fen")
    raise TestFailure(f"could not construct valid {piece_count}-piece position")


def reference_positions() -> list[str]:
    rng = random.Random(0x4B4F54484E4E5545)
    positions = [
        chess.STARTING_FEN,
        "8/7k/8/8/8/4K3/8/8 w - - 0 1",
        "8/7k/8/8/8/4K3/8/8 b - - 0 1",
    ]
    for bucket in range(8):
        lower = 2 if bucket == 0 else bucket * 4 + 1
        upper = min(32, bucket * 4 + 4)
        for piece_count in range(lower, upper + 1):
            positions.append(synthetic_fen(rng, piece_count))
            positions.append(synthetic_fen(rng, piece_count))

    board = chess.Board()
    for _ in range(48):
        if board.is_game_over(claim_draw=False):
            break
        board.push(rng.choice(list(board.legal_moves)))
        if (
            len(board.move_stack) % 4 == 0
            and board.king(chess.WHITE) not in GOALS
            and board.king(chess.BLACK) not in GOALS
            and not board.is_game_over(claim_draw=False)
        ):
            positions.append(board.fen(en_passant="fen"))
    return positions


def assert_fail_closed(
    executable: Path, path: Path, expected_code: str, assertions: list[int]
) -> None:
    code, output = run_engine(
        executable,
        ["uci", f"setoption name EvalFile value {path}", "isready", "quit"],
    )
    assertions[0] += 4
    check(code != 0, f"invalid net did not abort ({expected_code}): {output}")
    check(expected_code in output, f"expected {expected_code}: {output}")
    check("network loaded=true" not in output, f"invalid net reported loaded: {output}")
    check("readyok" not in output, f"invalid net reached readyok: {output}")


def mutate_copy(source: Path, root: Path, case: str, offset: int) -> Path:
    case_dir = root / case
    case_dir.mkdir()
    target = case_dir / CANONICAL_NAME
    shutil.copyfile(source, target)
    with target.open("r+b") as stream:
        stream.seek(offset)
        original = stream.read(1)
        check(len(original) == 1, f"mutation offset outside network: {offset}")
        stream.seek(offset)
        stream.write(bytes([original[0] ^ 0x01]))
    return target


def run_suite(
    executable: Path, network_path: Path, historical_engine: Path | None
) -> tuple[int, int, int]:
    assertions = [0]
    digest = hashlib.sha256(network_path.read_bytes()).hexdigest().upper()
    check(network_path.name == CANONICAL_NAME, "legacy input must retain canonical basename")
    check(digest == SHA256, f"legacy input digest mismatch: {digest}")
    assertions[0] += 2

    positions = reference_positions()
    observed_buckets: set[int] = set()
    with LegacyNetwork(network_path) as reference:
        engine = EngineSession(executable, network_path)
        try:
            for fen in positions:
                expected = reference.evaluate(fen)
                observed = engine.evaluate(fen)
                expected_fields = {
                    "bucket": expected.bucket,
                    "psqt_raw": expected.psqt_raw,
                    "positional_raw": expected.positional_raw,
                    "psqt": expected.psqt,
                    "positional": expected.positional,
                    "total": expected.total,
                }
                assertions[0] += 1
                check(observed == expected_fields, f"scalar parity mismatch {fen}: {observed} != {expected_fields}")
                observed_buckets.add(observed["bucket"])
        finally:
            engine.close()
    assertions[0] += 1
    check(observed_buckets == set(range(8)), f"bucket coverage incomplete: {observed_buckets}")

    historical_positions = 0
    if historical_engine is not None:
        with LegacyNetwork(network_path) as reference:
            historical = HistoricalEngineSession(historical_engine, network_path)
            try:
                for fen in positions:
                    expected = reference.evaluate(fen).total
                    observed = historical.evaluate(fen)
                    assertions[0] += 1
                    check(
                        observed == expected,
                        f"historical parity mismatch {fen}: {observed} != {expected}",
                    )
                    historical_positions += 1
            finally:
                historical.close()

    with tempfile.TemporaryDirectory(prefix="koth-legacy-nnue-") as temporary:
        root = Path(temporary)
        alias = root / ALIAS_NAME
        link_or_copy(network_path, alias)
        assertions[0] += 1
        check(
            hashlib.sha256(alias.read_bytes()).hexdigest().upper() == SHA256,
            "alias is not byte-identical",
        )
        code, output = run_engine(
            executable,
            ["uci", f"setoption name EvalFile value {alias}", "kothnetstatus", "isready", "quit"],
        )
        assertions[0] += 3
        check(code == 0, f"exact alias failed: {output}")
        check("basename=KOTH_v1.nnue" in output, output)
        check("network loaded=true" in output and "readyok" in output, output)

        default_dir = root / "default"
        default_dir.mkdir()
        link_or_copy(network_path, default_dir / CANONICAL_NAME)
        code, output = run_engine(executable, ["uci", "isready", "quit"], cwd=default_dir)
        assertions[0] += 2
        check(code == 0, f"canonical default failed: {output}")
        check("network loaded=true" in output and "readyok" in output, output)

        wrong_name = root / "wrong.nnue"
        link_or_copy(network_path, wrong_name)
        assert_fail_closed(executable, wrong_name, "code=KOTH_NET_BASENAME", assertions)
        assert_fail_closed(
            executable,
            root / "missing" / CANONICAL_NAME,
            "code=KOTH_NET_MISSING",
            assertions,
        )

        short_dir = root / "short"
        short_dir.mkdir()
        short_path = short_dir / CANONICAL_NAME
        shutil.copyfile(network_path, short_path)
        with short_path.open("r+b") as stream:
            stream.truncate(short_path.stat().st_size - 1)
        assert_fail_closed(executable, short_path, "code=KOTH_NET_SIZE", assertions)

        extra_dir = root / "extra"
        extra_dir.mkdir()
        extra_path = extra_dir / CANONICAL_NAME
        shutil.copyfile(network_path, extra_path)
        with extra_path.open("ab") as stream:
            stream.write(b"\x00")
        assert_fail_closed(executable, extra_path, "code=KOTH_NET_SIZE", assertions)

        fixed_mutations = [
            ("version", 0, "code=KOTH_NET_VERSION"),
            ("architecture", 4, "code=KOTH_NET_ARCHITECTURE"),
            ("description_size", 8, "code=KOTH_NET_DESCRIPTION_SIZE"),
            ("description", 12, "code=KOTH_NET_DESCRIPTION_MISMATCH"),
            ("feature_hash", 87, "code=KOTH_NET_FEATURE_HASH"),
            ("payload", 1_000, "code=KOTH_NET_SHA256"),
        ]
        for case, offset, expected_code in fixed_mutations:
            target = mutate_copy(network_path, root, case, offset)
            assert_fail_closed(executable, target, expected_code, assertions)

        for bucket in range(8):
            target = mutate_copy(
                network_path,
                root,
                f"layer_hash_{bucket}",
                FIRST_LAYER_OFFSET + bucket * LAYER_BYTES,
            )
            assert_fail_closed(
                executable,
                target,
                f"code=KOTH_NET_LAYER_HASH bucket={bucket}",
                assertions,
            )

    return assertions[0], len(positions), historical_positions


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("engine", type=Path)
    parser.add_argument("network", type=Path)
    parser.add_argument("historical_engine", type=Path, nargs="?")
    args = parser.parse_args()
    executable = args.engine.resolve()
    network = args.network.resolve()
    if not executable.is_file():
        raise TestFailure(f"engine not found: {executable}")
    if not network.is_file():
        raise TestFailure(f"network not found: {network}")
    historical_engine = (
        args.historical_engine.resolve() if args.historical_engine is not None else None
    )
    if historical_engine is not None and not historical_engine.is_file():
        raise TestFailure(f"historical engine not found: {historical_engine}")

    assertions, positions, historical_positions = run_suite(
        executable, network, historical_engine
    )
    print(
        "PASS koth_legacy_nnue "
        f"assertions={assertions} scalar_positions={positions} "
        f"historical_positions={historical_positions} buckets=8 negative_cases=18"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, subprocess.SubprocessError, TestFailure, ValueError) as error:
        print(f"FAIL {error}", file=sys.stderr)
        raise SystemExit(1)
