#!/usr/bin/env python3
"""Validate and describe the project-owned KOTH runner canary EPD."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

import chess.variant

import koth_referee as referee_module


class BookFailure(RuntimeError):
    pass


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest().upper()


def load_book(path: Path) -> dict[str, object]:
    raw = path.read_bytes()
    if raw.startswith(b"\xef\xbb\xbf"):
        raise BookFailure("UTF-8 BOM is forbidden")
    if b"\r" in raw:
        raise BookFailure("book must use LF framing")
    if not raw.endswith(b"\n"):
        raise BookFailure("book must end with LF")
    try:
        lines = raw.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise BookFailure("book is not UTF-8") from error
    if not lines or any(not line.strip() for line in lines):
        raise BookFailure("book contains no records or a blank record")

    records: list[dict[str, object]] = []
    identifiers: set[str] = set()
    epds: set[str] = set()
    coverage: Counter[tuple[str, str, str]] = Counter()
    goal_names = {chess.square_name(square) for square in referee_module.GOAL_SQUARES}

    for index, line in enumerate(lines):
        try:
            board, operations = chess.variant.KingOfTheHillBoard.from_epd(line)
        except (ValueError, TypeError) as error:
            raise BookFailure(f"line {index + 1}: invalid EPD") from error
        identifier = operations.get("c0")
        kind = operations.get("c1")
        move = operations.get("c2")
        if not all(isinstance(value, str) for value in (identifier, kind, move)):
            raise BookFailure(f"line {index + 1}: c0/c1/c2 strings are required")
        assert isinstance(identifier, str) and isinstance(kind, str) and isinstance(move, str)
        if identifier in identifiers:
            raise BookFailure(f"duplicate record id: {identifier}")
        identifiers.add(identifier)

        epd = board.epd(en_passant="fen")
        if epd in epds:
            raise BookFailure(f"duplicate physical EPD: {epd}")
        epds.add(epd)
        if not board.is_valid():
            raise BookFailure(f"line {index + 1}: invalid physical board")

        root_fen = board.fen(en_passant="fen")
        referee = referee_module.Referee(root_fen)
        if referee.terminal:
            raise BookFailure(f"line {index + 1}: terminal root")
        legal = set(referee.physical_moves())
        color = "white" if board.turn else "black"
        target = move[2:4] if len(move) == 4 else "invalid"
        if target not in goal_names:
            raise BookFailure(f"line {index + 1}: move target is not a goal square")

        if kind == "GOAL":
            if move not in legal:
                raise BookFailure(f"line {index + 1}: expected goal move is not legal")
            result = referee.submit(move)
            if result.primary != "HILL" or result.winner != color:
                raise BookFailure(f"line {index + 1}: goal transition mismatch")
        elif kind == "ATTACKED":
            if move in legal:
                raise BookFailure(f"line {index + 1}: attacked move is legal")
        else:
            raise BookFailure(f"line {index + 1}: unknown class {kind!r}")

        coverage[(kind, color, target)] += 1
        records.append(
            {
                "index": index,
                "id": identifier,
                "kind": kind,
                "epd": epd,
                "fen": root_fen,
                "side_to_move": color,
                "probe_move": move,
                "goal_square": target,
                "legal_move_count": len(legal),
            }
        )

    expected = {
        (kind, color, goal)
        for kind in ("GOAL", "ATTACKED")
        for color in ("white", "black")
        for goal in goal_names
    }
    if set(coverage) != expected or any(count != 1 for count in coverage.values()):
        raise BookFailure(f"coverage mismatch: {dict(coverage)}")

    return {
        "schema": "KOTH_RUNNER_CANARY_BOOK_V1",
        "schema_version": 1,
        "role": "CORRECTNESS_AND_RUNNER_CANARY_ONLY",
        "strength_book": False,
        "source": "PROJECT_OWNED_GOLDEN_RULE_FIXTURES",
        "path_basename": path.name,
        "bytes": len(raw),
        "sha256": sha256_bytes(raw),
        "record_count": len(records),
        "white_to_move": sum(record["side_to_move"] == "white" for record in records),
        "black_to_move": sum(record["side_to_move"] == "black" for record in records),
        "goal_records": sum(record["kind"] == "GOAL" for record in records),
        "attacked_records": sum(record["kind"] == "ATTACKED" for record in records),
        "pairing_contract": "EACH_RECORD_MUST_BE_PLAYED_TWICE_WITH_ENGINE_COLORS_SWAPPED",
        "records": records,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("book", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    book = args.book.resolve()
    if not book.is_file():
        raise BookFailure(f"book not found: {book}")
    manifest = load_book(book)
    rendered = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode("utf-8")
    if args.output is None:
        sys.stdout.buffer.write(rendered)
    else:
        output = args.output.resolve()
        if output.exists():
            raise BookFailure(f"output already exists: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(rendered)
        print(
            "PASS koth_book "
            f"records={manifest['record_count']} sha256={manifest['sha256']} "
            f"output_sha256={sha256_bytes(rendered)}"
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, BookFailure, referee_module.RefereeError, ValueError) as error:
        print(f"FAIL {error}", file=sys.stderr)
        raise SystemExit(1)
