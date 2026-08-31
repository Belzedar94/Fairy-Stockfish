#!/usr/bin/env python3
"""Contract test for the project-owned KOTH runner canary book."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS))

import koth_book  # noqa: E402


class BookTestFailure(RuntimeError):
    pass


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("book", type=Path)
    args = parser.parse_args()
    book = args.book.resolve()
    if not book.is_file():
        raise BookTestFailure(f"book not found: {book}")
    manifest = koth_book.load_book(book)
    checks = (
        manifest["schema"] == "KOTH_RUNNER_CANARY_BOOK_V1",
        manifest["role"] == "CORRECTNESS_AND_RUNNER_CANARY_ONLY",
        manifest["strength_book"] is False,
        manifest["record_count"] == 16,
        manifest["white_to_move"] == 8,
        manifest["black_to_move"] == 8,
        manifest["goal_records"] == 8,
        manifest["attacked_records"] == 8,
        len({record["id"] for record in manifest["records"]}) == 16,
        len({record["epd"] for record in manifest["records"]}) == 16,
    )
    if not all(checks):
        raise BookTestFailure(str(manifest))
    print(
        "PASS koth_book_contract assertions=10 records=16 "
        f"sha256={manifest['sha256']}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, BookTestFailure, koth_book.BookFailure, ValueError) as error:
        print(f"FAIL {error}", file=sys.stderr)
        raise SystemExit(1)
