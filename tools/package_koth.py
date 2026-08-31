#!/usr/bin/env python3
"""Create a deterministic network-free KOTH Windows runtime package."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import stat
import sys
import zipfile
from pathlib import Path


NETWORK_SHA256 = "978B86D0E6A45E05F9F1375DCED129CEA0ACEA13041EA65960691632EDC47AF7"
NETWORK_BYTES = 47_721_371
PACKAGE_ROOT = "KOTH-Stockfish"
PAYLOAD_SOURCES = {
    "AUTHORS": "AUTHORS",
    "Copying.txt": "Copying.txt",
    "KOTH-RUNTIME.md": "KOTH-RUNTIME.md",
    "tools/koth_bench.py": "tools/koth_bench.py",
    "tools/koth_referee.py": "tools/koth_referee.py",
    "tools/koth-referee-requirements.txt": "tools/koth-referee-requirements.txt",
}


class PackageFailure(RuntimeError):
    pass


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest().upper()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def zip_info(name: str, epoch: int, executable: bool = False) -> zipfile.ZipInfo:
    timestamp = dt.datetime.fromtimestamp(epoch, tz=dt.timezone.utc)
    if timestamp.year < 1980 or timestamp.year > 2107:
        raise PackageFailure("source date epoch is outside the ZIP timestamp range")
    info = zipfile.ZipInfo(name, timestamp.timetuple()[:6])
    info.compress_type = zipfile.ZIP_DEFLATED
    info.create_system = 3
    mode = 0o755 if executable else 0o644
    info.external_attr = (stat.S_IFREG | mode) << 16
    return info


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--engine", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--engine-source-commit", required=True)
    parser.add_argument("--engine-source-tree", required=True)
    parser.add_argument("--package-source-commit", required=True)
    parser.add_argument("--package-source-tree", required=True)
    parser.add_argument("--source-date-epoch", type=int, required=True)
    args = parser.parse_args()

    engine = args.engine.resolve()
    source_root = args.source_root.resolve()
    output = args.output.resolve()
    if not engine.is_file():
        raise PackageFailure(f"engine not found: {engine}")
    if not source_root.is_dir():
        raise PackageFailure(f"source root not found: {source_root}")
    if output.exists():
        raise PackageFailure(f"output already exists: {output}")
    if output.suffix.lower() != ".zip":
        raise PackageFailure("output must have a .zip suffix")
    source_ids = (
        args.engine_source_commit,
        args.engine_source_tree,
        args.package_source_commit,
        args.package_source_tree,
    )
    if not all(re_full_hex(value, 40) for value in source_ids):
        raise PackageFailure("source commits and trees must be full lowercase Git object IDs")

    payload: dict[str, bytes] = {
        "KOTH-Stockfish.exe": engine.read_bytes(),
    }
    for target, relative_source in PAYLOAD_SOURCES.items():
        source = source_root / relative_source
        if not source.is_file():
            raise PackageFailure(f"required source file missing: {source}")
        payload[target] = source.read_bytes()

    if any(name.lower().endswith(".nnue") for name in payload):
        raise PackageFailure("network bytes must not be bundled")

    manifest = {
        "schema_version": 1,
        "format": "KOTH_RUNTIME_V1",
        "project": "KOTH-Stockfish",
        "target": "windows-x86-64-sse2",
        "source": {
            "upstream": "official-stockfish/Stockfish",
            "upstream_base_commit": "8bc5caa2e4b1d4c189b1428e93158b10d3edb0b6",
            "engine_commit": args.engine_source_commit,
            "engine_tree": args.engine_source_tree,
            "package_commit": args.package_source_commit,
            "package_tree": args.package_source_tree,
            "source_date_epoch": args.source_date_epoch,
        },
        "engine": {
            "path": "KOTH-Stockfish.exe",
            "bytes": len(payload["KOTH-Stockfish.exe"]),
            "sha256": sha256_bytes(payload["KOTH-Stockfish.exe"]),
            "runtime_dependencies": ["KERNEL32.dll", "msvcrt.dll", "SHELL32.dll"],
        },
        "rule_profile": {
            "name": "KOTH_LICHESS_V1",
            "uci_variant": "kingofthehill",
            "goal_squares": ["d4", "e4", "d5", "e5"],
        },
        "network": {
            "included": False,
            "external": True,
            "canonical_filename": "kingofthehill-978b86d0e6a4.nnue",
            "compatibility_alias": "KOTH_v1.nnue",
            "alias_is_default": False,
            "bytes": NETWORK_BYTES,
            "sha256": NETWORK_SHA256,
        },
        "search_contract": {
            "name": "KOTH_SEARCH_T0",
            "transposition_table_authority": "disabled",
            "syzygy": "disabled",
            "bench": "tools/koth_bench.py",
        },
        "referee": {
            "path": "tools/koth_referee.py",
            "requirements": "tools/koth-referee-requirements.txt",
        },
    }
    payload["manifest.json"] = (
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    ).encode("utf-8")
    checksums = "".join(
        f"{sha256_bytes(data)} *{name}\n" for name, data in sorted(payload.items())
    ).encode("ascii")
    payload["SHA256SUMS"] = checksums

    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(
        output, mode="x", compression=zipfile.ZIP_DEFLATED, compresslevel=9
    ) as archive:
        for name, data in sorted(payload.items()):
            archive.writestr(
                zip_info(
                    f"{PACKAGE_ROOT}/{name}",
                    args.source_date_epoch,
                    executable=name.endswith(".exe") or name.endswith(".py"),
                ),
                data,
            )

    result = {
        "package": output.name,
        "bytes": output.stat().st_size,
        "sha256": sha256_file(output),
        "entries": len(payload),
        "network_included": False,
        "engine_sha256": manifest["engine"]["sha256"],
        "engine_source_commit": args.engine_source_commit,
        "engine_source_tree": args.engine_source_tree,
        "package_source_commit": args.package_source_commit,
        "package_source_tree": args.package_source_tree,
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


def re_full_hex(value: str, length: int) -> bool:
    return len(value) == length and all(character in "0123456789abcdef" for character in value)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, PackageFailure, ValueError, zipfile.BadZipFile) as error:
        print(f"FAIL {error}", file=sys.stderr)
        raise SystemExit(1)
