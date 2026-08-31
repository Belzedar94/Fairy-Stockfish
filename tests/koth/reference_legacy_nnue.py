#!/usr/bin/env python3
"""Independent scalar reference for KOTH_LEGACY_NNUE_V1.

The container and feature geometry are derived from official Stockfish commit
e8d64af1230fdac65bb0da246df3e7abe82e0838. This module intentionally shares no
production C++ evaluator code.
"""

from __future__ import annotations

import hashlib
import mmap
import struct
from dataclasses import dataclass
from pathlib import Path


FILE_BYTES = 47_721_371
SHA256 = "978B86D0E6A45E05F9F1375DCED129CEA0ACEA13041EA65960691632EDC47AF7"
VERSION = 0x7AF32F20
ARCHITECTURE = 0x3C103E72
FEATURE_HASH = 0x5F2348B8
LAYER_HASH = 0x633376CA
DESCRIPTION = b"Network trained with the https://github.com/glinscott/nnue-pytorch trainer."
INPUT_DIMENSIONS = 45_056
L1 = 512
BUCKETS = 8
FIRST_LAYER_OFFSET = 47_580_251
LAYER_BYTES = 17_640
OUTPUT_SCALE = 16
LAZY_THRESHOLD = 1400

PIECE_CODE = {
    "P": 1,
    "N": 2,
    "B": 3,
    "R": 4,
    "Q": 5,
    "K": 6,
    "p": 9,
    "n": 10,
    "b": 11,
    "r": 12,
    "q": 13,
    "k": 14,
}

PIECE_SQUARE_INDEX = (
    (0, 0, 128, 256, 384, 512, 640, 0, 0, 64, 192, 320, 448, 576, 640, 0),
    (0, 64, 192, 320, 448, 576, 640, 0, 0, 0, 128, 256, 384, 512, 640, 0),
)


@dataclass(frozen=True)
class Position:
    pieces: tuple[tuple[int, int], ...]
    side_to_move: int
    king_squares: tuple[int, int]


@dataclass(frozen=True)
class Evaluation:
    bucket: int
    psqt_raw: int
    positional_raw: int

    @property
    def lazy(self) -> bool:
        return abs(self.psqt_raw) > LAZY_THRESHOLD * OUTPUT_SCALE

    @property
    def psqt(self) -> int:
        return trunc_div(self.psqt_raw, OUTPUT_SCALE)

    @property
    def positional(self) -> int:
        return self.total - self.psqt

    @property
    def total(self) -> int:
        if self.lazy:
            return self.psqt
        return trunc_div(self.psqt_raw + self.positional_raw, OUTPUT_SCALE)


def trunc_div(value: int, divisor: int) -> int:
    """C++ integer division, including negative values."""
    return value // divisor if value >= 0 else -((-value) // divisor)


def wrap_signed(value: int, bits: int) -> int:
    """Two's-complement lane addition used by the historical SIMD accumulator."""
    modulus = 1 << bits
    sign = 1 << (bits - 1)
    return (value + sign) % modulus - sign


def parse_fen(fen: str) -> Position:
    fields = fen.split()
    if len(fields) != 6:
        raise ValueError("expected canonical six-field FEN")
    if fields[1] not in {"w", "b"}:
        raise ValueError("invalid side to move")

    pieces: list[tuple[int, int]] = []
    kings: list[int | None] = [None, None]
    ranks = fields[0].split("/")
    if len(ranks) != 8:
        raise ValueError("invalid board rank count")

    for fen_rank, encoded in enumerate(ranks):
        board_rank = 7 - fen_rank
        file_index = 0
        for token in encoded:
            if token.isdigit():
                file_index += int(token)
                continue
            if token not in PIECE_CODE or file_index >= 8:
                raise ValueError("invalid board token")
            square = board_rank * 8 + file_index
            piece = PIECE_CODE[token]
            pieces.append((square, piece))
            if token == "K":
                if kings[0] is not None:
                    raise ValueError("duplicate white king")
                kings[0] = square
            elif token == "k":
                if kings[1] is not None:
                    raise ValueError("duplicate black king")
                kings[1] = square
            file_index += 1
        if file_index != 8:
            raise ValueError("invalid board width")

    if any(square is None for square in kings):
        raise ValueError("missing king")
    if not 1 <= len(pieces) <= 32:
        raise ValueError("unsupported piece count")

    return Position(
        tuple(pieces),
        0 if fields[1] == "w" else 1,
        (int(kings[0]), int(kings[1])),
    )


class LegacyNetwork:
    def __init__(self, path: Path):
        self.path = path
        if path.stat().st_size != FILE_BYTES:
            raise ValueError("wrong network size")
        if hashlib.sha256(path.read_bytes()).hexdigest().upper() != SHA256:
            raise ValueError("wrong network digest")

        self._file = path.open("rb")
        self.data = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)
        self._validate_layout()

    def close(self) -> None:
        self.data.close()
        self._file.close()

    def __enter__(self) -> "LegacyNetwork":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def _u32(self, offset: int) -> int:
        return struct.unpack_from("<I", self.data, offset)[0]

    def _validate_layout(self) -> None:
        if self._u32(0) != VERSION:
            raise ValueError("wrong container version")
        if self._u32(4) != ARCHITECTURE:
            raise ValueError("wrong architecture hash")
        description_size = self._u32(8)
        if description_size != len(DESCRIPTION):
            raise ValueError("wrong description size")
        if self.data[12 : 12 + description_size] != DESCRIPTION:
            raise ValueError("wrong description")
        if self._u32(12 + description_size) != FEATURE_HASH:
            raise ValueError("wrong feature-transformer hash")
        for bucket in range(BUCKETS):
            offset = FIRST_LAYER_OFFSET + bucket * LAYER_BYTES
            if self._u32(offset) != LAYER_HASH:
                raise ValueError(f"wrong layer hash for bucket {bucket}")
        if FIRST_LAYER_OFFSET + BUCKETS * LAYER_BYTES != FILE_BYTES:
            raise AssertionError("layout does not reach exact EOF")

    @staticmethod
    def _orient(perspective: int, square: int) -> int:
        return square ^ (56 if perspective else 0)

    @staticmethod
    def _clip(value: int) -> int:
        return max(0, min(127, trunc_div(value, 64)))

    def evaluate(self, fen: str, forced_bucket: int | None = None) -> Evaluation:
        position = parse_fen(fen)
        piece_count = len(position.pieces)
        bucket = (piece_count - 1) // 4 if forced_bucket is None else forced_bucket
        if not 0 <= bucket < BUCKETS:
            raise ValueError("invalid bucket")

        feature_bias_offset = 91
        feature_weight_offset = feature_bias_offset + L1 * 2
        psqt_weight_offset = feature_weight_offset + INPUT_DIMENSIONS * L1 * 2
        biases = struct.unpack_from(f"<{L1}h", self.data, feature_bias_offset)

        accumulations: list[list[int]] = []
        psqt: list[list[int]] = []
        for perspective in (0, 1):
            accumulation = list(biases)
            psqt_values = [0] * BUCKETS
            king_square = self._orient(perspective, position.king_squares[perspective])
            for square, piece in position.pieces:
                index = (
                    self._orient(perspective, square)
                    + PIECE_SQUARE_INDEX[perspective][piece]
                    + 704 * king_square
                )
                if not 0 <= index < INPUT_DIMENSIONS:
                    raise AssertionError("feature index out of bounds")
                weights = struct.unpack_from(
                    f"<{L1}h", self.data, feature_weight_offset + index * L1 * 2
                )
                for output, weight in enumerate(weights):
                    accumulation[output] = wrap_signed(
                        accumulation[output] + weight, 16
                    )
                values = struct.unpack_from(
                    f"<{BUCKETS}i", self.data, psqt_weight_offset + index * BUCKETS * 4
                )
                for output, weight in enumerate(values):
                    psqt_values[output] = wrap_signed(
                        psqt_values[output] + weight, 32
                    )
            accumulations.append(accumulation)
            psqt.append(psqt_values)

        perspectives = (position.side_to_move, 1 - position.side_to_move)
        transformed = bytearray()
        for perspective in perspectives:
            transformed.extend(max(0, min(127, value)) for value in accumulations[perspective])

        psqt_raw = trunc_div(
            psqt[perspectives[0]][bucket] - psqt[perspectives[1]][bucket], 2
        )

        offset = FIRST_LAYER_OFFSET + bucket * LAYER_BYTES + 4
        fc0_biases = struct.unpack_from("<16i", self.data, offset)
        offset += 16 * 4
        fc0_weights = self.data[offset : offset + 16 * 1024]
        offset += 16 * 1024
        fc0 = []
        for output in range(16):
            row = struct.unpack_from("<1024b", fc0_weights, output * 1024)
            fc0.append(fc0_biases[output] + sum(w * x for w, x in zip(row, transformed)))
        ac0 = [self._clip(value) for value in fc0]

        fc1_biases = struct.unpack_from("<32i", self.data, offset)
        offset += 32 * 4
        fc1_weights = self.data[offset : offset + 32 * 32]
        offset += 32 * 32
        fc1 = []
        for output in range(32):
            row = struct.unpack_from("<32b", fc1_weights, output * 32)
            fc1.append(fc1_biases[output] + sum(w * x for w, x in zip(row[:16], ac0)))
        ac1 = [self._clip(value) for value in fc1]

        output_bias = struct.unpack_from("<i", self.data, offset)[0]
        offset += 4
        output_weights = struct.unpack_from("<32b", self.data, offset)
        offset += 32
        if offset != FIRST_LAYER_OFFSET + (bucket + 1) * LAYER_BYTES:
            raise AssertionError("layer parsing did not reach exact boundary")
        positional_raw = output_bias + sum(w * x for w, x in zip(output_weights, ac1))

        return Evaluation(bucket, psqt_raw, positional_raw)
