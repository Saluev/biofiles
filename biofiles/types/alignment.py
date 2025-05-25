from dataclasses import dataclass


__all__ = [
    "Alignment",
    "BAMFlag",
    "BAMTag",
    "CIGAR",
    "CIGAROpKind",
    "CIGAROperation",
    "ReferenceSequence",
]

from enum import IntFlag

from typing import Any, Literal


@dataclass(frozen=True)
class ReferenceSequence:
    id: str
    length: int


@dataclass(frozen=True, slots=True)
class BAMTag:
    tag: str
    value: Any


CIGAROpKind = Literal["M", "I", "D", "N", "S", "H", "P", "=", "X"]


@dataclass(frozen=True, slots=True)
class CIGAROperation:
    kind: CIGAROpKind
    count: int


@dataclass(frozen=True)
class CIGAR:
    operations: tuple[CIGAROperation, ...]

    def __repr__(self) -> str:
        return f'CIGAR("{self}")'

    def __str__(self) -> str:
        return "".join(f"{op.count}{op.kind}" for op in self.operations)


class BAMFlag(IntFlag):
    MULTIPLE_SEGMENTS = 1 << 0
    EACH_SEGMENT_PROPERLY_ALIGNED = 1 << 1
    SEGMENT_UNMAPPED = 1 << 2
    NEXT_SEGMENT_UNMAPPED = 1 << 3
    READ_SEQUENCE_REVERSE_COMPLEMENTED = 1 << 4
    NEXT_SEGMENT_READ_SEQUENCE_REVERSE_COMPLEMENTED = 1 << 5
    FIRST_SEGMENT = 1 << 6
    LAST_SEGMENT = 1 << 7
    SECONDARY_SEGMENT = 1 << 8
    NOT_PASSING_QUALITY_CONTROL = 1 << 9
    DUPLICATE = 1 << 10
    SUPPLEMENTARY_ALIGNMENT = 1 << 11


@dataclass(frozen=True)
class Alignment:
    reference_sequence: ReferenceSequence | None

    start_c: int
    # 0-based leftmost coordinate.
    read_name: str
    mapping_quality: int
    bai_index_bin: int

    next_reference_sequence: ReferenceSequence | None
    next_start_c: int
    template_length: int
    cigar: CIGAR
    read_sequence: str
    quality: str

    bam_flags: int
    bam_tags: tuple[BAMTag, ...]
