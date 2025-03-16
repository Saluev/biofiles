from dataclasses import dataclass


__all__ = ["Sequence"]


@dataclass(frozen=True)
class Sequence:
    id: str
    description: str
    sequence: str


@dataclass(frozen=True)
class SequenceDescription:
    id: str
    length: int
    byte_offset: int
