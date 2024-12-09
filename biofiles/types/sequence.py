from dataclasses import dataclass


__all__ = ["Sequence"]


@dataclass(frozen=True)
class Sequence:
    id: str
    description: str
    sequence: str
