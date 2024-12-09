from dataclasses import dataclass

from biofiles.common import Strand


__all__ = ["Repeat"]


@dataclass(frozen=True)
class Repeat:
    sw_score: int
    divergence_percent: float
    deletion_percent: float
    insertion_percent: float
    sequence_id: str
    sequence_start_original: int
    sequence_end_original: int
    sequence_start_c: int
    sequence_end_c: int
    sequence_left: int
    strand: Strand
    repeat_name: str
    repeat_class: str
    repeat_family: str | None
    repeat_start_original: int
    repeat_end_original: int
    repeat_start_c: int
    repeat_end_c: int
    repeat_left: int
    repeat_id: str | None
