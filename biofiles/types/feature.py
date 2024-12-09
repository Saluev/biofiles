from dataclasses import dataclass

from biofiles.common import Strand


__all__ = ["Feature", "Gene", "Exon"]


@dataclass(frozen=True)
class Feature:
    sequence_id: str
    source: str
    type_: str

    start_original: int
    end_original: int
    # Original, 1-based inclusive values.

    start_c: int
    end_c: int
    # Standardized ("C-style") 0-based values, start inclusive, end exclusive.

    score: float | None
    strand: Strand | None
    phase: int | None
    attributes: dict[str, str]

    parent: "GFFFeature | None"
    children: tuple["Feature", ...]


# Custom types for particular kinds of features:


@dataclass(frozen=True)
class Gene(Feature):
    name: str
    biotype: str
    exons: tuple["Exon", ...]


@dataclass(frozen=True)
class Exon(Feature):
    gene: Gene
    # TODO transcript, mRNA
