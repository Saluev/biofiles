from dataclasses import dataclass
from enum import StrEnum


class TranscriptClassCode(StrEnum):
    """Describes query transcript layout relative to reference transcript."""

    MATCH = "="
    """ Complete, exact match of intron chain. """

    CONTAINED_WITHIN_EXONS = "c"
    """ Contained in reference (intron-compatible). """

    CONTAINS_WITHIN_EXONS = "k"
    """ Containment of reference (inverse containment). """

    RETAINED_INTRONS = "m"
    """ Retained intron(s), all introns are matched or retained. """

    RETAINED_INTRONS_INCOMPLETE = "n"
    """ Retained intron(s), not all introns matched/covered."""

    JUNCTION_MATCH = "j"
    """ Multi-exon with at least one junction match, """

    SINGLE_EXON_TRANSFRAG = "e"
    """ Single exon transfrag partially covering an intron, possible pre-mRNA fragment. """

    SAME_STRAND_OVERLAP = "o"
    """ Other same strand overlap with reference exons. """

    OPPOSITE_STRAND_INTRON_MATCH = "s"
    """ Intron match on the opposite strand (likely a mapping error). """

    OPPOSITE_STRAND_EXON_OVERLAP = "x"
    """ Exonic overlap on the opposite strand (like o or e but on the opposite strand). """

    CONTAINED_WITHIN_INTRON = "i"
    """ Fully contained within reference intron. """

    CONTAINS_WITHIN_INTRON = "i"
    """ Contains a reference within its intron(s). """

    POLYMERASE_RUN_ON = "p"
    """ Possible polymerase run-on (no actual overlap). """

    REPEAT = "r"
    """ Repeat (at least 50% bases soft-masked). """

    UNKNOWN = "u"
    """ None of the above (unknown, intergenic). """


@dataclass(frozen=True)
class TranscriptMapping:
    reference_gene_id: str | None
    reference_transcript_id: str | None
    class_code: TranscriptClassCode
    query_gene_id: str
    query_transcript_id: str
    query_number_of_exons: int
    query_fpkm: float
    query_tpm: float
    query_coverage: float
    query_length_bp: int
    query_major_isoform_gene_id: str
    match_length_bp: int | None
