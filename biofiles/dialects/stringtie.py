"""Feature dialect for StringTie .gtf/.gff output."""

from biofiles.types.feature import Feature, relation, id_field, field, Dialect

exon_transcript, transcript_exons = relation(source="transcript_id")


class Transcript(Feature, type="transcript"):
    id: str = id_field(source="transcript_id")
    gene_id: str = field(source="gene_id")
    exons: list["Exon"] = transcript_exons
    coverage: float = field(source="cov")
    fpkm: float = field(source="FPKM")
    tpm: float = field(source="TPM")


class Exon(Feature, type="exon"):
    id: tuple[str, int] = id_field(source=("transcript_id", "exon_number"))
    transcript: Transcript = exon_transcript
    coverage: float = field(source="cov")


STRINGTIE_FEATURE_TYPES = [Transcript, Exon]
STRINGTIE_DIALECT = Dialect(name="StringTie", feature_types=STRINGTIE_FEATURE_TYPES)
