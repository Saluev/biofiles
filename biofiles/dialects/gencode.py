"""Feature dialect for GENCODE .gtf/.gff3 files."""

from enum import StrEnum
from types import NoneType

from biofiles.types.feature_v2 import Feature, id_field, field, relation, no_id_field


class GeneType(StrEnum):
    LNC_RNA = "lncRNA"
    PROTEIN_CODING = "protein_coding"


class TranscriptType(StrEnum):
    LNC_RNA = "lncRNA"
    PROTEIN_CODING = "protein_coding"


transcript_gene, gene_transcripts = relation(source="gene_id")
selenocysteine_gene, _ = relation(source="gene_id")
selenocysteine_transcript, _ = relation(source="transcript_id")
exon_transcript, transcript_exons = relation(source="transcript_id")
exon_gene, _ = relation(source="gene_id")
cds_exon, exon_cds = relation(source="exon_id", one_to_one=True)
utr_transcript, transcript_utrs = relation(source="transcript_id")
utr_gene, _ = relation(source="gene_id")
five_prime_utr_transcript, transcript_five_prime_utr = relation(
    source="transcript_id", one_to_one=True
)
five_prime_utr_gene, _ = relation(source="gene_id")
three_prime_utr_transcript, transcript_three_prime_utr = relation(
    source="transcript_id", one_to_one=True
)
three_prime_utr_gene, _ = relation(source="gene_id")
start_codon_transcript, transcript_start_codon = relation(
    source="transcript_id", one_to_one=True
)
start_codon_exon, _ = relation(source="exon_id", one_to_one=True)
stop_codon_transcript, transcript_stop_codon = relation(
    source="transcript_id", one_to_one=True
)
stop_codon_exon, _ = relation(source="exon_id", one_to_one=True)


class Gene(Feature, type="gene"):
    id: str = id_field(source="gene_id")
    type: GeneType = field(source="gene_type")
    name: str = field(source="gene_name")
    transcripts: list["Transcript"] = gene_transcripts


class Transcript(Feature, type="transcript"):
    id: str = id_field(source="transcript_id")
    type: TranscriptType = field(source="transcript_type")
    name: str = field(source="transcript_name")
    gene: Gene = transcript_gene
    exons: list["Exon"] = transcript_exons
    utrs: list["UTR"] = transcript_utrs
    start_codon: "StartCodon | None" = transcript_start_codon
    stop_codon: "StopCodon | None" = transcript_stop_codon


class Selenocysteine(Feature, type="selenocysteine"):
    id: str = no_id_field()
    gene: Gene = selenocysteine_gene
    transcript: Transcript = selenocysteine_transcript


class Exon(Feature, type="exon"):
    id: tuple[str, int] = id_field(source=("transcript_id", "exon_number"))
    number: int = field(source="exon_number")
    transcript: Transcript = exon_transcript
    gene: Gene = exon_gene
    cds: "CDS | None" = exon_cds


class CDS(Feature, type="cds"):
    id: tuple[str, int] = id_field(source=("transcript_id", "exon_number"))
    exon: Exon = cds_exon


class UTR(Feature, type="utr"):
    id: NoneType = no_id_field()
    transcript: Transcript = utr_transcript
    gene: Gene = utr_gene


# GENCODE doesn't distinguish between 5' and 3' UTRs at the annotation level.


class StartCodon(Feature, type="start_codon"):
    id: tuple[str, int] = id_field(source=("transcript_id", "exon_number"))
    transcript: Transcript = start_codon_transcript
    exon: Exon = start_codon_exon


class StopCodon(Feature, type="stop_codon"):
    id: tuple[str, int] = id_field(source=("transcript_id", "exon_number"))
    transcript: Transcript = stop_codon_transcript
    exon: Exon = stop_codon_exon


GENCODE_FEATURE_TYPES = [
    Gene,
    Transcript,
    Selenocysteine,
    Exon,
    CDS,
    UTR,
    StartCodon,
    StopCodon,
]
