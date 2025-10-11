"""Feature dialect for GENCODE .gtf/.gff3 files."""

from enum import StrEnum
from types import NoneType

from biofiles.types.feature_v2 import Feature, id_field, field, relation, no_id_field


class GeneType(StrEnum):
    ARTIFACT = "artifact"
    IG_C_GENE = "IG_C_gene"
    IG_C_PSEUDOGENE = "IG_C_pseudogene"
    IG_D_GENE = "IG_D_gene"
    IG_D_PSEUDOGENE = "IG_D_pseudogene"
    IG_J_GENE = "IG_J_gene"
    IG_J_PSEUDOGENE = "IG_J_pseudogene"
    IG_PSEUDOGENE = "IG_pseudogene"
    IG_V_GENE = "IG_V_gene"
    IG_V_PSEUDOGENE = "IG_V_pseudogene"
    LNCRNA = "lncRNA"
    MIRNA = "miRNA"
    MISC_RNA = "misc_RNA"
    MT_RRNA = "Mt_rRNA"
    MT_TRNA = "Mt_tRNA"
    PROCESSED_PSEUDOGENE = "processed_pseudogene"
    PROTEIN_CODING = "protein_coding"
    RIBOZYME = "ribozyme"
    RRNA = "rRNA"
    RRNA_PSEUDOGENE = "rRNA_pseudogene"
    SRNA = "sRNA"
    SCRNA = "scRNA"
    SCARNA = "scaRNA"
    SNRNA = "snRNA"
    SNORNA = "snoRNA"
    TEC = "TEC"
    TR_C_GENE = "TR_C_gene"
    TR_C_PSEUDOGENE = "TR_C_pseudogene"
    TR_D_GENE = "TR_D_gene"
    TR_D_PSEUDOGENE = "TR_D_pseudogene"
    TR_J_GENE = "TR_J_gene"
    TR_J_PSEUDOGENE = "TR_J_pseudogene"
    TR_V_GENE = "TR_V_gene"
    TR_V_PSEUDOGENE = "TR_V_pseudogene"
    TRANSCRIBED_PROCESSED_PSEUDOGENE = "transcribed_processed_pseudogene"
    TRANSCRIBED_UNITARY_PSEUDOGENE = "transcribed_unitary_pseudogene"
    TRANSCRIBED_UNPROCESSED_PSEUDOGENE = "transcribed_unprocessed_pseudogene"
    TRANSLATED_PROCESSED_PSEUDOGENE = "translated_processed_pseudogene"
    UNITARY_PSEUDOGENE = "unitary_pseudogene"
    UNPROCESSED_PSEUDOGENE = "unprocessed_pseudogene"
    VAULT_RNA = "vault_RNA"


class TranscriptType(StrEnum):
    ARTIFACT = "artifact"
    IG_C_GENE = "IG_C_gene"
    IG_C_PSEUDOGENE = "IG_C_pseudogene"
    IG_D_GENE = "IG_D_gene"
    IG_D_PSEUDOGENE = "IG_D_pseudogene"
    IG_J_GENE = "IG_J_gene"
    IG_J_PSEUDOGENE = "IG_J_pseudogene"
    IG_PSEUDOGENE = "IG_pseudogene"
    IG_V_GENE = "IG_V_gene"
    IG_V_PSEUDOGENE = "IG_V_pseudogene"
    LNCRNA = "lncRNA"
    MIRNA = "miRNA"
    MISC_RNA = "misc_RNA"
    MT_RRNA = "Mt_rRNA"
    MT_TRNA = "Mt_tRNA"
    PROCESSED_PSEUDOGENE = "processed_pseudogene"
    PROTEIN_CODING = "protein_coding"
    RIBOZYME = "ribozyme"
    RRNA = "rRNA"
    RRNA_PSEUDOGENE = "rRNA_pseudogene"
    SRNA = "sRNA"
    SCRNA = "scRNA"
    SCARNA = "scaRNA"
    SNRNA = "snRNA"
    SNORNA = "snoRNA"
    TEC = "TEC"
    TR_C_GENE = "TR_C_gene"
    TR_C_PSEUDOGENE = "TR_C_pseudogene"
    TR_D_GENE = "TR_D_gene"
    TR_D_PSEUDOGENE = "TR_D_pseudogene"
    TR_J_GENE = "TR_J_gene"
    TR_J_PSEUDOGENE = "TR_J_pseudogene"
    TR_V_GENE = "TR_V_gene"
    TR_V_PSEUDOGENE = "TR_V_pseudogene"
    TRANSCRIBED_PROCESSED_PSEUDOGENE = "transcribed_processed_pseudogene"
    TRANSCRIBED_UNITARY_PSEUDOGENE = "transcribed_unitary_pseudogene"
    TRANSCRIBED_UNPROCESSED_PSEUDOGENE = "transcribed_unprocessed_pseudogene"
    TRANSLATED_PROCESSED_PSEUDOGENE = "translated_processed_pseudogene"
    UNITARY_PSEUDOGENE = "unitary_pseudogene"
    UNPROCESSED_PSEUDOGENE = "unprocessed_pseudogene"
    VAULT_RNA = "vault_RNA"

    # Transcript-specific:
    NON_STOP_DECAY = "non_stop_decay"
    NONSENSE_MEDIATED_DECAY = "nonsense_mediated_decay"
    PROCESSED_TRANSCRIPT = "processed_transcript"
    PROTEIN_CODING_CDS_NOT_DEFINED = "protein_coding_CDS_not_defined"
    PROTEIN_CODING_LOF = "protein_coding_LoF"
    RETAINED_INTRON = "retained_intron"


transcript_gene, gene_transcripts = relation(source="gene_id")
selenocysteine_gene, _ = relation(source="gene_id")
selenocysteine_transcript, _ = relation(source="transcript_id")
exon_transcript, transcript_exons = relation(source="transcript_id")
exon_gene, _ = relation(source="gene_id")
cds_exon, exon_cds = relation(source=("transcript_id", "exon_number"), one_to_one=True)
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
start_codon_exon, _ = relation(source=("transcript_id", "exon_number"), one_to_one=True)
stop_codon_transcript, transcript_stop_codon = relation(
    source="transcript_id", one_to_one=True
)
stop_codon_exon, _ = relation(source=("transcript_id", "exon_number"), one_to_one=True)


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


class Selenocysteine(
    Feature, type=("selenocysteine", "stop_codon_redefined_as_selenocysteine")
):
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


class FivePrimeUTR(UTR, type="five_prime_utr"):
    id: NoneType = no_id_field()
    transcript: Transcript = five_prime_utr_transcript
    gene: Gene = five_prime_utr_gene


class ThreePrimeUTR(UTR, type="three_prime_utr"):
    id: NoneType = no_id_field()
    transcript: Transcript = three_prime_utr_transcript
    gene: Gene = three_prime_utr_gene


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
    FivePrimeUTR,
    ThreePrimeUTR,
    StartCodon,
    StopCodon,
]
