"""Feature dialect for RefSeq .gtf/.gff3 files."""

from enum import StrEnum

from biofiles.dialects.genomic_base import (
    Gene as BaseGene,
    Transcript as BaseTranscript,
    Exon as BaseExon,
    CDS as BaseCDS,
)
from biofiles.types.feature_v2 import Feature, id_field, field, relation, no_id_field


class GeneType(StrEnum):
    ANTISENSE_RNA = "antisense_RNA"
    C_REGION = "C_region"
    C_REGION_PSEUDOGENE = "C_region_pseudogene"
    D_SEGMENT = "D_segment"
    D_SEGMENT_PSEUDOGENE = "D_segment_pseudogene"
    J_SEGMENT = "J_segment"
    J_SEGMENT_PSEUDOGENE = "J_segment_pseudogene"
    LNCRNA = "lncRNA"
    MIRNA = "miRNA"
    MISC_RNA = "misc_RNA"
    NCRNA = "ncRNA"
    NCRNA_PSEUDOGENE = "ncRNA_pseudogene"
    OTHER = "other"
    PROTEIN_CODING = "protein_coding"
    PSEUDOGENE = "pseudogene"
    RNASE_MRP_RNA = "RNase_MRP_RNA"
    RNASE_P_RNA = "RNase_P_RNA"
    RRNA = "rRNA"
    SCARNA = "scaRNA"
    SCRNA = "scRNA"
    SNORNA = "snoRNA"
    SNRNA = "snRNA"
    TELOMERASE_RNA = "telomerase_RNA"
    TRANSCRIBED_PSEUDOGENE = "transcribed_pseudogene"
    TRNA = "tRNA"
    VAULT_RNA = "vault_RNA"
    V_SEGMENT = "V_segment"
    V_SEGMENT_PSEUDOGENE = "V_segment_pseudogene"
    Y_RNA = "Y_RNA"


class TranscriptType(StrEnum):
    ANTISENSE_RNA = "antisense_RNA"
    C_GENE_SEGMENT = "C_gene_segment"
    D_GENE_SEGMENT = "D_gene_segment"
    J_GENE_SEGMENT = "J_gene_segment"
    LNC_RNA = "lnc_RNA"
    MIRNA = "miRNA"
    MRNA = "mRNA"
    PRIMARY_TRANSCRIPT = "primary_transcript"
    RNASE_MRP_RNA = "RNase_MRP_RNA"
    RNASE_P_RNA = "RNase_P_RNA"
    RRNA = "rRNA"
    SCARNA = "scaRNA"
    SCRNA = "scRNA"
    SNORNA = "snoRNA"
    SNRNA = "snRNA"
    TELOMERASE_RNA = "telomerase_RNA"
    TRANSCRIPT = "transcript"
    TRNA = "tRNA"
    VAULT_RNA = "vault_RNA"
    V_GENE_SEGMENT = "V_gene_segment"
    Y_RNA = "Y_RNA"


transcript_gene, gene_transcripts = relation(source="gene_id")
exon_transcript, transcript_exons = relation(source="transcript_id")
exon_gene, _ = relation(source="gene_id")
cds_exon, exon_cds = relation(source=("transcript_id", "exon_number"))
start_codon_transcript, transcript_start_codon = relation(
    source="transcript_id", one_to_one=True
)
start_codon_exon, _ = relation(source=("transcript_id", "exon_number"), one_to_one=True)
stop_codon_transcript, transcript_stop_codon = relation(
    source="transcript_id", one_to_one=True
)
stop_codon_exon, _ = relation(source=("transcript_id", "exon_number"), one_to_one=True)


class Gene(BaseGene, type="gene"):
    id: str = id_field(source="gene_id")
    type: GeneType = field(source="gene_biotype")
    name: str = field(source="gene")
    synonyms: list[str] = field(source="gene_synonym", default_factory=list)
    transcripts: list["Transcript"] = gene_transcripts


class Transcript(BaseTranscript, type="transcript"):
    id: str = id_field(source="transcript_id")
    type: TranscriptType = field(source="transcript_biotype")
    product: str | None = field(source="product", default=None)
    gene: Gene = transcript_gene
    exons: list["Exon"] = transcript_exons
    start_codon: "StartCodon | None" = transcript_start_codon
    stop_codon: "StopCodon | None" = transcript_stop_codon


class Exon(BaseExon, type="exon"):
    id: tuple[str, int] = id_field(source=("transcript_id", "exon_number"))
    number: int = field(source="exon_number")
    transcript: Transcript = exon_transcript
    gene: Gene = exon_gene
    cdss: list["CDS"] = exon_cds


class CDS(BaseCDS, type="cds"):
    id: tuple[str, int] = no_id_field()
    exon: Exon = cds_exon


class StartCodon(Feature, type="start_codon"):
    id: tuple[str, int] = id_field(source=("transcript_id", "exon_number"))
    transcript: Transcript = start_codon_transcript
    exon: Exon = start_codon_exon


class StopCodon(Feature, type="stop_codon"):
    id: tuple[str, int] = id_field(source=("transcript_id", "exon_number"))
    transcript: Transcript = stop_codon_transcript
    exon: Exon = stop_codon_exon


REFSEQ_FEATURE_TYPES = [
    Gene,
    Transcript,
    Exon,
    CDS,
    StartCodon,
    StopCodon,
]
