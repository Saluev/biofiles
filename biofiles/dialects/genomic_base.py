from biofiles.types.feature import Feature


class Gene(Feature):
    biotype: str
    transcripts: list["Transcript"]


class Transcript(Feature):
    biotype: str
    gene: Gene
    exons: list["Exon"]


class Exon(Feature):
    transcript: Transcript
    cdss: list["CDS"]


class CDS(Feature):
    exon: Exon


class UTR(Feature):
    transcript: Transcript
