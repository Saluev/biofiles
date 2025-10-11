from biofiles.types.feature import Feature


class Gene(Feature):
    id: str
    biotype: str
    name: str
    transcripts: list["Transcript"]


class Transcript(Feature):
    id: str
    biotype: str
    gene: Gene
    exons: list["Exon"]


class Exon(Feature):
    number: int
    transcript: Transcript
    cdss: list["CDS"]


class CDS(Feature):
    exon: Exon


class UTR(Feature):
    transcript: Transcript


class FivePrimeUTR(UTR):
    pass


class ThreePrimeUTR(UTR):
    pass
