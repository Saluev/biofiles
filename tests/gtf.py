import pathlib

from biofiles.dialects.genomic_base import Gene, Transcript, Exon
from biofiles.dialects.refseq import REFSEQ_FEATURE_TYPES
from biofiles.gtf import GTFReader


def test_parse_refseq_annotation() -> None:
    path = pathlib.Path(__file__).parent / "files" / "refseq_annotation.gtf"
    with GTFReader(path, REFSEQ_FEATURE_TYPES) as r:
        features = [*r]
    assert sum(1 for f in features if isinstance(f, Gene)) == 1
    assert sum(1 for f in features if isinstance(f, Transcript)) == 1
    assert sum(1 for f in features if isinstance(f, Exon)) == 3
