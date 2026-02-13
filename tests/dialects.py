from pathlib import Path

from biofiles.dialects.detector import detect_dialect
from biofiles.dialects.gencode import GENCODE_DIALECT, Transcript
from biofiles.dialects.refseq import REFSEQ_DIALECT
from biofiles.gff import GFFReader


def test_detect_gencode_dialect() -> None:
    path = Path(__file__).parent / "files" / "gencode_49_annotation.gff"
    assert detect_dialect(path) is GENCODE_DIALECT


def test_detect_refseq_dialect() -> None:
    path = Path(__file__).parent / "files" / "refseq_annotation.gtf"
    assert detect_dialect(path) is REFSEQ_DIALECT


def test_parse_gencode_spliced_utr() -> None:
    path = Path(__file__).parent / "files" / "gencode_spliced_utr.gff"
    with GFFReader(path, GENCODE_DIALECT) as r:
        features = [*r]
    (transcript,) = (f for f in features if isinstance(f, Transcript))
    assert len(transcript.five_prime_utrs) == 2
