import pathlib
from io import StringIO

import pytest

from biofiles.dialects.gencode import GENCODE_FEATURE_TYPES
from biofiles.gff import GFFReader, GFF3Writer
from biofiles.dialects.genomic_base import Gene, Transcript, Exon, Feature


def test_parse_gencode_annotation() -> None:
    path = pathlib.Path(__file__).parent / "files" / "gencode_49_annotation.gff"
    with GFFReader(path, GENCODE_FEATURE_TYPES) as r:
        features = [*r]
    assert sum(1 for f in features if isinstance(f, Gene)) == 1
    assert sum(1 for f in features if isinstance(f, Transcript)) == 1
    assert sum(1 for f in features if isinstance(f, Exon)) == 3


def test_parse_liftoff_annotation() -> None:
    pytest.skip("skipped until specific dialect for Liftoff is implemented")
    path = pathlib.Path(__file__).parent / "files" / "liftoff_annotation.gff"
    with GFFReader(path) as r:
        features = [*r]
    assert sum(1 for f in features if isinstance(f, Gene)) == 3
    assert sum(1 for f in features if isinstance(f, Exon)) == 11


def test_parse_gnomon_annotation() -> None:
    pytest.skip("skipped until specific dialect for Gnomon is implemented")
    path = pathlib.Path(__file__).parent / "files" / "gnomon_annotation.gff"
    with GFFReader(path) as r:
        features = [*r]
    assert sum(1 for f in features if isinstance(f, Gene)) == 2
    assert sum(1 for f in features if isinstance(f, Exon)) == 9


def test_write_gff3() -> None:
    io = StringIO()
    w = GFF3Writer(io)
    w.write(
        Feature(
            sequence_id="chr1",
            source="Imagination",
            type_="gene",
            start_original=70,
            end_original=420,
            start_c=69,
            end_c=420,
            score=None,
            strand="+",
            phase=None,
            attributes={"biotype": "protein_coding", "foo": "bar"},
            id=None,
            parent=None,
            children=(),
        )
    )
    assert (
        io.getvalue()
        == "##gff-version 3\nchr1\tImagination\tgene\t70\t420\t.\t+\t.\tbiotype=protein_coding;foo=bar\n"
    )
