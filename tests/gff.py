import pathlib

from biofiles.gff import GFFReader
from biofiles.types.feature import Gene, Exon


def test_parse_liftoff_annotation() -> None:
    path = pathlib.Path(__file__).parent / "files" / "liftoff_annotation.gff"
    with GFFReader(path) as r:
        features = [*r]
    assert sum(1 for f in features if isinstance(f, Gene)) == 3
    assert sum(1 for f in features if isinstance(f, Exon)) == 11


def test_parse_gnomon_annotation() -> None:
    path = pathlib.Path(__file__).parent / "files" / "gnomon_annotation.gff"
    with GFFReader(path) as r:
        features = [*r]
    assert sum(1 for f in features if isinstance(f, Gene)) == 2
    assert sum(1 for f in features if isinstance(f, Exon)) == 9
