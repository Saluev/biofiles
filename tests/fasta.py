import pathlib
from io import StringIO

from biofiles.fasta import FASTAReader, FASTASequence, FASTAWriter


def test_read_single_sequence_from_string() -> None:
    path = pathlib.Path(__file__).parent / "files" / "single_sequence.fasta"
    with FASTAReader(f"{path}") as r:
        sequences = [*r]

    assert sequences == [
        FASTASequence(
            id="SEQ", description="A very important sequence.", sequence="AT" * 400
        )
    ]


def test_read_single_sequence_from_path() -> None:
    path = pathlib.Path(__file__).parent / "files" / "single_sequence.fasta"
    with FASTAReader(path) as r:
        sequences = [*r]

    assert sequences == [
        FASTASequence(
            id="SEQ", description="A very important sequence.", sequence="AT" * 400
        )
    ]


def test_read_single_sequence_from_file() -> None:
    path = pathlib.Path(__file__).parent / "files" / "single_sequence.fasta"
    with open(path) as f, FASTAReader(f) as r:
        sequences = [*r]

    assert sequences == [
        FASTASequence(
            id="SEQ", description="A very important sequence.", sequence="AT" * 400
        )
    ]


def test_read_multiple_sequences_from_file() -> None:
    path = pathlib.Path(__file__).parent / "files" / "multiple_sequences.fasta"
    with open(path) as f, FASTAReader(f) as r:
        sequences = [*r]

    assert sequences == [
        FASTASequence(id="SEQ1", description="Goose", sequence="GAGAGA"),
        FASTASequence(id="SEQ2", description="Walker", sequence="ATAT"),
        FASTASequence(id="SEQ3", description="Moon landing", sequence="CG"),
    ]


def test_read_no_description() -> None:
    with FASTAReader(StringIO(">SEQ\nATGC")) as r:
        sequences = [*r]

    assert sequences == [FASTASequence(id="SEQ", description="", sequence="ATGC")]


def test_write_short_single_sequence() -> None:
    io = StringIO()
    w = FASTAWriter(io)
    w.write(FASTASequence(id="SEQ", description="A sequence.", sequence="ATGC"))
    assert io.getvalue() == ">SEQ A sequence.\nATGC\n"


def test_write_long_single_sequence() -> None:
    io = StringIO()
    w = FASTAWriter(io, width=2)
    w.write(FASTASequence(id="SEQ", description="A sequence.", sequence="ATGC"))
    assert io.getvalue() == ">SEQ A sequence.\nAT\nGC\n"


def test_write_multiple_sequences() -> None:
    io = StringIO()
    w = FASTAWriter(io)
    w.write(FASTASequence(id="SEQ1", description="Goose", sequence="GAGAGA"))
    w.write(FASTASequence(id="SEQ2", description="Walker", sequence="ATAT"))
    w.write(FASTASequence(id="SEQ3", description="Moon landing", sequence="CG"))
    assert (
        io.getvalue()
        == ">SEQ1 Goose\nGAGAGA\n>SEQ2 Walker\nATAT\n>SEQ3 Moon landing\nCG\n"
    )
