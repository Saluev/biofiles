import sys
from itertools import islice
from pathlib import Path

from biofiles.dialects.gencode import GENCODE_DIALECT
from biofiles.dialects.refseq import REFSEQ_DIALECT
from biofiles.types.feature import Dialect
from biofiles.utility.feature import RawFeatureReader


class CantDetectDialect(Exception):
    pass


class DialectDetector:
    def __init__(self, raw_reader: RawFeatureReader, num_samples: int = 1000) -> None:
        self._raw_reader = raw_reader
        self._num_samples = num_samples

    def detect(self) -> Dialect:
        gencode_rows = 0
        refseq_rows = 0
        total_rows = 0
        for fd in islice(self._raw_reader, self._num_samples):
            total_rows += 1
            source = fd.source.lower()
            if source in ("havana", "ensembl"):
                gencode_rows += 1
            elif source in ("bestrefseq", "bestrefseq%2cgnomon", "gnomon", "refseq"):
                refseq_rows += 1

        if gencode_rows > 0 and gencode_rows >= 0.9 * total_rows:
            return GENCODE_DIALECT
        if refseq_rows > 0 and refseq_rows >= 0.9 * total_rows:
            return REFSEQ_DIALECT

        raise CantDetectDialect(
            f"of {total_rows} read rows {gencode_rows} look like GENCODE, {refseq_rows} look like RefSeq"
        )


def detect_dialect(path: Path) -> Dialect:
    if path.suffix == ".gtf":
        from biofiles.gtf import RawGTFReader

        raw_reader = RawGTFReader(path)
    elif path.suffix in (".gff", ".gff3"):
        from biofiles.gff import RawGFFReader

        raw_reader = RawGFFReader(path)
    else:
        raise CantDetectDialect(f"unknown file extension {path.suffix}")
    detector = DialectDetector(raw_reader=raw_reader)
    return detector.detect()


if __name__ == "__main__":
    exit_code = 0
    for path_str in sys.argv[1:]:
        path = Path(path_str)
        try:
            dialect = detect_dialect(path)
            print(f"{path}\t{dialect.name}")
        except CantDetectDialect as exc:
            print(f"Failed to detect dialect for {path}: {exc}", file=sys.stderr)
            exit_code = 1
    sys.exit(exit_code)
