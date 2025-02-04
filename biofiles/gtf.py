__all__ = ["GTFReader"]

import sys
from typing import Iterator

from biofiles.gff import GFFReader
from biofiles.types.feature import Gene, Exon, Feature, UTR


class GTFReader(GFFReader):
    def __iter__(self) -> Iterator[Feature]:
        yield from self._read_gff3()

    def _parse_attributes(self, line: str, attributes_str: str) -> dict[str, str]:
        return {
            k: v.strip('"')
            for part in attributes_str.strip(";").split(";")
            for k, v in (part.strip().split(None, 1),)
        }


if __name__ == "__main__":
    for path in sys.argv[1:]:
        with GTFReader(path) as r:
            total_features = 0
            annotated_genes = 0
            annotated_exons = 0
            annotated_utrs = 0
            parsed_genes = 0
            parsed_exons = 0
            parsed_utrs = 0
            for feature in r:
                total_features += 1
                annotated_genes += feature.type_ == "gene"
                annotated_exons += feature.type_ == "exon"
                annotated_utrs += "utr" in feature.type_.lower()
                parsed_genes += isinstance(feature, Gene)
                parsed_exons += isinstance(feature, Exon)
                parsed_utrs += isinstance(feature, UTR)
        print(
            f"{path}: {total_features} features, "
            f"{parsed_genes} genes parsed out of {annotated_genes}, "
            f"{parsed_exons} exons parsed out of {annotated_exons}, "
            f"{parsed_utrs} UTRs parsed out of {annotated_utrs}"
        )
