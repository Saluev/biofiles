__all__ = ["GTFReader", "GTFWriter"]

import sys
from typing import Iterator

from biofiles.common import Writer
from biofiles.gff import GFFReader
from biofiles.types.feature import Gene, Exon, Feature, UTR


class GTFReader(GFFReader):
    def __iter__(self) -> Iterator[Feature]:
        yield from self._read_gff3()

    def _parse_attributes(self, line: str, attributes_str: str) -> dict[str, str]:
        return {
            k: v.removeprefix('"').removesuffix('"').replace(r"\"", '"')
            for part in attributes_str.strip(";").split(";")
            for k, v in (part.strip().split(None, 1),)
        }


class GTFWriter(Writer):
    def write(self, feature: Feature) -> None:
        fields = (
            feature.sequence_id,
            feature.source,
            feature.type_,
            str(feature.start_c + 1),
            str(feature.end_c),
            str(feature.score) if feature.score is not None else ".",
            str(feature.strand) if feature.strand is not None else ".",
            str(feature.phase) if feature.phase is not None else ".",
            "; ".join(
                f'{k} "' + v.replace('"', r"\"") + '"'
                for k, v in feature.attributes.items()
            ),
        )
        self._output.write("\t".join(fields))
        self._output.write("\n")


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
                annotated_genes += "gene" in feature.type_.lower()
                annotated_exons += feature.type_ == "exon"
                annotated_utrs += "utr" in feature.type_.lower()
                parsed_genes += isinstance(feature, Gene)
                parsed_exons += isinstance(feature, Exon)
                parsed_utrs += isinstance(feature, UTR)
        print(
            f"{path}: {total_features} features, "
            f"{parsed_genes} genes parsed out of {annotated_genes}, "
            f"{parsed_exons} exons parsed out of {annotated_exons}, "
            f"{parsed_utrs} UTRs parsed out of {annotated_utrs}",
            file=sys.stderr,
        )
