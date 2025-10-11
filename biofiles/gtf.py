__all__ = ["GTFReader", "GTFWriter"]

import sys
from pathlib import Path
from typing import Iterator

from biofiles.common import Writer
from biofiles.dialects.detector import detect_dialect
from biofiles.dialects.genomic_base import Gene, Exon, Feature, CDS, UTR
from biofiles.gff import RawGFFReader
from biofiles.utility.feature import FeatureReader, RawFeatureReader, FeatureDraft


class RawGTFReader(RawGFFReader):
    def __iter__(self) -> Iterator[FeatureDraft]:
        yield from self._read_gff3()

    def _parse_attributes(
        self, line: str, attributes_str: str
    ) -> dict[str, str | list[str]]:
        try:
            result: dict[str, str | list[str]] = {}
            for part in attributes_str.strip().strip(";").split(";"):
                k, v = part.strip().split(None, 1)
                v = v.removeprefix('"').removesuffix('"').replace(r"\"", '"')
                if k in result:
                    if not isinstance(result[k], list):
                        result[k] = [result[k]]
                    result[k].append(v)
                else:
                    result[k] = v
            return result
        except ValueError as exc:
            raise ValueError(
                f"failed to parse attribute string {attributes_str!r}: {exc}"
            ) from exc


class GTFReader(FeatureReader):

    def _make_raw_feature_reader(self) -> RawFeatureReader:
        return RawGTFReader(self._input)


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
        dialect = detect_dialect(Path(path))
        with GTFReader(path, dialect=dialect) as r:
            total_features = 0
            annotated_genes = 0
            annotated_exons = 0
            annotated_cds = 0
            annotated_utrs = 0
            parsed_genes = 0
            parsed_exons = 0
            parsed_cds = 0
            parsed_utrs = 0
            for feature in r:
                total_features += 1
                annotated_genes += "gene" in feature.type_.lower()
                annotated_exons += feature.type_ == "exon"
                annotated_cds += feature.type_.lower() == "cds"
                annotated_utrs += "utr" in feature.type_.lower()
                parsed_genes += isinstance(feature, Gene)
                parsed_exons += isinstance(feature, Exon)
                parsed_cds += isinstance(feature, CDS)
                parsed_utrs += isinstance(feature, UTR)
        print(
            f"{path}: {total_features} features, "
            f"{parsed_genes} genes parsed out of {annotated_genes}, "
            f"{parsed_exons} exons parsed out of {annotated_exons}, "
            f"{parsed_cds} CDS parsed out of {annotated_cds}, "
            f"{parsed_utrs} UTRs parsed out of {annotated_utrs}",
            file=sys.stderr,
        )
