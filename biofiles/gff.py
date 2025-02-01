import sys
from pathlib import Path
from typing import Iterator, cast, TextIO

from biofiles.common import Strand, Writer
from biofiles.feature import FeatureReader, FeatureDraft, FeatureDrafts
from biofiles.types.feature import Feature, Gene, Exon

__all__ = ["GFFReader", "GFF3Writer"]


class GFFReader(FeatureReader):
    def __iter__(self) -> Iterator[Feature]:
        for line in self._input:
            line = line.rstrip("\n")
            if line.startswith(_VERSION_PREFIX):
                version = line.removeprefix(_VERSION_PREFIX)
                if version == "3":
                    yield from self._read_gff3()
                    return
                raise ValueError(f"unsupported version {version!r}")
            if line.startswith("#"):
                continue
            raise ValueError(f"unexpected line {line!r}, expected version")

    def _read_gff3(self) -> Iterator[Feature]:
        drafts = FeatureDrafts()
        idx = 0
        for line in self._input:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            parts = line.split("\t", maxsplit=8)
            if len(parts) != 9:
                raise ValueError(f"unexpected line {line!r}, expected 9 columns")
            (
                sequence_id,
                source,
                type_,
                start_str,
                end_str,
                score_str,
                strand_str,
                phase_str,
                attributes_str,
            ) = parts
            score = self._parse_score(line, score_str)
            strand = self._parse_strand(line, strand_str)
            phase = self._parse_phase(line, phase_str)
            attributes = self._parse_attributes(line, attributes_str)

            parent_id = attributes.get("Parent", None)
            # if parent_id is None:
            #     yield from self._finalize_drafts(drafts)
            #     drafts = _FeatureDrafts()
            if parent_id is not None and parent_id not in drafts.by_id:
                raise ValueError(
                    f"unexpected line {line!r}, parent ID not among recent feature IDs"
                )

            draft = FeatureDraft(
                idx=idx,
                sequence_id=sequence_id,
                source=source,
                type_=type_,
                start_original=int(start_str),
                end_original=int(end_str),
                score=score,
                strand=strand,
                phase=phase,
                attributes=attributes,
            )
            drafts.add(draft)
            idx += 1

            # yield from self._finalize_drafts(drafts, self._streaming_window)

        yield from self._finalize_drafts(drafts, None)

    def _parse_score(self, line: str, score_str: str) -> float | None:
        if score_str == ".":
            return None
        try:
            return float(score_str)
        except ValueError as exc:
            raise ValueError(
                f"unexpected line {line!r}, score should be a number or '.'"
            ) from exc

    def _parse_strand(self, line: str, strand_str: str) -> Strand | None:
        if strand_str in ("-", "+"):
            return cast(Strand, strand_str)
        if strand_str == ".":
            return None
        raise ValueError(f"unexpected line {line!r}, strand should be '-', '+' or '.'")

    def _parse_phase(self, line: str, phase_str: str) -> int | None:
        if phase_str == ".":
            return None
        try:
            return int(phase_str)
        except ValueError as exc:
            raise ValueError(
                f"unexpected line {line!r}, phase should be an integer or '.'"
            ) from exc

    def _parse_attributes(self, line: str, attributes_str: str) -> dict[str, str]:
        return {
            k: v
            for part in attributes_str.strip(";").split(";")
            for k, v in (part.split("=", 1),)
        }


class GFF3Writer(Writer):
    def __init__(self, output: TextIO | Path | str) -> None:
        super().__init__(output)
        self._output.write(f"{_VERSION_PREFIX}3\n")

    def write(self, feature: Feature) -> None:
        fields = (
            feature.sequence_id,
            feature.source,
            feature.type_,
            str(feature.start_original),
            str(feature.end_original),
            str(feature.score) if feature.score is not None else ".",
            str(feature.strand) if feature.strand is not None else ".",
            str(feature.phase) if feature.phase is not None else ".",
            ";".join(f"{k}={v}" for k, v in feature.attributes.items()),
        )
        self._output.write("\t".join(fields))
        self._output.write("\n")


_VERSION_PREFIX = "##gff-version "


if __name__ == "__main__":
    for path in sys.argv[1:]:
        with GFFReader(path) as r:
            total_features = 0
            annotated_genes = 0
            annotated_exons = 0
            parsed_genes = 0
            parsed_exons = 0
            for feature in r:
                total_features += 1
                annotated_genes += feature.type_ == "gene"
                annotated_exons += feature.type_ == "exon"
                parsed_genes += isinstance(feature, Gene)
                parsed_exons += isinstance(feature, Exon)
        print(
            f"{path}: {total_features} features, {parsed_genes} genes parsed out of {annotated_genes}, {parsed_exons} exons parsed out of {annotated_exons}"
        )
