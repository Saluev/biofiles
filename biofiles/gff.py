import sys
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, cast, TextIO

from biofiles.common import Strand, Reader, Writer
from biofiles.types.feature import Feature, Gene, Exon

__all__ = ["GFFReader", "GFF3Writer"]


@dataclass
class _FeatureDraft:
    idx: int
    sequence_id: str
    source: str
    type_: str
    start_original: int
    end_original: int
    score: float | None
    strand: Strand | None
    phase: int | None
    attributes: dict[str, str]

    def pick_attribute(self, *keys: str) -> str | None:
        for key in keys:
            if (value := self.attributes.get(key, None)) is not None:
                return value
        return None


@dataclass
class _FeatureDrafts:
    drafts: deque[_FeatureDraft] = field(default_factory=deque)
    by_id: dict[str, _FeatureDraft] = field(default_factory=dict)
    # deps: dict[int, deque[int]] = field(default_factory=lambda: defaultdict(deque))

    def add(self, draft: _FeatureDraft) -> None:
        self.drafts.append(draft)
        if id_ := draft.attributes.get("ID", None):
            self.by_id[id_] = draft
        # if parent_id := draft.attributes.get("Parent", None):
        #     parent = self.by_id[parent_id]
        #     self.deps[parent.idx].append(draft.idx)

    # def remove_first_n(self, n: int) -> None:
    #     for _ in range(n):
    #         draft = self.drafts.popleft()
    #         if id_ := draft.attributes.get("ID", None):
    #             del self.by_id[id_]
    #         self.deps.pop(draft.idx, None)


@dataclass
class _Features:
    features: list[Feature] = field(default_factory=list)
    by_id: dict[str, Feature] = field(default_factory=dict)

    def add(self, feature: Feature):
        self.features.append(feature)
        if id_ := feature.attributes.get("ID", None):
            self.by_id[id_] = feature


class GFFReader(Reader):
    def __init__(
        self, input_: TextIO | Path | str, /, streaming_window: int | None = 1000
    ):
        super().__init__(input_)
        self._streaming_window = streaming_window

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
        drafts = _FeatureDrafts()
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

            draft = _FeatureDraft(
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

    def _finalize_drafts(
        self, drafts: _FeatureDrafts, w: int | None
    ) -> Iterator[Feature]:
        # TODO streaming version!
        #      code below is already tracking
        # if not drafts.drafts:
        #     return
        # if w is not None and len(drafts.drafts) <= w:
        #     return
        #
        # end_idx = drafts.drafts[-w].idx if w is not None else drafts.drafts[-1].idx
        #
        # i = 0
        # while i < len(drafts.drafts) and (
        #     not drafts.deps[drafts.drafts[i].idx]
        #     or drafts.deps[drafts.drafts[i].idx][-1] <= end_idx
        # ):
        #     i += 1
        #
        # print(f"FINALIZING {i} DRAFTS OUT OF {len(drafts.drafts)}")
        #
        # result = _Features()
        # for j in range(i):
        #     draft = drafts.drafts[j]
        #     feature = self._finalize_draft(draft, result)
        #     result.add(feature)
        # drafts.remove_first_n(i)
        # yield from result.features

        result = _Features()
        for draft in drafts.drafts:
            feature = self._finalize_draft(draft, result)
            result.add(feature)
        yield from result.features

    def _finalize_draft(self, draft: _FeatureDraft, result: _Features) -> Feature:
        match draft.type_:
            case "gene":
                feature = self._finalize_gene(draft, result)
            case "exon":
                feature = self._finalize_exon(draft, result)
            case _:
                feature = self._finalize_other(draft, result)
        if feature.parent:
            new_children = feature.parent.children + (feature,)
            object.__setattr__(feature.parent, "children", new_children)
        return feature

    def _finalize_gene(self, draft: _FeatureDraft, result: _Features) -> Feature:
        feature = self._finalize_other(draft, result)
        name = draft.pick_attribute("gene_name", "Name")
        biotype = draft.pick_attribute("gene_biotype", "biotype")
        if name is None or biotype is None:
            return feature
        return Gene(**feature.__dict__, name=name, biotype=biotype, exons=())

    def _finalize_exon(self, draft: _FeatureDraft, result: _Features) -> Feature:
        feature = self._finalize_other(draft, result)

        gene = feature.parent
        while gene and not isinstance(gene, Gene):
            gene = gene.parent

        if gene is None:
            return feature
        exon = Exon(**feature.__dict__, gene=gene)
        object.__setattr__(gene, "exons", gene.exons + (exon,))
        return exon

    def _finalize_other(self, draft: _FeatureDraft, result: _Features) -> Feature:
        parent_id = draft.attributes.get("Parent", None)
        parent = result.by_id[parent_id] if parent_id is not None else None

        return Feature(
            sequence_id=draft.sequence_id,
            source=draft.source,
            type_=draft.type_,
            start_original=draft.start_original,
            end_original=draft.end_original,
            start_c=draft.start_original - 1,
            end_c=draft.end_original,
            score=draft.score,
            strand=draft.strand,
            phase=draft.phase,
            attributes=draft.attributes,
            parent=parent,
            children=(),
        )

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
