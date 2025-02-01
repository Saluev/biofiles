from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, TextIO

from biofiles.common import Reader, Strand
from biofiles.types.feature import Feature, Gene, ThreePrimeUTR, Exon


@dataclass
class FeatureDraft:
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
class FeatureDrafts:
    drafts: deque[FeatureDraft] = field(default_factory=deque)
    by_id: dict[str, FeatureDraft] = field(default_factory=dict)
    # deps: dict[int, deque[int]] = field(default_factory=lambda: defaultdict(deque))

    def add(self, draft: FeatureDraft) -> None:
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
class Features:
    features: list[Feature] = field(default_factory=list)
    by_id: dict[str, Feature] = field(default_factory=dict)

    def add(self, feature: Feature):
        self.features.append(feature)
        if id_ := feature.id:
            self.by_id[id_] = feature


class FeatureReader(Reader):
    def __init__(
        self, input_: TextIO | Path | str, /, streaming_window: int | None = 1000
    ):
        super().__init__(input_)
        self._streaming_window = streaming_window

    def _finalize_drafts(
        self, drafts: FeatureDrafts, w: int | None
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

        result = Features()
        for draft in drafts.drafts:
            feature = self._finalize_draft(draft, result)
            result.add(feature)
        yield from result.features

    def _finalize_draft(self, draft: FeatureDraft, result: Features) -> Feature:
        match draft.type_.lower():
            case "gene":
                feature = self._finalize_gene(draft, result)
            case "exon":
                feature = self._finalize_exon(draft, result)
            case "three_prime_utr":
                feature = self._finalize_three_prime_utr(draft, result)
            case _:
                feature = self._finalize_other(draft, result)
        if feature.parent:
            new_children = feature.parent.children + (feature,)
            object.__setattr__(feature.parent, "children", new_children)
        return feature

    def _finalize_gene(self, draft: FeatureDraft, result: Features) -> Feature:
        feature = self._finalize_other(draft, result)
        name = draft.pick_attribute("gene_name", "Name")
        biotype = draft.pick_attribute("gene_biotype", "biotype", "gene_type")
        if name is None or biotype is None:
            return feature
        return Gene(**feature.__dict__, name=name, biotype=biotype, exons=())

    def _finalize_exon(self, draft: FeatureDraft, result: Features) -> Feature:
        feature = self._finalize_other(draft, result)

        gene = feature.parent
        while gene and not isinstance(gene, Gene):
            gene = gene.parent

        if gene is None:
            return feature
        exon = Exon(**feature.__dict__, gene=gene)
        object.__setattr__(gene, "exons", gene.exons + (exon,))
        return exon

    def _finalize_three_prime_utr(
        self, draft: FeatureDraft, result: Features
    ) -> Feature:
        feature = self._finalize_other(draft, result)

        gene = feature.parent
        while gene and not isinstance(gene, Gene):
            gene = gene.parent

        if gene is None:
            return feature
        return ThreePrimeUTR(**feature.__dict__, gene=gene)

    def _finalize_other(self, draft: FeatureDraft, result: Features) -> Feature:
        parent_id = self._extract_parent_id(draft)
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
            id=self._extract_id(draft),
            parent=parent,
            children=(),
        )

    def _extract_id(self, draft: FeatureDraft) -> str | None:
        if (id_ := draft.attributes.get("ID")) is not None:
            return id_
        if draft.type_ == "gene" and (id_ := draft.attributes.get("gene_id")):
            return id_
        if draft.type_ == "transcript" and (
            id_ := draft.attributes.get("transcript_id")
        ):
            return id_
        return None

    def _extract_parent_id(self, draft: FeatureDraft) -> str | None:
        if (id_ := draft.attributes.get("Parent")) is not None:
            return id_
        if draft.type_ == "transcript" and (id_ := draft.attributes.get("gene_id")):
            return id_
        if draft.type_ == "exon" and (id_ := draft.attributes.get("transcript_id")):
            return id_
        return None
