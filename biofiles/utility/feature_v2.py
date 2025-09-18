from collections import deque, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator, TextIO

from biofiles.common import Strand, Reader
from biofiles.types.feature_v2 import Feature, FeatureMetaclass, Relation


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

    class_: type | None = None
    id: Any = None
    finalized: Feature | None = None


class FeatureTypes:
    def __init__(self, feature_types: list[FeatureMetaclass]) -> None:
        for ft in feature_types:
            if not ft.__id_field_name__:
                raise ValueError(
                    f"{ft.__name__} is not proper feature type - has no id_field()"
                )

        self.ambiguous_type_mapping: dict[str, list[FeatureMetaclass]] = defaultdict(
            list
        )
        self.unique_type_mapping: dict[str, FeatureMetaclass] = {}

        for ft in feature_types:
            self.ambiguous_type_mapping[ft.__filter_type__].append(ft)

        for key, fts in [*self.ambiguous_type_mapping.items()]:
            if len(fts) == 1:
                self.unique_type_mapping[key] = fts[0]
                del self.ambiguous_type_mapping[key]
                continue
            self.ambiguous_type_mapping[key] = _sort_by_filter_specificity(fts)


def _sort_by_filter_specificity(fts: list[FeatureMetaclass]) -> list[FeatureMetaclass]:
    """Sort feature classes by their filter specificity, most specific -> least specific."""
    key = lambda ft: bool(ft.__filter_starts__) + bool(ft.__filter_ends__)
    return sorted(fts, key=key, reverse=True)


@dataclass
class FeatureDrafts:
    feature_types: FeatureTypes
    drafts: list[FeatureDraft] = field(default_factory=deque)
    by_class_and_id: dict[tuple[type, Any], FeatureDraft] = field(default_factory=dict)

    def add(self, draft: FeatureDraft) -> None:
        self.drafts.append(draft)
        if class_ := self.feature_types.unique_type_mapping.get(draft.type_):
            draft.class_ = class_
            draft.id = draft.attributes[class_.__id_field_name__]
            self.register(draft)

    def register(self, draft: FeatureDraft) -> None:
        if (key := (draft.class_, draft.id)) in self.by_class_and_id:
            raise ValueError(
                f"duplicate feature ID {draft.id} for class {class_.__name__}"
            )
        self.by_class_and_id[key] = draft


class FeatureReader(Reader):

    def __init__(
        self, input_: TextIO | Path | str, feature_types: list[FeatureMetaclass]
    ) -> None:
        super().__init__(input_)
        self._feature_types = FeatureTypes(feature_types)

    def __iter__(self) -> Iterator[Feature]:
        raise NotImplementedError

    def _finalize_drafts(self, fds: FeatureDrafts) -> Iterator[Feature]:
        self._choose_classes(fds)
        pass

    def _choose_classes(self, fds: FeatureDrafts) -> Iterator[Feature]:
        for fd in fds.drafts:
            if fd.class_:
                continue

            fts = self._feature_types.ambiguous_type_mapping[fd.type_]
            matching_fts = [ft for ft in fts if self._check_filters(fd, ft)]
            if not matching_fts:
                raise ValueError(
                    f"no matching classes (out of {len(fts)}) for "
                    f"feature with type {fd.type_!r}, attributes {fd.attributes!r}"
                )
            if len(matching_fts) > 1:
                raise ValueError(
                    f"too many matching classes ({len(matching_fts)}) for "
                    f"feature with type {fd.type_!r}, attributes {fd.attributes!r}"
                )
            ft = matching_fts[0]
            fd.class_ = ft
            fd.id = fd.attributes[ft.__id_field_name__]
            fds.register(fd)

    def _check_filters(
        self, fds: FeatureDrafts, fd: FeatureDraft, ft: FeatureMetaclass
    ) -> bool:
        if r := ft.__filter_starts__:
            related_fd = self._get_related_feature_draft(fds, fd, r)
            if fd.strand != related_fd.strand:
                return False
            if fd.strand == "+" and fd.start_original != related_fd.start_original:
                return False
            if fd.strand == "-" and fd.end_original != related_fd.end_original:
                return False
        if r := ft.__filter_ends__:
            related_fd = self._get_related_feature_draft(fds, fd, r)
            if fd.strand != related_fd.strand:
                return False
            if fd.strand == "+" and fd.end_original != related_fd.end_original:
                return False
            if fd.strand == "-" and fd.start_original != related_fd.start_original:
                return False
        return True

    def _get_related_feature_draft(
        self, fds: FeatureDrafts, fd: FeatureDraft, r: Relation
    ) -> FeatureDraft:
        related_class = r.inverse.class_
        related_id = fd.attributes[r.id_field_name]
        try:
            return fds.by_class_and_id[related_class, related_id]
        except KeyError as exc:
            raise ValueError(
                f"can't find related {related_class.__name__} for "
                f"{fd.class_.__name__} with attributes {fd.attributes!r}"
            ) from exc
