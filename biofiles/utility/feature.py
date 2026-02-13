from collections import deque, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator, TextIO, Type

from biofiles.common import Strand, Reader
from biofiles.types.feature import (
    Feature,
    FeatureMetaclass,
    Relation,
    get_composite_field,
    Dialect,
)


@dataclass(slots=True)
class FeatureDraft:
    idx: int
    sequence_id: str
    source: str
    type_: str
    start_original: int
    end_original: int
    start_c: int
    end_c: int
    score: float | None
    strand: Strand | None
    phase: int | None
    attributes: dict[str, str | list[str]]

    class_: Type[Feature] | None = None
    id: Any = None
    finalized: Feature | None = None


class FeatureTypes:
    ambiguous_type_mapping: dict[str, list[FeatureMetaclass]]
    unique_type_mapping: dict[str, FeatureMetaclass]

    def __init__(self, feature_types: list[FeatureMetaclass]) -> None:
        self.ambiguous_type_mapping = defaultdict(list)
        self.unique_type_mapping = {}

        for ft in feature_types:
            for type in ft.__filter_type__:
                self.ambiguous_type_mapping[type].append(ft)

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


@dataclass(slots=True)
class FeatureDrafts:
    feature_types: FeatureTypes
    drafts: list[FeatureDraft] = field(default_factory=deque)
    by_class_and_id: dict[tuple[type, Any], FeatureDraft] = field(default_factory=dict)

    def add(self, draft: FeatureDraft) -> None:
        self.drafts.append(draft)
        if class_ := self.feature_types.unique_type_mapping.get(draft.type_.lower()):
            draft.class_ = class_
            draft.id = get_composite_field(
                draft.attributes, class_.__id_attribute_source__
            )
            self.register(draft)

    def register(self, draft: FeatureDraft) -> None:
        if draft.id is None:
            return
        if (key := (draft.class_, draft.id)) in self.by_class_and_id:
            raise ValueError(
                f"duplicate feature ID {draft.id} for class {draft.class_.__name__}"
            )
        self.by_class_and_id[key] = draft


class RawFeatureReader(Reader):
    def __init__(self, input_: TextIO | Path) -> None:
        super().__init__(input_)

    def __iter__(self) -> Iterator[FeatureDraft]:
        raise NotImplementedError


class FeatureReader(Reader):

    def __init__(self, input_: TextIO | Path | str, dialect: Dialect) -> None:
        super().__init__(input_)
        self._feature_types = FeatureTypes(dialect.feature_types)
        self._raw_reader = self._make_raw_feature_reader()

    def _make_raw_feature_reader(self) -> RawFeatureReader:
        raise NotImplementedError

    def __iter__(self) -> Iterator[Feature]:
        fds = FeatureDrafts(self._feature_types)
        for draft in self._raw_reader:
            fds.add(draft)
        yield from self._finalize_drafts(fds)

    def _finalize_drafts(self, fds: FeatureDrafts) -> Iterator[Feature]:
        self._choose_classes(fds)
        self._instantiate_objects(fds)
        self._fill_relations(fds)
        for fd in fds.drafts:
            yield fd.finalized

    def _choose_classes(self, fds: FeatureDrafts) -> None:
        for fd in fds.drafts:
            if fd.class_:
                continue

            fts = self._feature_types.ambiguous_type_mapping[fd.type_]
            matching_fts = [ft for ft in fts if self._check_filters(fds, fd, ft)]
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
            fd.id = get_composite_field(fd.attributes, ft.__id_attribute_source__)
            fds.register(fd)

    def _instantiate_objects(self, fds: FeatureDrafts) -> None:
        for fd in fds.drafts:
            fd.finalized = fd.class_(
                sequence_id=fd.sequence_id,
                source=fd.source,
                type_=fd.type_,
                start_original=fd.start_original,
                end_original=fd.end_original,
                start_c=fd.start_c,
                end_c=fd.end_c,
                score=fd.score,
                strand=fd.strand,
                phase=fd.phase,
                attributes=fd.attributes,
            )

    def _fill_relations(self, fds: FeatureDrafts) -> None:
        for fd in fds.drafts:
            for relation in fd.class_.__relations__:
                related_id = get_composite_field(
                    fd.attributes, relation.id_attribute_source
                )
                related_class = relation.inverse.class_
                try:
                    related_fd = fds.by_class_and_id[related_class, related_id]
                except KeyError as exc:
                    raise ValueError(
                        f"can't find related {related_class.__name__} {related_id} for {fd.finalized}"
                    ) from exc
                setattr(fd.finalized, relation.attribute_name, related_fd.finalized)
                if relation.inverse.attribute_name is None:
                    pass
                elif relation.inverse.one_to_one:
                    if (
                        getattr(
                            related_fd.finalized, relation.inverse.attribute_name, None
                        )
                        is not None
                    ):
                        raise ValueError(
                            f"too many related {relation.class_.__name__}s for {related_fd.finalized} (at least 2), "
                            f"expected one-to-one relation"
                        )
                    setattr(
                        related_fd.finalized,
                        relation.inverse.attribute_name,
                        fd.finalized,
                    )
                else:
                    getattr(
                        related_fd.finalized, relation.inverse.attribute_name
                    ).append(fd.finalized)

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
        related_id = fd.attributes[r.id_attribute_source]
        try:
            return fds.by_class_and_id[related_class, related_id]
        except KeyError as exc:
            raise ValueError(
                f"can't find related {related_class.__name__} for "
                f"{fd.class_.__name__} with attributes {fd.attributes!r}"
            ) from exc
