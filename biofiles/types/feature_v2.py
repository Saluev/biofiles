from dataclasses import dataclass, Field, field as dataclass_field
from typing import dataclass_transform

from biofiles.common import Strand


@dataclass
class Relation:
    """Equivalent of SQL foreign key — declarative description
    of a relation between two types of features."""

    id_attribute_name: str
    """ Name of GTF/GFF attribute which contains related feature ID. """

    inverse: "InverseRelation | None" = None

    class_: type | None = None
    """ Python class for the related feature. """


@dataclass
class InverseRelation:
    inverse: Relation
    one_to_one: bool
    class_: type | None = None


@dataclass_transform()
class FeatureMetaclass(type):
    __id_attribute_name__: str
    """ Name of GTF/GFF attribute which contains the type-unique ID. """

    __filter_type__: str
    """ Filter by feature type ("gene", "transcript", etc.). """

    __filter_starts__: Relation | None
    """ Filter by start position — feature starts at the same point as related feature. """

    __filter_ends__: Relation | None
    """ Filter by end position — feature ends at the same point as related feature. """

    def __new__(
        cls,
        name,
        bases,
        namespace,
        type: str | None = None,
        starts: Field | None = None,
        ends: Field | None = None,
    ):
        result = super().__new__(cls, name, bases, namespace)
        result.__id_attribute_name__ = cls._find_id_attribute(namespace)
        result._fill_relation_classes(namespace)
        result._fill_filters(type=type, starts=starts, ends=ends)

        # TODO generate dataclass-like __init__ method,
        #      keep all relations optional

        return result

    @staticmethod
    def _find_id_attribute(namespace) -> str:
        result = ""
        for key, value in namespace.items():
            match value:
                case Field(metadata={"id_attribute_name": id_attribute_name}):
                    if result:
                        raise TypeError(
                            f"should specify exactly one id_field() in class {result.__name__}"
                        )
                    result = id_attribute_name
        return result

    def _fill_relation_classes(cls, namespace) -> None:
        for key, value in namespace.items():
            match value:
                case Field(metadata={"relation": Relation() as r}):
                    r.class_ = cls
                    if key in cls.__annotations__:
                        # TODO handle optionality and forward refs
                        r.inverse.class_ = cls.__annotations__[key]
                case Field(metadata={"relation": InverseRelation() as r}):
                    r.class_ = cls
                    # TODO calculating r.inverse.class_ based on type annotation

    def _fill_filters(
        cls,
        *,
        type: str | None = None,
        starts: Field | None = None,
        ends: Field | None = None,
    ) -> None:
        if type is not None:
            cls.__filter_type__ = type
        cls.__filter_starts__ = None
        if starts is not None:
            cls.__filter_starts__ = starts.metadata["relation"]
        cls.__filter_ends__ = None
        if ends is not None:
            cls.__filter_ends__ = ends.metadata["relation"]


class Feature(metaclass=FeatureMetaclass):
    sequence_id: str
    source: str
    type_: str

    start_original: int
    end_original: int
    # Original values as they were present in the file (1-based inclusive for .gff and .gtf).

    start_c: int
    end_c: int
    # Standardized ("C-style") 0-based values, start inclusive, end exclusive.

    score: float | None
    strand: Strand | None
    phase: int | None
    attributes: dict[str, str]


def id_field(source: str) -> Field:
    return dataclass_field(metadata={"id_attribute_name": source})


def field(source: str) -> Field:
    return dataclass_field(metadata={"attribute_name": source})


def relation(source: str, *, one_to_one: bool = False) -> tuple[Field, Field]:
    forward = Relation(id_attribute_name=source)
    inverse = InverseRelation(inverse=forward, one_to_one=one_to_one)
    forward.inverse = inverse

    return dataclass_field(metadata={"relation": forward}), dataclass_field(
        metadata={"relation": inverse}
    )
