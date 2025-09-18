from dataclasses import dataclass, Field, field as dataclass_field
from typing import dataclass_transform

from biofiles.common import Strand


@dataclass
class Relation:
    id_field_name: str
    inverse: "InverseRelation | None" = None
    class_: type | None = None


@dataclass
class InverseRelation:
    inverse: Relation
    one_to_one: bool
    class_: type | None = None


@dataclass_transform()
class FeatureMetaclass(type):
    __id_field_name__: str
    __filter_type__: str
    __filter_starts__: Relation | None
    __filter_ends__: Relation | None

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
        result.__id_field_name__ = ""
        for key, value in namespace.items():
            match value:
                case Field(metadata={"id_field_name": id_field_name}):
                    if result.__id_field_name__:
                        raise TypeError(
                            f"should specify exactly one id_field() in class {result.__name__}"
                        )
                    result.__id_field_name__ = id_field_name
                case Field(metadata={"relation": Relation() as r}):
                    r.class_ = result
                    if key in result.__annotations__:
                        # TODO handle optionality and forward refs
                        r.inverse.class_ = result.__annotations__[key]
                case Field(metadata={"relation": InverseRelation() as r}):
                    r.class_ = result
                    # TODO calculating r.inverse.class_ based on type annotation

        if type is not None:
            result.__filter_type__ = type
        result.__filter_starts__ = None
        if starts is not None:
            result.__filter_starts__ = starts.metadata["relation"]
        result.__filter_ends__ = None
        if ends is not None:
            result.__filter_ends__ = ends.metadata["relation"]

        # TODO generate dataclass-like __init__ method,
        #      keep all relations optional

        return result


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
    return dataclass_field(metadata={"id_field_name": source})


def field(source: str) -> Field:
    return dataclass_field(metadata={"field_name": source})


def relation(source: str, *, one_to_one: bool = False) -> tuple[Field, Field]:
    forward = Relation(id_field_name=source)
    inverse = InverseRelation(inverse=forward, one_to_one=one_to_one)
    forward.inverse = inverse

    return dataclass_field(metadata={"relation": forward}), dataclass_field(
        metadata={"relation": inverse}
    )
