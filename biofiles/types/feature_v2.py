from dataclasses import dataclass, Field, field as dataclass_field
from typing import dataclass_transform, Type, Any, TypeAlias

from biofiles.common import Strand

Source: TypeAlias = str | tuple[str, ...]


@dataclass
class Relation:
    """Equivalent of SQL foreign key — declarative description
    of a relation between two types of features."""

    id_attribute_source: Source
    """ Name of GTF/GFF attribute(s) which contains related feature ID. """

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
    __id_attribute_source__: Source
    """ Name of GTF/GFF attribute(s) which contains the type-unique ID. """

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
        result.__id_attribute_source__ = cls._find_id_attribute_source(namespace)
        result._fill_relation_classes(namespace)
        result._fill_filters(type=type, starts=starts, ends=ends)
        result._fill_slots()
        result._fill_init_method(namespace)

        # TODO generate dataclass-like __init__ method,
        #      keep all relations optional

        return result

    @staticmethod
    def _find_id_attribute_source(namespace) -> str:
        result = ""
        for key, value in namespace.items():
            match value:
                case Field(metadata={"id_attribute_name": id_attribute_source}):
                    if result:
                        raise TypeError(
                            f"should specify exactly one id_field() in class {result.__name__}"
                        )
                    result = id_attribute_source
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

    def _fill_slots(cls) -> None:
        cls.__slots__ = [
            key
            for ancestor in cls.__mro__[::-1][1:]
            for key in ancestor.__annotations__
        ]

    def _fill_init_method(cls, namespace) -> None:
        default_arguments: list[str] = []
        non_default_arguments: list[str] = []
        assignments: list[str] = []

        key_to_ancestor: dict[str, Type] = {}
        for ancestor in cls.__mro__[:-1]:
            for key, value in ancestor.__annotations__.items():
                key_to_ancestor.setdefault(key, ancestor)

        for ancestor in cls.__mro__[::-1][1:]:
            for key, value in ancestor.__annotations__.items():
                if key_to_ancestor[key] is not ancestor:
                    # Overridden in a descendant class.
                    continue

                field_value = getattr(cls, key, None)
                argument, assignment = cls._compose_field(key, value, field_value)

                if argument and argument.endswith(" = None"):
                    default_arguments.append(argument)
                elif argument:
                    non_default_arguments.append(argument)
                assignments.append(assignment)

        body = "\n    ".join(assignments)
        all_arguments = [*non_default_arguments, *default_arguments]
        source_code = f"def __init__(self, {', '.join(all_arguments)}):\n    {body}"
        locals = {}
        exec(source_code, {}, locals)
        cls.__init__ = locals["__init__"]

    def _compose_field(
        cls, field_name: str, field_annotation: Any, field_value: Field | None
    ) -> tuple[str | None, str]:
        argument: str | None
        assignment: str
        match field_value:
            case Field(metadata={"relation": _}):
                argument = f"{field_name}: {cls._format_type_arg(field_annotation, optional=True)} = None"
                assignment = f"self.{field_name} = {field_name}"
            case Field(metadata={"attribute_name": attribute_name}) | Field(
                metadata={"id_attribute_name": attribute_name}
            ):
                argument = None
                assignment = f"self.{field_name} = attributes[{repr(attribute_name)}]"
                # TODO necessary conversions, proper exceptions
            case None:
                argument = f"{field_name}: {cls._format_type_arg(field_annotation, optional=False)}"
                assignment = f"self.{field_name} = {field_name}"
            case other:
                raise TypeError(f"unsupported field: {field_value}")
        return argument, assignment

    def _format_type_arg(cls, type: str | Type, optional: bool) -> str:
        if isinstance(type, str):
            return f'"{type} | None"' if optional else type
        try:
            if type.__module__ == "builtins":
                return f"{type.__name__} | None" if optional else type.__name__
            return f'"{type.__module__}.{type.__name__}"'
        except AttributeError:
            # TODO Properly support Optional, Union, etc., especially with built-in types
            return f'"{str(type)} | None"' if optional else repr(str(type))


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

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.sequence_id}:{self.start_c}-{self.end_c})"


def id_field(source: Source) -> Field:
    return dataclass_field(metadata={"id_attribute_name": source})


def field(source: Source) -> Field:
    return dataclass_field(metadata={"attribute_name": source})


def relation(source: Source, *, one_to_one: bool = False) -> tuple[Field, Field]:
    forward_relation = Relation(id_attribute_source=source)
    inverse_relation = InverseRelation(inverse=forward_relation, one_to_one=one_to_one)
    forward_relation.inverse = inverse_relation

    forward_field = dataclass_field(metadata={"relation": forward_relation})
    inverse_field = dataclass_field(metadata={"relation": inverse_relation})
    return forward_field, inverse_field
