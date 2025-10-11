from dataclasses import dataclass, Field, field as dataclass_field
from enum import Enum
from typing import dataclass_transform, Type, Any, TypeAlias
from uuid import uuid4

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
    attribute_name: str | None = None


@dataclass
class InverseRelation:
    inverse: Relation
    one_to_one: bool
    class_: type | None = None
    attribute_name: str | None = None


def get_composite_field(
    attributes: dict[str, str], source: Source
) -> str | tuple[str, ...] | None:
    if source is None:
        return None
    if isinstance(source, str):
        return attributes[source]
    return tuple(attributes[attribute_name] for attribute_name in source)


@dataclass_transform()
class FeatureMetaclass(type):
    __id_attribute_source__: Source | None
    """ Name of GTF/GFF attribute(s) which contains the type-unique ID. """

    __filter_type__: tuple[str, ...]
    """ Filter by feature type ("gene", "transcript", etc.). """

    __filter_starts__: Relation | None
    """ Filter by start position — feature starts at the same point as related feature. """

    __filter_ends__: Relation | None
    """ Filter by end position — feature ends at the same point as related feature. """

    __relations__: list[Relation]
    """ All direct relations for this type, for faster parsing. """

    def __new__(
        cls,
        name,
        bases,
        namespace,
        type: str | tuple[str, ...] | None = None,
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
        result: str | None = None
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
        cls.__relations__ = []
        for key, value in namespace.items():
            match value:
                case Field(metadata={"relation": Relation() as r}):
                    r.class_ = cls
                    r.attribute_name = key
                    if key in cls.__annotations__:
                        # TODO handle optionality and forward refs
                        r.inverse.class_ = cls.__annotations__[key]
                    cls.__relations__.append(r)
                case Field(metadata={"relation": InverseRelation() as r}):
                    r.class_ = cls
                    r.attribute_name = key
                    # TODO calculating r.inverse.class_ based on type annotation

    def _fill_filters(
        cls,
        *,
        type: str | tuple[str, ...] | None = None,
        starts: Field | None = None,
        ends: Field | None = None,
    ) -> None:
        if type is not None:
            cls.__filter_type__ = (type,) if isinstance(type, str) else type

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
        globals: dict[str, Any] = {}

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
                argument, assignment = cls._compose_field(
                    key, value, field_value, globals
                )

                if argument and argument.endswith(" = None"):
                    default_arguments.append(argument)
                elif argument:
                    non_default_arguments.append(argument)
                assignments.append(assignment)

        body = "\n    ".join(assignments)
        all_arguments = [*non_default_arguments, *default_arguments]
        source_code = f"def __init__(self, {', '.join(all_arguments)}):\n    {body}"
        locals = {}
        exec(source_code, globals, locals)
        cls.__init__ = locals["__init__"]

    def _compose_field(
        cls,
        field_name: str,
        field_annotation: Any,
        field_value: Field | None,
        globals: dict[str, Any],
    ) -> tuple[str | None, str]:
        argument: str | None
        assignment: str
        match field_value:
            case Field(metadata={"relation": r}):
                argument = f"{field_name}: {cls._format_type_arg(field_annotation, optional=True)} = None"
                if isinstance(r, InverseRelation) and not r.one_to_one:
                    assignment = f"self.{field_name} = {field_name} if {field_name} is not None else []"
                else:
                    assignment = f"self.{field_name} = {field_name}"
            case Field(metadata={"id_attribute_name": None}):
                argument = None
                assignment = f"self.{field_name} = None"
            case Field(metadata={"attribute_name": attribute_name}) | Field(
                metadata={"id_attribute_name": attribute_name}
            ):
                default = field_value.metadata.get("attribute_default", _no_default)
                default_factory = field_value.metadata.get(
                    "attribute_default_factory", _no_default
                )
                default_variable_name = f"default_{uuid4().hex}"
                argument = None
                if isinstance(attribute_name, str):
                    if default is not _no_default:
                        globals[default_variable_name] = default
                        getter = f"attributes.get({repr(attribute_name)}, {default_variable_name})"
                    elif default_factory is not _no_default:
                        globals[default_variable_name] = default_factory
                        getter = f"attributes.get({repr(attribute_name)}, {default_variable_name}())"
                    else:
                        getter = f"attributes[{repr(attribute_name)}]"
                else:
                    if default is not _no_default or default_factory is not _no_default:
                        raise NotImplementedError()
                    globals["get_composite_field"] = get_composite_field
                    getter = f"get_composite_field(attributes, {repr(attribute_name)})"
                if isinstance(field_annotation, type) and issubclass(
                    field_annotation, (int, float)
                ):
                    getter = f"{field_annotation.__name__}({getter})"
                elif isinstance(field_annotation, type) and issubclass(
                    field_annotation, Enum
                ):
                    globals[field_annotation.__name__] = field_annotation
                    getter = f"{field_annotation.__name__}({getter})"
                # TODO int | None, list[Enum], etc.
                # TODO ensure it's a list if annotated as list
                assignment = f"self.{field_name} = {getter}"
                # TODO necessary conversions, proper exceptions
            case None:
                argument = f"{field_name}: {cls._format_type_arg(field_annotation, optional=False)}"
                assignment = f"self.{field_name} = {field_name}"
            case property():
                argument = None
                assignment = ""
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


def no_id_field() -> Field:
    return dataclass_field(metadata={"id_attribute_name": None})


_no_default = object()


def field(
    source: Source, *, default: Any = _no_default, default_factory: Any = _no_default
) -> Field:
    metadata = {"attribute_name": source}
    if default is not _no_default:
        metadata["attribute_default"] = default
    if default_factory is not _no_default:
        metadata["attribute_default_factory"] = default_factory
    return dataclass_field(metadata=metadata)


def relation(source: Source, *, one_to_one: bool = False) -> tuple[Field, Field]:
    forward_relation = Relation(id_attribute_source=source)
    inverse_relation = InverseRelation(inverse=forward_relation, one_to_one=one_to_one)
    forward_relation.inverse = inverse_relation

    forward_field = dataclass_field(metadata={"relation": forward_relation})
    inverse_field = dataclass_field(metadata={"relation": inverse_relation})
    return forward_field, inverse_field


@dataclass(frozen=True)
class Dialect:
    name: str
    feature_types: list[Type[Feature]]
