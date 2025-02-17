from dataclasses import dataclass
from pathlib import Path
from typing import TypeAlias, Callable, Any, Literal, Type

from biofiles.types.feature import Feature, Gene, Transcript, UTR, Exon

FeatureFilter: TypeAlias = Callable[[Feature], bool]
FeatureMapper: TypeAlias = Callable[[Feature], Any]


@dataclass
class Pipeline:
    inputs: list[Path]
    filters: list[FeatureFilter]
    mapper: FeatureMapper | None

    def filter(self, feature: Feature) -> bool:
        for f in self.filters:
            if not f(feature):
                return False
        return True

    def map(self, feature: Feature) -> Any:
        if not self.mapper:
            return feature
        return self.mapper(feature)


Mode: TypeAlias = Literal["inputs", "filters", "done"]


def parse_pipeline_args(argv: list[str]) -> Pipeline:
    pipeline = Pipeline(inputs=[], filters=[], mapper=None)

    mode: Mode = "inputs"
    i = 0
    while i < len(argv):
        match mode, argv[i:]:
            case "inputs", [str_path, *_] if (path := Path(str_path)).is_file():
                pipeline.inputs.append(path)
                i += 1
            case "inputs", ["--filter", *_]:
                mode = "filters"
                i += 1
            case "inputs" | "filters", ["--attr", key]:
                path = key.split(".")
                pipeline.mapper = _produce_attr_mapper(path)
                mode = "done"
                i += 2
            case "filters", [filter_str, *_]:
                filter_ = _parse_filter(filter_str)
                pipeline.filters.append(filter_)
                i += 1
            case other:
                raise ValueError(f"can't parse command line arguments {argv[i:]}")

    return pipeline


def _parse_filter(filter_str: str) -> FeatureFilter:
    if "=" not in filter_str:
        # --filter gene,transcript
        type_strs = filter_str.split(",")
        types = tuple(_parse_feature_type(t) for t in type_strs)
        return lambda f: isinstance(f, types)

    # --filter attr=value1,value2
    key, value = filter_str.split("=", maxsplit=1)
    values = value.split(",")
    match key:
        case "chromosome":
            return lambda f: f.sequence_id in values
        case "type":
            return lambda f: f.type_ in values
        case "strand":
            return lambda f: f.strand in values
        case _:
            path = key.split(".")
            return _produce_attr_filter(path, values)

    raise ValueError(f"can't parse filter {filter_str!r}")


def _parse_feature_type(t: str) -> Type[Feature]:
    if t not in _FEATURE_TYPES:
        raise ValueError(f"unknown feature type {t!r}")
    return _FEATURE_TYPES[t]


def _produce_attr_filter(path: list[str], values: list[str]) -> FeatureFilter:
    assert path
    if len(path) == 1:
        (key,) = path
        match key:
            case "chromosome" | "type" | "strand" | "id":
                return lambda f: getattr(f, key) in values
            # TODO other attributes
            case _:
                return lambda f: f.attributes.get(key) in values

    if path[0] not in ("gene", "transcript", "parent"):
        raise ValueError(f"unknown attribute {path[-2]!r}")

    nested = _produce_attr_filter(path[1:], values)
    return lambda f: (nested(nf) if (nf := getattr(f, path[0], None)) else False)


def _produce_attr_mapper(path: list[str]) -> FeatureMapper:
    assert path
    if len(path) == 1:
        (key,) = path
        match key:
            case "chromosome" | "type" | "strand" | "id":
                return lambda f: getattr(f, key)
            # TODO other attributes
            case _:
                return lambda f: f.attributes.get(key, "")

    if path[0] not in ("gene", "transcript", "parent"):
        raise ValueError(f"unknown attribute {path[-2]!r}")

    nested = _produce_attr_mapper(path[1:])
    return lambda f: (nested(nf) if (nf := getattr(f, path[0], None)) else None)


_FEATURE_TYPES = {"gene": Gene, "transcript": Transcript, "exon": Exon, "utr": UTR}
