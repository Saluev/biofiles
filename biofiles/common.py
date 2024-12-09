from pathlib import Path
from types import TracebackType
from typing import TypeAlias, Literal, TextIO

Strand: TypeAlias = Literal["+", "-"]


class Reader:
    def __init__(self, input_: TextIO | Path | str) -> None:
        if isinstance(input_, Path | str):
            input_ = open(input_)
        self._input = input_

    def __enter__(self):
        self._input.__enter__()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> None:
        self._input.__exit__(exc_type, exc_val, exc_tb)


class Writer:
    def __init__(self, output: TextIO | Path | str) -> None:
        if isinstance(output, Path | str):
            output = open(output, "w")
        self._output = output

    def __enter__(self):
        self._output.__enter__()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> None:
        self._output.__exit__(exc_type, exc_val, exc_tb)
