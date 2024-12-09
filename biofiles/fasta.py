from dataclasses import dataclass, field
from pathlib import Path
from types import TracebackType
from typing import TextIO, Iterator, Self


__all__ = ["FASTASequence", "FASTAReader", "FASTAWriter"]


@dataclass(frozen=True)
class FASTASequence:
    id: str
    description: str
    sequence: str


@dataclass
class _FASTASequenceDraft:
    id: str
    description: str
    sequence_parts: list[str] = field(default_factory=list)

    def finalize(self) -> FASTASequence:
        return FASTASequence(
            id=self.id,
            description=self.description,
            sequence="".join(self.sequence_parts),
        )


class FASTAReader:
    def __init__(self, input_: TextIO | Path | str) -> None:
        if isinstance(input_, Path | str):
            input_ = open(input_)
        self._input = input_

    def __enter__(self) -> "FASTAReader":
        self._input.__enter__()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> None:
        self._input.__exit__(exc_type, exc_val, exc_tb)

    def __iter__(self) -> Iterator[FASTASequence]:
        draft: _FASTASequenceDraft | None = None
        for line in self._input:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if draft:
                    yield draft.finalize()
                line = line.removeprefix(">").lstrip()
                match line.split(maxsplit=1):
                    case [id_, desc]:
                        pass
                    case [id_]:
                        desc = ""
                    case []:
                        raise ValueError(
                            f"unexpected line {line!r}, expected a non-empty sequence identifier"
                        )
                draft = _FASTASequenceDraft(id=id_, description=desc)
            elif line:
                if not draft:
                    raise ValueError(f"unexpected line {line!r}, expected >")
                draft.sequence_parts.append(line)
        if draft:
            yield draft.finalize()


class FASTAWriter:
    def __init__(self, output: TextIO | Path | str, width: int = 80) -> None:
        if isinstance(output, Path | str):
            output = open(output, "w")
        self._output = output
        self._width = width

    def __enter__(self) -> "FASTAWriter":
        self._output.__enter__()
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_val: BaseException | None,
        exc_tb: TracebackType | None,
    ) -> None:
        self._output.__exit__(exc_type, exc_val, exc_tb)

    def write(self, sequence: FASTASequence) -> None:
        self._output.write(f">{sequence.id} {sequence.description}\n")
        sequence_len = len(sequence.sequence)
        for offset in range(0, sequence_len, self._width):
            self._output.write(
                sequence.sequence[offset : min(offset + self._width, sequence_len)]
            )
            self._output.write("\n")
