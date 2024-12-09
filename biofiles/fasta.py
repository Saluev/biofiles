from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO, Iterator

from biofiles.common import Reader, Writer
from biofiles.types.sequence import Sequence


__all__ = ["FASTAReader", "FASTAWriter"]


@dataclass
class _SequenceDraft:
    id: str
    description: str
    sequence_parts: list[str] = field(default_factory=list)

    def finalize(self) -> Sequence:
        return Sequence(
            id=self.id,
            description=self.description,
            sequence="".join(self.sequence_parts),
        )


class FASTAReader(Reader):
    def __iter__(self) -> Iterator[Sequence]:
        draft: _SequenceDraft | None = None
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
                draft = _SequenceDraft(id=id_, description=desc)
            elif line:
                if not draft:
                    raise ValueError(f"unexpected line {line!r}, expected >")
                draft.sequence_parts.append(line)
        if draft:
            yield draft.finalize()


class FASTAWriter(Writer):
    def __init__(self, output: TextIO | Path | str, width: int = 80) -> None:
        super().__init__(output)
        self._width = width

    def write(self, sequence: Sequence) -> None:
        self._output.write(f">{sequence.id} {sequence.description}\n")
        sequence_len = len(sequence.sequence)
        for offset in range(0, sequence_len, self._width):
            self._output.write(
                sequence.sequence[offset : min(offset + self._width, sequence_len)]
            )
            self._output.write("\n")
