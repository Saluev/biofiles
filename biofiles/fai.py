import sys
from typing import Iterator

from biofiles.common import Reader
from biofiles.types.sequence import SequenceDescription


class FAIReader(Reader):
    def __iter__(self) -> Iterator[SequenceDescription]:
        for line in self._input:
            sequence_id, length_str, byte_offset_str, _, _ = line.rstrip("\n").split(
                "\t"
            )
            yield SequenceDescription(
                id=sequence_id, length=int(length_str), byte_offset=int(byte_offset_str)
            )


if __name__ == "__main__":
    for path in sys.argv[1:]:
        with FAIReader(path) as reader:
            for seq_desc in reader:
                print(seq_desc)
