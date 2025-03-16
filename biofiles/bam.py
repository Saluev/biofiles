import gzip
import struct
import sys
from io import BytesIO
from pathlib import Path
from types import TracebackType
from typing import Iterator, Any

from biofiles.types.alignment import (
    ReferenceSequence,
    Alignment,
    BAMTag,
    CIGAR,
    CIGAROpKind,
    CIGAROperation,
)


class BAMReader:
    def __init__(self, input_: BytesIO | Path | str) -> None:
        if isinstance(input_, Path | str):
            input_ = open(input_, "rb")
        self._input = input_
        self._ungzipped_input = gzip.open(input_)

        self._header_text: str | None = None
        self._ref_seqs: list[ReferenceSequence] = []

        self._read_header()

    def _read_header(self) -> None:
        magic_bytes = self._ungzipped_input.read(8)
        magic_data = struct.unpack("<ccccI", magic_bytes)
        if b"".join(magic_data[:4]) != b"BAM\1":
            raise ValueError("not a BAM file, invalid magic bytes")

        header_text_length = magic_data[-1]
        self._header_text = self._ungzipped_input.read(header_text_length)
        (num_ref_seqs,) = struct.unpack("<I", self._ungzipped_input.read(4))

        for _ in range(num_ref_seqs):
            (ref_seq_name_length,) = struct.unpack("<I", self._ungzipped_input.read(4))
            ref_seq_name = self._ungzipped_input.read(ref_seq_name_length)
            (ref_seq_length,) = struct.unpack("<I", self._ungzipped_input.read(4))
            ref_seq = ReferenceSequence(
                id=ref_seq_name.rstrip(b"\0").decode("ascii"), length=ref_seq_length
            )
            self._ref_seqs.append(ref_seq)

    def __iter__(self) -> Iterator[Alignment]:
        return self

    def __next__(self) -> Alignment:
        block_size_bytes = self._ungzipped_input.read(4)
        if not block_size_bytes:
            raise StopIteration

        (block_length,) = struct.unpack("<I", block_size_bytes)

        body_format = "<iiBBHHHIiii"
        body_bytes = self._ungzipped_input.read(struct.calcsize(body_format))
        (
            ref_seq_idx,
            pos,
            read_name_length,
            mapping_quality,
            bai_index_bin,
            num_cigar_ops,
            flags,
            seq_length,
            next_ref_seq_idx,
            next_pos,
            template_length,
        ) = struct.unpack(body_format, body_bytes)
        read_name_bytes = self._ungzipped_input.read(read_name_length)

        cigar_format = "<" + "I" * num_cigar_ops
        cigar_bytes = self._ungzipped_input.read(struct.calcsize(cigar_format))
        encoded_cigar = struct.unpack(cigar_format, cigar_bytes)

        seq_bytes = self._ungzipped_input.read((seq_length + 1) // 2)
        encoded_seq = struct.unpack("<" + "B" * len(seq_bytes), seq_bytes)

        quality = self._ungzipped_input.read(seq_length).decode("ascii")

        remaining_length = (
            block_length
            - len(body_bytes)
            - len(read_name_bytes)
            - len(cigar_bytes)
            - len(seq_bytes)
            - len(quality)
        )

        tags: list[BAMTag] = []
        while remaining_length > 0:
            tag, used_length = self._read_tag()
            tags.append(tag)
            remaining_length -= used_length
        if remaining_length < 0:
            raise ValueError("invalid BAM file, wrong tag length")

        ref_seq = self._ref_seqs[ref_seq_idx] if ref_seq_idx >= 0 else None
        next_ref_seq = (
            self._ref_seqs[next_ref_seq_idx] if next_ref_seq_idx >= 0 else None
        )
        return Alignment(
            reference_sequence=ref_seq,
            start_c=pos,
            read_name=read_name_bytes.rstrip(b"\0").decode("utf-8"),
            mapping_quality=mapping_quality,
            bai_index_bin=bai_index_bin,
            next_reference_sequence=next_ref_seq,
            next_start_c=next_pos,
            template_length=template_length,
            cigar=self._decode_cigar(encoded_cigar),
            read_sequence=self._decode_seq(encoded_seq),
            quality=quality,
            bam_flags=flags,
            bam_tags=tuple(tags),
        )

    def _decode_cigar(self, encoded_cigar: tuple[int, ...]) -> CIGAR:
        return CIGAR(
            operations=tuple(
                CIGAROperation(kind=_BAM_CIGAR_OP_KINDS[item & 0b1111], count=item >> 4)
                for item in encoded_cigar
            )
        )

    def _decode_seq(self, encoded_seq: tuple[int, ...]) -> str:
        return "".join(
            f"{_BAM_SEQUENCE_LETTERS[b >> 4]}{_BAM_SEQUENCE_LETTERS[b & 15]}"
            for b in encoded_seq
        )

    def _read_tag(self) -> tuple[BAMTag, int]:
        tag = self._ungzipped_input.read(2).decode("ascii")
        value_type = self._ungzipped_input.read(1)
        value, value_length = self._read_tag_value(value_type)
        return BAMTag(tag=tag, value=value), 3 + value_length

    def _read_tag_value(self, value_type: bytes) -> tuple[Any, int]:
        if value_type in (b"Z", b"H"):
            characters: list[bytes] = []
            last_character = b""
            while last_character != b"\0":
                characters.append(last_character)
                last_character = self._ungzipped_input.read(1)
            value = b"".join(characters).decode("utf-8")
            return value, len(characters)

        elif value_type == b"B":
            subtype, count = struct.unpack("<cI", self._ungzipped_input.read(5))
            format_ = "<" + _BAM_FORMAT_TO_STRUCT_FORMAT[subtype] * count
            length = struct.calcsize(format_)
            value = struct.unpack(format_, self._ungzipped_input.read(length))
            return value, 5 + length

        else:
            format_ = "<" + _BAM_FORMAT_TO_STRUCT_FORMAT[value_type]
            length = struct.calcsize(format_)
            (value,) = struct.unpack(format_, self._ungzipped_input.read(length))
            return value, length

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


_BAM_FORMAT_TO_STRUCT_FORMAT = {
    b"A": "c",
    b"c": "b",
    b"C": "B",
    b"s": "h",
    b"S": "H",
    b"i": "i",
    b"I": "I",
    b"f": "f",
}

_BAM_CIGAR_OP_KINDS: list[CIGAROpKind] = ["M", "I", "D", "N", "S", "H", "P", "=", "X"]
_BAM_SEQUENCE_LETTERS = "=ACMGRSVTWYHKDBN"

if __name__ == "__main__":
    for path in sys.argv[1:]:
        num_alignments = 0
        with BAMReader(path) as reader:
            for record in reader:
                num_alignments += 1
        print(f"Parsed {num_alignments} alignments from {path}")
