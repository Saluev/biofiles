import sys
from collections import Counter
from typing import Iterator

from biofiles.common import Reader
from biofiles.types.repeat import Repeat


__all__ = ["RepeatMaskerReader"]


class RepeatMaskerReader(Reader):
    def __iter__(self) -> Iterator[Repeat]:
        for line in self._input:
            parts = line.split("\t")
            if not (14 <= len(parts) <= 15):
                # Probably some metainfo. No way to tell.
                continue

            (
                sw_score_str,
                div_str,
                del_str,
                ins_str,
                seq_id,
                seq_start_str,
                seq_end_str,
                seq_left_str,
                strand_str,
                repeat_name,
                repeat_class_family,
                repeat_start_str,
                repeat_end_str,
                repeat_left_str,
                *repeat_id_or_none,
            ) = parts

            sw_score = int(sw_score_str)
            div_percent = float(div_str)
            del_percent = float(del_str)
            ins_percent = float(ins_str)
            seq_start = int(seq_start_str)
            seq_end = int(seq_end_str)
            seq_left = int(seq_left_str[1:-1])
            strand = {"+": "+", "C": "-"}[strand_str]

            if "/" in repeat_class_family:
                repeat_class, repeat_family = repeat_class_family.split("/", 1)
            else:
                repeat_class, repeat_family = repeat_class_family, None
            if strand_str == "C":
                repeat_start_str, repeat_left_str = (repeat_left_str, repeat_start_str)
            repeat_start = int(repeat_start_str)
            repeat_end = int(repeat_end_str)
            repeat_left = int(repeat_left_str[1:-1])
            repeat_id = repeat_id_or_none[0] if repeat_id_or_none else None
            yield Repeat(
                sw_score=sw_score,
                divergence_percent=div_percent,
                insertion_percent=ins_percent,
                deletion_percent=del_percent,
                sequence_id=seq_id,
                sequence_start_original=seq_start,
                sequence_end_original=seq_end,
                sequence_start_c=seq_start - 1,
                sequence_end_c=seq_end,
                sequence_left=seq_left,
                strand=strand,
                repeat_name=repeat_name,
                repeat_class=repeat_class,
                repeat_family=repeat_family,
                repeat_start_original=repeat_start,
                repeat_end_original=repeat_end,
                repeat_start_c=repeat_start - 1,
                repeat_end_c=repeat_end,
                repeat_left=repeat_left,
                repeat_id=repeat_id,
            )


if __name__ == "__main__":
    for path in sys.argv[1:]:
        with RepeatMaskerReader(path) as r:
            repeats_per_class = Counter(repeat.repeat_class for repeat in r)
        print(f"Repeat classes in {path}:")
        for k, v in repeats_per_class.most_common():
            print(f"    {k}: {v} repeats")
