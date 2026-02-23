import csv
import sys
from typing import Iterable

from biofiles.common import Reader
from biofiles.types.transcript_mapping import TranscriptMapping, TranscriptClassCode

_TMAP_HEADER = [
    "ref_gene_id",
    "ref_id",
    "class_code",
    "qry_gene_id",
    "qry_id",
    "num_exons",
    "FPKM",
    "TPM",
    "cov",
    "len",
    "major_iso_id",
    "ref_match_len",
]


class TranscriptMappingReader(Reader):
    def __iter__(self) -> Iterable[TranscriptMapping]:
        r = csv.reader(self._input, delimiter="\t")

        it = iter(r)
        header = next(it)
        if header[:12] != _TMAP_HEADER:
            raise ValueError(f"unsupported .tmap header: {header}")

        for row in it:
            (
                ref_gene_id,
                ref_id,
                class_code,
                qry_gene_id,
                qry_id,
                num_exons,
                fpkm,
                tpm,
                cov,
                len,
                major_iso_id,
                ref_match_len,
            ) = row
            yield TranscriptMapping(
                reference_gene_id=ref_gene_id if ref_gene_id != "-" else None,
                reference_transcript_id=ref_id if ref_id != "-" else None,
                class_code=TranscriptClassCode(class_code),
                query_gene_id=qry_gene_id,
                query_transcript_id=qry_id,
                query_number_of_exons=int(num_exons),
                query_fpkm=float(fpkm),
                query_tpm=float(tpm),
                query_coverage=float(cov),
                query_length_bp=int(len),
                query_major_isoform_gene_id=major_iso_id,
                match_length_bp=int(ref_match_len) if ref_match_len != "-" else None,
            )


if __name__ == "__main__":
    for path in sys.argv[1:]:
        with TranscriptMappingReader(path) as tmr:
            for tm in tmr:
                print(tm)
