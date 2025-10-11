# biofiles

Pure-Python, zero-dependency collection of bioinformatics-related 
file readers and writers.

## Installation

```shell
python -m pip install biofiles
```

## Usage

Reading FASTA files:

```python
from biofiles.fasta import FASTAReader

with FASTAReader("sequences.fasta") as r:
    for seq in r:
        print(seq.id, len(seq.sequence))

# or

with open("sequences.fasta") as f:
    r = FASTAReader(f)
    for seq in r:
        print(seq.id, len(seq.sequence))
```

Writing FASTA files:

```python
from biofiles.fasta import FASTAWriter
from biofiles.types.sequence import Sequence

seq = Sequence(id="SEQ", description="Important sequence", sequence="GAGAGA")

with FASTAWriter("output.fasta") as w:
    w.write(seq)
```

Reading GFF genome annotations:

```python
from biofiles.gff import GFFReader
from biofiles.dialects.gencode import GENCODE_DIALECT
from biofiles.dialects.genomic_base import Gene

with GFFReader("GCF_009914755.1_T2T-CHM13v2.0_genomic.gff", dialect=GENCODE_DIALECT) as r:
    for feature in r:
        if isinstance(feature, Gene):
            print(feature.name, len(feature.exons))
```

Currently three dialects are supported:
* `biofiles.dialects.gencode.GENCODE_DIALECT` for GENCODE genome annotation;
* `biofiles.dialects.refseq.REFSEQ_DIALECT` for RefSeq genome annotation;
* `biofiles.dialects.stringtie.STRINGTIE_DIALECT` for StringTie output files.

## License 

MIT license, see [License](LICENSE).
