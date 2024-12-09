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
from biofiles.types.feature import Gene

with GFFReader("GCF_009914755.1_T2T-CHM13v2.0_genomic.gff") as r:
    for feature in r:
        if isinstance(feature, Gene):
            print(feature.name, len(feature.exons))
```

## License 

MIT license, see [License](LICENSE).
