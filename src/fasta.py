from typing import (
    TextIO
)


def read_fasta(f: TextIO) -> dict[str, str]:
    # This is a bit flaky of a fasta parser, but it is okay
    # for this class...
    genome: dict[str, str] = {}
    chromosomes = f.read().split('>')
    for i in range(1, len(chromosomes)):
        name, *seq = chromosomes[i].split()
        genome[name] = ''.join(seq)
    return genome
