"""Read-mapper."""


import argparse
import sys

from preprocess import (
    preprocess,
    load_preprocessed
)
from fasta import read_fasta
from fastq import scan_reads
from sam import ssam_record


def main() -> None:
    """FM-index + Li & Durbin based approximative pattern matching."""
    argparser = argparse.ArgumentParser(
        description="Approximative pattern matching",
        usage="\n\treadmap -p genome\n\treadmap genome reads"
    )
    argparser.add_argument(
        "-p", action="store_true",
        help="preprocess the genome."
    )
    argparser.add_argument(
        "-d", type=int, metavar="integer",
        default=1, help="max edit distance."
    )
    argparser.add_argument(
        "genome",
        help="Simple-FASTA file containing the genome.",
        type=argparse.FileType('r')
    )
    argparser.add_argument(
        "reads", nargs="?",
        help="Simple-FASTQ file containing the reads.",
        type=argparse.FileType('r')
    )
    args = argparser.parse_args()

    if args.p:
        preprocess(read_fasta(args.genome), args.genome.name+".readmap")
    else:
        # here we need the optional argument reads
        if args.reads is None:
            argparser.print_help()
            sys.exit(1)

        genome_searchers = load_preprocessed(args.genome.name+".readmap")
        for read_name, read_seq in scan_reads(args.reads):
            for chr_name, search in genome_searchers.items():
                for i, cigar in search(read_seq, args.d):
                    ssam_record(sys.stdout,
                                read_name, chr_name,
                                i, cigar,
                                read_seq)


if __name__ == '__main__':
    main()
