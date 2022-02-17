from typing import (
    TextIO, Iterator
)


def scan_reads(f: TextIO) -> Iterator[tuple[str, str]]:
    """Read sequences from a SimpleFASTQ format."""
    # This is a fucking hack, but Python won't let you check for EOF
    # and when parsing four lines at a time, that makes things harder
    # than they have any good reason for being.
    itr = iter(f)
    try:
        while True:
            name = next(itr).strip()[1:]
            seq = next(itr).strip()
            # removed the next two lines to get "Simple FASTQ"
            # next(itr)
            # qual = next(itr).strip()
            yield name, seq
    except StopIteration:
        return
