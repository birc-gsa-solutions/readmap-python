"""Implementatin of the Burrows-Wheeler transform and related algorithms."""

from typing import (
    Iterator,
    Callable,
    NamedTuple
)
import numpy as np
import numpy.typing as npt
from array import array
from collections import Counter
from alphabet import Alphabet
from sais import sais_alphabet
# from prefix_dub import prefix_doubling
from approx import (
    Edit,
    edits_to_cigar
)

ApproxSearchFunc = Callable[
    # p + dist
    [str, int],
    # pos + CIGAR
    Iterator[tuple[int, str]]
]


def burrows_wheeler_transform(x: str) -> \
        tuple[bytearray, Alphabet, array]:
    """
    Construct the Burrows-Wheeler transform.

    Build the bwt string from a string x, first mapping
    x to a new alphabet.

    Returns the transformed string, the mapping alphabet,
    and the suffix array over x.
    """
    x_, alpha = Alphabet.mapped_string_with_sentinel(x)
    sa = sais_alphabet(x_, alpha)
    bwt = bytearray(x_[j - 1] for j in sa)
    return bwt, alpha, sa


def build_ctab(bwt: bytearray, asize: int) -> npt.NDArray[np.int32]:
    """
    Construct a C-table.

    Compute the C-table from the bwt transformed string and
    the alphabet size.

    For CTable ctab, ctab[⍺] is the number of occurrences
    of letters a < ⍺ in the bwt string (or the orignal string,
    since they have the same letters).
    """
    # Count occurrences of characters in bwt
    counts = Counter(bwt)
    tab = np.zeros(asize, dtype=np.int32)
    for a in range(1, asize):
        tab[a] = counts[a-1] + tab[a-1]
    return tab


def build_otab(bwt: bytearray, asize: int) -> npt.NDArray[np.int32]:
    """
    Create O-table.

    Compute the O-table from the bwt transformed string and the size
    of the alphabet the bwt string is over.
    """
    nrow = asize
    ncol = len(bwt) + 1

    tbl = np.zeros((nrow, ncol), dtype=np.int32)
    char_add = np.eye(nrow, nrow, dtype=np.int32)

    for i in range(1, ncol):
        tbl[:, i] = tbl[:, i-1] + char_add[:, bwt[i-1]]

    return tbl


class FMIndexTables(NamedTuple):
    """Preprocessed FMIndex tables."""

    alpha: Alphabet
    sa: array
    ctab: npt.NDArray[np.int32]
    otab: npt.NDArray[np.int32]
    rotab: npt.NDArray[np.int32]
    dtab: list[int]
    edit_ops: list[Edit]
    p: bytearray


def preprocess_tables(x: str) -> \
    tuple[Alphabet, array,
          npt.NDArray[np.int32],
          npt.NDArray[np.int32], npt.NDArray[np.int32]]:
    """Preprocess tables for exact FM/bwt search."""
    bwt, alpha, sa = burrows_wheeler_transform(x)
    ctab = build_ctab(bwt, len(alpha))
    otab = build_otab(bwt, len(alpha))
    bwt, alpha, _ = burrows_wheeler_transform(x[::-1])
    rotab = build_otab(bwt, len(alpha))
    return alpha, sa, ctab, otab, rotab


def build_dtab(
    p: bytearray,
    sa: array,
    ctab: npt.NDArray[np.int32],
    rotab: npt.NDArray[np.int32]
) -> list[int]:
    """Build the D table for the approximative search."""
    dtab = [0] * (len(p) + 1)  # one extra so we have a zero at -1
    min_edits = 0
    left, right = 0, len(sa)
    for i, a in enumerate(p):
        left = int(ctab[a]) + int(rotab[a, left])
        right = int(ctab[a]) + int(rotab[a, right])
        if left == right:
            min_edits += 1
            left, right = 0, len(sa)
        dtab[i] = min_edits
    return dtab


def do_m(tbls: FMIndexTables,
         i: int, left: int, right: int,
         edits: int) -> Iterator[tuple[int, str]]:
    """Perform a match/mismatch operation in the approx search."""
    tbls.edit_ops.append(Edit.MATCH)
    for a in range(1, len(tbls.alpha)):
        next_left = tbls.ctab[a] + tbls.otab[a, left]
        next_right = tbls.ctab[a] + tbls.otab[a, right]
        if next_left >= next_right:
            continue

        next_edits = edits - (a != tbls.p[i])
        yield from rec_search(tbls, i - 1, next_left, next_right, next_edits)
    tbls.edit_ops.pop()


def do_i(tbls: FMIndexTables,
         i: int, left: int, right: int,
         edits: int) -> Iterator[tuple[int, str]]:
    """Perform an insertion operation in the approx search."""
    tbls.edit_ops.append(Edit.INSERT)
    yield from rec_search(tbls, i - 1, left, right, edits - 1)
    tbls.edit_ops.pop()


def do_d(tbls: FMIndexTables,
         i: int, left: int, right: int,
         edits: int) -> Iterator[tuple[int, str]]:
    """Perform a deletion operation in the approx search."""
    tbls.edit_ops.append(Edit.DELETE)
    for a in range(1, len(tbls.alpha)):
        next_left = tbls.ctab[a] + tbls.otab[a, left]
        next_right = tbls.ctab[a] + tbls.otab[a, right]
        if next_left >= next_right:
            continue
        yield from rec_search(tbls, i, next_left, next_right, edits - 1)
    tbls.edit_ops.pop()


def rec_search(tbls: FMIndexTables,
               i: int, left: int, right: int,
               edits: int) \
        -> Iterator[tuple[int, str]]:
    """Handle recursive operations in approx search."""
    # Do we have a match here?
    if i < 0 <= edits:
        # Remember to reverse the operations, since
        # we did the backwards in the bwt search
        cigar = edits_to_cigar(tbls.edit_ops[::-1])
        for j in range(left, right):
            yield tbls.sa[j], cigar
        return

    # Can we get to a match with the edits we have left?
    if edits < tbls.dtab[i]:
        return

    yield from do_m(tbls, i, left, right, edits)
    yield from do_i(tbls, i, left, right, edits)
    yield from do_d(tbls, i, left, right, edits)


def approx_searcher_from_tables(
        alpha: Alphabet,
        sa: array,
        ctab: npt.NDArray[np.int32],
        otab: npt.NDArray[np.int32],
        rotab: npt.NDArray[np.int32]
) -> ApproxSearchFunc:
    """Build an exact search function from preprocessed tables."""
    # @profile
    def search(p_: str, edits: int) -> Iterator[tuple[int, str]]:
        assert p_, "We can't do approx search with an empty pattern!"
        try:
            p = alpha.map(p_)
        except KeyError:
            return  # can't map, so no matches

        dtab = build_dtab(p, sa, ctab, rotab)
        tbls = FMIndexTables(alpha, sa,
                             ctab, otab, rotab, dtab,
                             list[Edit](), p)

        # Do the first operation in this function to avoid
        # deletions in the beginning (end) of the search
        left, right, i = 0, len(sa), len(p) - 1
        yield from do_m(tbls, i, left, right, edits)
        yield from do_i(tbls, i, left, right, edits)

    return search


def approx_preprocess(x: str) -> ApproxSearchFunc:
    """Build an approximative search function for searching in string x."""
    return approx_searcher_from_tables(*preprocess_tables(x))
