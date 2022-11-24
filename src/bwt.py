"""Implementatin of the Burrows-Wheeler transform and related algorithms."""

import numpy as np
import numpy.typing as npt
from typing import (
    Iterator, Callable,
    NamedTuple
)
from array import array
from alphabet import Alphabet
from sais import sais_alphabet
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


def burrows_wheeler_transform(x: str) -> tuple[bytearray, Alphabet, array]:
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


class CTable:
    """
    C-table for other bwt/fm-index search algorithms.

    for CTable ctab, ctab[⍺] is the number of occurrences
    of letters a < ⍺ in the bwt string (or the orignal string,
    since they have the same letters).
    """

    _cumsum: list[int]

    def __init__(self, bwt: bytearray, asize: int) -> None:
        """
        Construct a C-table.

        Compute the C-table from the bwt transformed string and
        the alphabet size.
        """
        # Count occurrences of characters in bwt
        counts = [0] * asize
        for a in bwt:
            counts[a] += 1
        # Get the cumulative sum
        n = 0
        for a, count in enumerate(counts):
            counts[a] = n
            n += count
        # That is all we need...
        self._cumsum = counts

    def __getitem__(self, a: int) -> int:
        """Get the number of occurrences of letters in the bwt less than a."""
        return self._cumsum[a]


class OTable:
    """
    O-table for the FM-index based search.

    For OTable otab, otab[a,i] is the number of occurrences j < i
    where bwt[j] == a.
    """

    # _tbl: list[list[int]]
    _tbl = npt.DTypeLike

    def __init__(self, bwt: bytearray, asize: int) -> None:
        """
        Create O-table.

        Compute the O-table from the bwt transformed string and the size
        of the alphabet the bwt string is over.
        """
        nrow = asize
        ncol = len(bwt) + 1

        self._tbl = np.zeros((nrow, ncol), dtype='i')

        char_inc = np.eye(asize, asize)
        for i in range(1, ncol):
            prev = self._tbl[:, i-1]
            add = char_inc[bwt[i-1], :]
            self._tbl[:, i] = self._tbl[:, i-1] + char_inc[bwt[i-1], :]

    def __getitem__(self, idx: tuple[int, int]) -> int:
        """
        Get the number of occurrences j < i where bwt[j] == a.

        a is the first and i the second value in the idx tuple.
        """
        a, i = idx
        return self._tbl[a, i]


class FMIndexTables(NamedTuple):
    """Preprocessed FMIndex tables."""

    alpha: Alphabet
    sa: list[int]
    ctab: CTable
    otab: OTable
    rotab: OTable


def preprocess_tables(x: str) -> FMIndexTables:
    """Preprocess tables for exact FM/bwt search."""
    bwt, alpha, sa = burrows_wheeler_transform(x)
    ctab = CTable(bwt, len(alpha))
    otab = OTable(bwt, len(alpha))
    bwt, alpha, _ = burrows_wheeler_transform(x[::-1])
    rotab = OTable(bwt, len(alpha))
    return FMIndexTables(alpha, sa, ctab, otab, rotab)


class LiDurbinState(NamedTuple):
    """Full state when running the Li & Durbin algorithm."""

    alpha: Alphabet
    sa: list[int]
    ctab: CTable
    otab: OTable
    rotab: OTable
    dtab: list[int]
    edit_ops: list[Edit]
    p: bytearray


def do_m(
    tbls: LiDurbinState,
    i: int, left: int, right: int,
    edits: int
) -> Iterator[tuple[int, str]]:
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


def do_i(
    tbls: LiDurbinState,
    i: int, left: int, right: int,
    edits: int
) -> Iterator[tuple[int, str]]:
    """Perform an insertion operation in the approx search."""
    tbls.edit_ops.append(Edit.INSERT)
    yield from rec_search(tbls, i - 1, left, right, edits - 1)
    tbls.edit_ops.pop()


def do_d(
    tbls: LiDurbinState,
    i: int, left: int, right: int,
    edits: int
) -> Iterator[tuple[int, str]]:
    """Perform a deletion operation in the approx search."""
    tbls.edit_ops.append(Edit.DELETE)
    for a in range(1, len(tbls.alpha)):
        next_left = tbls.ctab[a] + tbls.otab[a, left]
        next_right = tbls.ctab[a] + tbls.otab[a, right]
        if next_left >= next_right:
            continue
        yield from rec_search(tbls, i, next_left, next_right, edits - 1)
    tbls.edit_ops.pop()


def rec_search(
    tbls: LiDurbinState,
    i: int, left: int, right: int,
    edits: int
) -> Iterator[tuple[int, str]]:
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


def build_dtab(
    p: bytearray, sa: list[int],
    ctab: CTable, rotab: OTable
) -> list[int]:
    """Build the D table for the approximative search."""
    dtab = [0] * len(p)
    min_edits = 0
    left, right, i = 0, len(sa), len(p) - 1
    for i, a in enumerate(p):
        left = ctab[a] + rotab[a, left]
        right = ctab[a] + rotab[a, right]
        if left == right:
            min_edits += 1
            left, right = 0, len(sa)
        dtab[i] = min_edits
    return dtab


def approx_searcher_from_tables(
        alpha: Alphabet,
        sa: list[int],
        ctab: CTable,
        otab: OTable,
        rotab: OTable
) -> ApproxSearchFunc:
    """Build an exact search function from preprocessed tables."""

    def search(p_: str, edits: int) -> Iterator[tuple[int, str]]:
        assert p_, "We can't do approx search with an empty pattern!"
        try:
            p = alpha.map(p_)
        except KeyError:
            return  # can't map, so no matches

        dtab = build_dtab(p, sa, ctab, rotab)
        tbls = LiDurbinState(
            alpha, sa,
            ctab, otab, rotab, dtab,
            list[Edit](), p
        )

        # Do the first operation in this function to avoid
        # deletions in the beginning (end) of the search
        left, right, i = 0, len(sa), len(p) - 1
        yield from do_m(tbls, i, left, right, edits)
        yield from do_i(tbls, i, left, right, edits)

    return search


def approx_preprocess(x: str) -> ApproxSearchFunc:
    """Build an approximative search function for searching in string x."""
    return approx_searcher_from_tables(*preprocess_tables(x))
