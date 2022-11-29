"""Implementatin of the Burrows-Wheeler transform and related algorithms."""

import enum
from typing import (
    Iterator, Callable,
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


def preprocess_tables(x: str) -> FMIndexTables:
    """Preprocess tables for exact FM/bwt search."""
    bwt, alpha, sa = burrows_wheeler_transform(x)
    ctab = build_ctab(bwt, len(alpha))
    otab = build_otab(bwt, len(alpha))
    bwt, alpha, _ = burrows_wheeler_transform(x[::-1])
    rotab = build_otab(bwt, len(alpha))
    return FMIndexTables(alpha, sa, ctab, otab, rotab)


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


# State machine
class MatchState(enum.Enum):
    """Edit operations."""

    DONE = enum.auto()
    REC = enum.auto()
    MATCH = enum.auto()
    INSERT = enum.auto()
    DELETE = enum.auto()
    POP_NEXT = enum.auto()


class MatchFrame(NamedTuple):
    """Stack frame when storing processing states for later."""

    state: MatchState  # The kind of state we have
    left: int          # left in search interval
    right: int         # right in search interval
    pos: int           # position into pattern
    a: int             # character
    edits_idx: int     # index into the edits stack
    max_edits: int     # edits we have left


def approx_searcher_from_tables(
        alpha: Alphabet,
        sa: array,
        ctab: npt.NDArray[np.int32],
        otab: npt.NDArray[np.int32],
        rotab: npt.NDArray[np.int32]
) -> ApproxSearchFunc:
    """Build an exact search function from preprocessed tables."""
    # @profile
    def search(p_: str, max_edits: int) -> Iterator[tuple[int, str]]:
        assert p_, "We can't do approx search with an empty pattern!"
        try:
            p = alpha.map(p_)
        except KeyError:
            return  # can't map, so no matches

        dtab = build_dtab(p, sa, ctab, rotab)
        left, right, pos = 0, len(sa), len(p) - 1
        a, edits_idx = -1, 0
        edits: list[Edit] = [Edit.MATCH] * (len(p) + max_edits)
        stack: list[MatchFrame] = [
            MatchFrame(MatchState.DONE, left, right, pos, -1, -1, -1)
        ]

        state = MatchState.REC
        while True:
            match state:
                case MatchState.DONE:
                    # If we get to the bottom of the stack and pop a 'DONE'
                    # state we are done searching.
                    return

                case MatchState.REC:
                    if left == right or max_edits < dtab[pos]:
                        state = MatchState.POP_NEXT
                        continue

                    if pos < 0:
                        # We have a hit!
                        cigar = edits_to_cigar(edits[:edits_idx][::-1])
                        for j in range(left, right):
                            yield sa[j], cigar
                        state = MatchState.POP_NEXT  # go to the next state
                        continue

                    # If we haven't reached the end or have something to report
                    # we continue with a match...
                    a = 1  # start after the sentinel
                    state = MatchState.MATCH
                    continue

                case MatchState.MATCH:
                    if a == len(alpha):
                        # done with matching, so move on to insertion
                        state = MatchState.INSERT
                        continue

                    # remember to do the next match later
                    stack.append(MatchFrame(
                        MatchState.MATCH,
                        left, right, pos, a+1,
                        edits_idx, max_edits
                    ))
                    # but do a recursion here...
                    state = MatchState.REC
                    left = int(ctab[a]) + int(otab[a, left])
                    right = int(ctab[a]) + int(otab[a, right])
                    max_edits -= int(a != p[pos])
                    pos -= 1
                    edits[edits_idx] = Edit.MATCH
                    edits_idx += 1
                    continue

                case MatchState.INSERT:
                    # remember to do a deletion later
                    stack.append(MatchFrame(
                        MatchState.DELETE,
                        left, right, pos, 1,
                        edits_idx, max_edits
                    ))
                    # but for now, continue with a rec
                    state = MatchState.REC
                    pos -= 1
                    max_edits -= 1
                    edits[edits_idx] = Edit.INSERT
                    edits_idx += 1
                    continue

                case MatchState.DELETE:
                    if a == len(alpha) or edits_idx == 0:
                        # we are either done with deletions or shouldn't
                        # start because we are at the beginning of the
                        # recursion
                        state = MatchState.POP_NEXT
                        continue

                    # Push another deletion for later
                    stack.append(MatchFrame(
                        MatchState.DELETE,
                        left, right, pos, a + 1,
                        edits_idx, max_edits
                    ))
                    # and do a recursion for now
                    state = MatchState.REC
                    left = int(ctab[a]) + int(otab[a, left])
                    right = int(ctab[a]) + int(otab[a, right])
                    max_edits -= 1
                    edits[edits_idx] = Edit.DELETE
                    edits_idx += 1
                    continue

                case MatchState.POP_NEXT:
                    state, left, right, pos, a, edits_idx, max_edits \
                        = stack.pop()
                    continue

    return search


def approx_preprocess(x: str) -> ApproxSearchFunc:
    """Build an approximative search function for searching in string x."""
    return approx_searcher_from_tables(*preprocess_tables(x))
