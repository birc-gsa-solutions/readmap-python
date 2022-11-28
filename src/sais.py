"""Implementation of the SAIS algorithm."""

from array import array
from typing import (
    Callable,
    Optional,
    TypeVar
)
import numpy as np
import numpy.typing as npt
from collections import Counter

from alphabet import Alphabet
from bitarray import bitarray

T = TypeVar('T')
UNDEFINED = -1  # Undefined val in SA


def classify_sl(is_s: bitarray, x: memoryview) -> None:
    """Classify positions into S or L."""
    last = len(x) - 1
    is_s[last] = True
    for i in range(len(x)-2, -1, -1):
        is_s[i] = x[i] < x[i + 1] or (x[i] == x[i + 1] and is_s[i + 1])


def compute_buckets(counts: Counter, asize: int) -> list[int]:
    buckets = [0] * (asize + 1)  # np.zeros(asize+1, dtype=np.int32)
    for a in range(1, asize+1):
        buckets[a] = buckets[a-1] + counts[a-1]
    return buckets


def bucket_lms(x: memoryview, asize: int,
               sa: memoryview, counts: Counter,
               is_s: bitarray) \
        -> None:
    """Place LMS strings in their correct buckets."""
    buckets = compute_buckets(counts, asize)
    for i in range(len(sa)):
        sa[i] = UNDEFINED
    for i, _ in enumerate(x):
        if is_s[i] and not is_s[i-1]:
            buckets[x[i]+1] -= 1
            sa[buckets[x[i]+1]] = i


def induce_l(x: memoryview, asize: int,
             sa: memoryview, counts: Counter,
             is_s: bitarray) \
        -> None:
    """Induce L suffixes from the LMS strings."""
    buckets = compute_buckets(counts, asize)
    for i in range(len(x)):
        j = sa[i] - 1
        if sa[i] in (0, UNDEFINED) or is_s[j]:
            continue
        sa[buckets[x[j]]] = j
        buckets[x[j]] += 1


def induce_s(x: memoryview, asize: int,
             sa: memoryview, counts: Counter,
             is_s: bitarray) \
        -> None:
    """Induce S suffixes from the L suffixes."""
    buckets = compute_buckets(counts, asize)
    for i in reversed(range(len(x))):
        j = sa[i] - 1
        if sa[i] == 0 or not is_s[j]:
            continue
        buckets[x[j]+1] -= 1
        sa[buckets[x[j]+1]] = j


def equal_lms(x: memoryview, is_s: bitarray, i: int, j: int) -> bool:
    """Test if two LMS strings are identical."""
    if i == j:
        # This happens as a special case in the beginning of placing them.
        return True

    k = 0
    while True:
        ik, jk = i + k, j + k
        i_lms = is_s[ik] and not is_s[ik - 1]
        j_lms = is_s[jk] and not is_s[jk - 1]
        if k > 0 and i_lms and j_lms:
            return True
        if i_lms != j_lms or x[ik] != x[jk]:
            return False
        k += 1

    # This assert is only hear to help the linter...
    # (checker doesn't understand infinite generators yet)
    assert False, "We only leave the loop with a return."  # pragma: no cover
    return False  # just for the linter


def reduce_lms(x: memoryview, sa: memoryview, is_s: bitarray) \
        -> tuple[memoryview, memoryview, int]:
    """Construct reduced string from LMS strings."""
    # Compact all the LMS indices in the first
    # part of the suffix array...
    k = 0
    for i in sa:
        if is_s[i] and not is_s[i-1]:
            sa[k] = i
            k += 1

    # Create the alphabet and write the translation
    # into the buffer in the right order
    compact, buffer = sa[:k], sa[k:]
    for i in range(len(buffer)):
        buffer[i] = UNDEFINED
    prev, letter = compact[0], 0
    for j in compact:
        if not equal_lms(x, is_s, prev, j):
            letter += 1
        buffer[j // 2] = letter
        prev = j

    # Then compact the buffer into the reduced string
    kk = 0
    for i in buffer:
        if i != UNDEFINED:
            buffer[kk] = i
            kk += 1

    return buffer[:k], compact, letter + 1


def reverse_reduction(x: memoryview, asize: int,
                      sa: memoryview,
                      offsets: memoryview,
                      red_sa: memoryview,
                      counts: Counter,
                      is_s: bitarray) -> None:
    """Get the LMS string order back from the reduced suffix array."""
    # Work out where the LMS strings are in the
    # original string. Compact those indices
    # into the buffer offsets
    k = 0
    for i in range(len(x)):
        if is_s[i] and not is_s[i-1]:
            offsets[k] = i
            k += 1

    # Compact the original indices into sa
    for i, j in enumerate(red_sa):
        sa[i] = offsets[j]

    # Mark the sa after the LMS indices as undefined
    for i in range(len(red_sa), len(sa)):
        sa[i] = UNDEFINED

    buckets = compute_buckets(counts, asize)
    for i in reversed(range(len(red_sa))):
        j, red_sa[i] = red_sa[i], UNDEFINED
        buckets[x[j]+1] -= 1
        sa[buckets[x[j]+1]] = j


@profile
def sais_rec(x: memoryview, sa: memoryview,
             asize: int, is_s: bitarray) -> None:
    """Recursive SAIS algorithm."""
    if len(x) == asize:
        # base case...
        for i, a in enumerate(x):
            sa[a] = i

    else:  # recursive case...
        classify_sl(is_s, x)
        counts = Counter(x)
        bucket_lms(x, asize, sa, counts, is_s)
        induce_l(x, asize, sa, counts, is_s)
        induce_s(x, asize, sa, counts, is_s)

        red, red_sa, red_asize = reduce_lms(x, sa, is_s)

        sais_rec(red, red_sa, red_asize, is_s)
        # restore state...
        classify_sl(is_s, x)

        reverse_reduction(x, asize, sa, red, red_sa, counts, is_s)
        induce_l(x, asize, sa, counts, is_s)
        induce_s(x, asize, sa, counts, is_s)


def sais_alphabet(x: array, alpha: Alphabet) -> array:
    """Run the sais algorithm from a subsequence and an alphabet."""
    sa = array('l', [0] * len(x))
    is_s = bitarray(len(x))
    sais_rec(memoryview(x), memoryview(sa), len(alpha), is_s)
    return sa


def sais(x: str) -> array:
    """Run the sais algorithm from a string."""
    x_, alpha = Alphabet.mapped_subseq_with_sentinel(x)
    return sais_alphabet(x_, alpha)
