from typing import (
    Iterator, NamedTuple, Optional
)
from dataclasses import (
    dataclass
)


@dataclass
class SARange:
    sa: list[int]
    start: int
    stop: int

    def __iter__(self) -> Iterator[int]:
        for i in range(self.start, self.stop):
            yield self.sa[i]


@dataclass
class InfiniteString:
    """Turns out-of-bounds index into sentinels."""
    _x: str

    def __getitem__(self, i: int) -> str:
        return self._x[i] if i < len(self._x) else chr(0)


class SearchSpace(NamedTuple):
    x: InfiniteString
    sa: list[int]


class SearchRange(NamedTuple):
    offset: int
    lo: int
    hi: int


def search_space(x: str, sa: list[int]) -> SearchSpace:
    return SearchSpace(InfiniteString(x), sa)


def search_range(offset: int, lo: int, hi: int) -> SearchRange:
    return SearchRange(offset, lo, hi)


def lower(a: str, srange: SearchRange, space: SearchSpace) -> int:
    """Finds the lower bound of `a` in the block defined by `srange`."""
    offset, lo, hi = srange
    x, sa = space
    while lo < hi:
        m = (lo + hi) // 2
        if x[sa[m] + offset] < a:
            lo = m + 1
        else:
            hi = m
    return lo


def upper(a: str, srange: SearchRange, space: SearchSpace) -> int:
    """Finds the upper bound of `a` in the block defined by `srange`."""
    return lower(chr(ord(a) + 1), srange, space)


def block(a: str, srange: SearchRange, space: SearchSpace) -> SearchRange:
    """Updates srange by finding the sub-block with `a` at `offset`. The
    result is a search range with an updated `offset` and `[lo,hi)` interval.
    Returns None if the block is empty."""
    return SearchRange(
        srange.offset + 1, lower(a, srange, space), upper(a, srange, space)
    )


def sa_bsearch(p: str, x: str, sa: list[int]) -> SARange:
    space = search_space(x, sa)
    srange = search_range(0, 0, len(sa))
    for a in p:
        srange = block(a, srange, space)
        if srange.lo == srange.hi:
            break
    return SARange(sa, srange.lo, srange.hi)
