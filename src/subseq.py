"""Sequence wrappers for constant time slicing."""

from __future__ import annotations

from typing import (
    MutableSequence,
    TypeVar,
    Protocol,
    Generic,
    Sequence,
    Optional,
    Iterator,
    cast,
    overload
)

# Type specifications...
T = TypeVar('T')
S = TypeVar('S', bound='SubSeq')  # type: ignore
C = TypeVar('C', covariant=True)


class Ordered(Protocol[C]):
    """Types that can be ordered with < and >."""

    def __lt__(self, other: object) -> bool:
        """Must support <."""
        ...  # pragma: no cover

    def __gt__(self, other: object) -> bool:
        """Must support >."""
        ...  # pragma: no cover


# Then the functional stuff...
class SubSeq(Generic[T], Sequence[T]):
    """A wrapper around lists and strings for constant time slicing.
    When a sequence x is wrapped in this, x[i:j] takes constant time, since
    we do not copy the sub-sequence but merely wrap a pointer into it.
    There is a slight overhead in dispatching indexing, though, since
    Python isn't statically typed and compiled.
    """

    _x: Sequence[T]
    _i: int
    _j: int

    @staticmethod
    def _fix_index(x: Sequence[T],
                   start: Optional[int],
                   stop: Optional[int]) -> tuple[int, int]:
        """Adjusts indices in the usual slice protocol."""
        start = start if start is not None else 0
        stop = stop if stop is not None else len(x)
        if start < 0:
            start += len(x)
        if stop < 0:
            stop += len(x)

        assert start <= stop, "Start must come before end."
        assert 0 <= start <= len(x), \
            "Indices must be within the sequence's range."
        assert 0 <= stop <= len(x), \
            "Indices must be within the sequence's range."

        return start, stop

    def __init__(self,
                 x: Sequence[T],
                 start: Optional[int] = None,
                 stop: Optional[int] = None
                 ) -> None:
        """Construct new subseq from a slice of x."""
        self._x = x
        self._i, self._j = SubSeq._fix_index(x, start, stop)

    def __repr__(self) -> str:  # pragma: no cover
        """Debugging output."""
        cls_name = self.__class__.__name__
        return f"{cls_name}(x={repr(self._x)}, start={self._i}, stop={self._j})"  # noqa: E501

    def __iter__(self) -> Iterator[T]:
        """Iterate through all elements in the subsequence."""
        return (self._x[i] for i in range(self._i, self._j))

    def __len__(self) -> int:
        """Length of the subsequence."""
        return self._j - self._i

    def __bool__(self) -> bool:
        """Is the subsequence non-empty."""
        return self._i < self._j

    def __str__(self) -> str:
        """Get string representation of subsequence."""
        return str(self._x[self._i:self._j])

    def __eq__(self, other: object) -> bool:
        """Test if this object is equilvalent to other."""
        # duck typing from here on... __eq__ must accept all objects
        # but we can only handle sequences.
        other = cast(Sequence[T], other)
        return len(self) == len(other) and \
            all(a == b for a, b in zip(self, other))

    # You can move this to a mixin if you need to deal with
    # types that do not have an ordering.
    def __lt__(self, other: object) -> bool:
        """Test if self is lexicographically less than other."""
        # duck typing from here on... __eq__ must accept all objects
        # but we can only handle sequences.
        other = cast(Sequence[Ordered[T]], other)
        for a, b in zip(self, other):
            if a < b:
                return True    # noqal
            if a > b:
                return False   # noqal
        return len(self) < len(other)

    @overload
    def __getitem__(self: S, idx: int) -> T:
        """Get a single element."""
        ...  # pragma: no cover

    @overload
    def __getitem__(self: S, idx: slice) -> S:
        """Get a slice of elements as a subseq."""
        ...  # pragma: no cover

    def __getitem__(self: S, idx: int | slice) -> T | S:
        """Get the value at an index, or a new subseq for a slice."""
        if isinstance(idx, int):
            return cast(T, self._x[self._i + idx])

        if isinstance(idx, slice):
            assert idx.step is None, \
                "Subsequences do not handle steps in slices"
            i, j = SubSeq._fix_index(self, idx.start, idx.stop)
            return self.__class__(self._x, self._i + i, self._i + j)

        assert False, "idx of invalid type"  # pragma: no cover
        return None  # just for the stupid linter...

    # FIXME: this shouldn't be here if the sequence is mutable, but the type annotations
    # are a hell to implement a mutable version explicitly... I have done so in pystr,
    # but I won't be bothering here.
    def __setitem__(self, idx: int | slice, val: T) -> None:
        """Set an index or a slice to a value."""
        x = cast(MutableSequence, self._x)
        if isinstance(idx, int):
            x[self._i + idx] = val
        else:
            # I don't handle steps
            assert isinstance(idx, slice) and idx.step is None
            start, stop = SubSeq._fix_index(self, idx.start, idx.stop)
            for i in range(start, stop):
                x[self._i + i] = val
