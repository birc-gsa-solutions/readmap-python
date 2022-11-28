"""Test of sais algorithm."""

from array import array
from test_helpers import (
    check_sorted,
    fibonacci_string,
    random_string
)
from alphabet import Alphabet
from bitarray import bitarray
from sais import classify_sl, sais


def test_remap() -> None:
    """Test that remapping works."""
    x = "mississippi"
    mapped, alpha = Alphabet.mapped_string_with_sentinel(x)
    assert mapped == bytearray([2, 1, 4, 4, 1, 4, 4, 1, 3, 3, 1, 0])
    assert len(alpha) == 5

    for _ in range(10):
        x = random_string(1000)
        mapped, alpha = Alphabet.mapped_string_with_sentinel(x)
        assert set(mapped) == set(range(len(alpha)))


def test_classify() -> None:
    """Test that the SL classification works."""
    # mississippi$
    # LSLLSLLSLLLS
    x, _ = Alphabet.mapped_subseq_with_sentinel("mississippi")
    assert len(x) == len("mississippi") + 1

    is_s = bitarray(len(x))
    assert len(is_s) == len(x)

    classify_sl(is_s, memoryview(x))

    expected = [
        # L    S     L      L      S     L      L
        False, True, False, False, True, False, False,
        # S   L      L      L      S
        True, False, False, False, True
    ]
    for i, b in enumerate(is_s):
        assert b == expected[i]


def test_base_case() -> None:
    """Test that we can sort base cases."""
    assert sais("abc") == array('l', [3, 0, 1, 2])
    assert sais("cba") == array('l', [3, 2, 1, 0])
    assert sais("acb") == array('l', [3, 0, 2, 1])


def test_mississippi() -> None:
    """Test on mississippi."""
    x = "mississippi"
    sa = sais(x)
    assert len(x) == len(sa) - 1
    check_sorted(x, sa)


def test_sais_sorted() -> None:
    """Test that the sais suffix array is sorted."""
    for _ in range(10):
        x = random_string(1000)
        sa = sais(x)
        check_sorted(x, sa)

    for n in range(10, 15):
        x = fibonacci_string(n)
        sa = sais(x)
        check_sorted(x, sa)


if __name__ == '__main__':
    globs = list(globals().items())
    for name, f in globs:
        if name.startswith("test_"):
            print(name)
            f()
