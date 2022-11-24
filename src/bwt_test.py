"""Test bwt."""

from numpy.testing import assert_array_equal
from test_helpers import check_matches
import alphabet
import approx
import bwt


def test_ctable() -> None:
    """Test C-table."""
    x, alpha = alphabet.Alphabet.mapped_string_with_sentinel("aabca")
    ctab = bwt.CTable(x, len(alpha))
    assert ctab[0] == 0, "Nothing is smaller than the sentinel"
    assert ctab[1] == 1, "$ is smaller than 'a'"
    assert ctab[2] == 4, "$ + three 'a'"
    assert ctab[3] == 5, "$ + three 'a' + one 'b'"


def test_otable() -> None:
    """Test O-table."""
    x = "aabca"
    transformed, alpha, _ = bwt.burrows_wheeler_transform(x)
    assert transformed == bytearray([1, 3, 0, 1, 1, 2])

    # we shouldn't look at private members, of course, but
    # we are only testing...
    # pylint: disable=protected-access
    otab = bwt.OTable(transformed, len(alpha))
    assert otab._tbl.shape == (len(alpha), len(transformed)+1)
    assert_array_equal(otab._tbl[0, :], [0, 0, 0, 1, 1, 1, 1], "$ counts")
    assert_array_equal(otab._tbl[1, :], [0, 1, 1, 1, 2, 3, 3], "a counts")
    assert_array_equal(otab._tbl[2, :], [0, 0, 0, 0, 0, 0, 1], "b counts")
    assert_array_equal(otab._tbl[3, :], [0, 0, 1, 1, 1, 1, 1], "c counts")


def test_mississippi_aprox_0() -> None:
    """Test approximative algorithm with edit 0."""
    x = "mississippi"
    search = bwt.approx_preprocess(x)
    for p in ("si", "ppi", "ssi", "pip", "x"):
        matches = list(search(p, 0))
        print(p, matches)
        check_matches(x, p, [idx for idx, _ in matches])


def test_mississippi_aprox_edit() -> None:
    """Test approximative matching."""
    x = "mississippi"
    search = bwt.approx_preprocess(x)
    for edits in [1, 2, 3]:
        for p in ("si", "ppi", "ssi", "pip", "x"):
            for pos, cigar in search(p, edits):
                align = approx.extract_alignment(x, p, pos, cigar)
                assert approx.count_edits(align) <= edits
