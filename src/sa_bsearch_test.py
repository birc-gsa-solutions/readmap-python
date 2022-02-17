from sa_bsearch import (
    search_space,
    search_range,
    lower, upper, block,
    sa_bsearch
)

from sais import sais


def test_mississippi():
    x = "mississippi"
    sa = sais(x)

    print("Suffixes")
    for i, j in enumerate(sa):
        print(f"{i:>2} {j:>2}", x[j:])
    print()

    space = search_space(x, sa)
    srange = search_range(0, 1, len(sa))

    # we don't find "a" but find the lower and upper bound at 1
    assert lower("a", srange, space) == 1
    assert upper("a", srange, space) == 1
    res = block("a", srange, space)
    assert res.offset == 1 and res.lo == 1 and res.hi == 1

    # we don't find "x" but find the lower and upper bound at len(sa)
    assert lower("x", srange, space) == len(sa)
    assert upper("x", srange, space) == len(sa)
    res = block("x", srange, space)
    assert res.offset == 1 and res.lo == len(sa) and res.hi == len(sa)

    assert lower("i", srange, space) == 1
    assert upper("i", srange, space) == 5
    res = block("i", srange, space)
    assert res.offset == 1 and res.lo == 1 and res.hi == 5

    assert lower("p", srange, space) == 6
    assert upper("p", srange, space) == 8
    res = block("p", srange, space)
    assert res.offset == 1 and res.lo == 6 and res.hi == 8

    res = sa_bsearch("a", x, sa)
    assert list(res) == []

    res = sa_bsearch("x", x, sa)
    assert list(res) == []

    res = sa_bsearch("i", x, sa)
    assert list(res) == [sa[j] for j in range(1, 5)]

    res = sa_bsearch("p", x, sa)
    assert list(res) == [sa[j] for j in range(6, 8)]

    res = sa_bsearch("si", x, sa)
    assert list(res) == [sa[j] for j in range(8, 10)]


if __name__ == '__main__':
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            print(name)
            f()
