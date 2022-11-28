"""
Suffix array constructing using prefix-doubling.

If you are here, you must want to learn about prefix doubling.
Good for you! You make me proud.

This is a suffix array construction algorithm that, like many others,
are based on bucket sorting, but it doesn't do divide-and-conquer.
Instead, it sorts the suffixes according to progressively longer
prefixes, doubling their length in each iteration (thus the name).

Take the string 'mississippi$' (it is as good an example as any).
The suffixes are shown below, with their length-1 prefixes in square
brackets:

[m]ississippi$
[i]ssissippi$
[s]sissippi$
[s]issippi$
[i]ssippi$
[s]sippi$
[s]ippi$
[i]ppi$
[p]pi$
[p]i$
[i]$
[$]

If we sort them with respect to these length-1 prefixes we get

[$]
[i]ssissippi$
[i]ssippi$
[i]ppi$
[i]$
[m]ississippi$
[p]pi$
[p]i$
[s]sissippi$
[s]issippi$
[s]sippi$
[s]ippi$

The suffixes are not completely sorted, but they are sorted with
respect to their first letter. That's all we get with one sweep of
bucket sort.

Then we sort them with respect to their length-2 prefixes:

[$ ]
[i$]
[ip]pi$
[is]sissippi$
[is]sippi$
[mi]ssissippi$
[pi]$
[pp]i$
[si]ssippi$
[si]ppi$
[ss]issippi$
[ss]ippi$

Again, they are not quite sorted. We can recognise that there might
be more work to do by the prefixes not being unique, similar to
how it works with skew or SAIS; so if the prefixes we have sorted
with respect to aren't unique we need to do more, and we will double
the length of prefixes we look at, now to length 4, and sort again.

[$   ]
[i$  ]
[ippi]$
[issi]ssippi$
[issi]ppi$
[miss]issippi$
[pi$ ]
[ppi$]
[sipp]i$
[siss]ippi$
[ssip]pi$
[ssis]sippi$

Once again the prefixes aren't unique (we have 'issi' twice) and thus
we need another doubling of prefix length, now to length 8.

[$       ]
[i$      ]
[ippi$   ]
[issippi$]
[ississip]pi$
[mississi]ppi$
[pi$     ]
[ppi$    ]
[sippi$  ]
[sissippi]$
[ssippi$ ]
[ssissipp]i$

Now, the prefixes are unique, so we know that we have sorted all
the strings.

What is the worst-case number of times that we need to sort? Well,
if we double the length of the prefixes every time, we will cover
the entire string in O(log n) iterations, and then the prefix must
be unique, so we need at most O(log n) iterations.

Usually, though, prefixes are unique much before. If we had random
strings, the chance that two strings are identical drops off
exponetially, and in that case we only expect a constant number of
iterations. So, in practice, the approach can be quite efficient.

Now, what about the sorting? We can obviously sort the length-1 prefix
in O(n) using a bucket sort, but how do we handle the longer prefixes?

Here, we use a trick similar to what we have seen in skew and sais.
We replace strings with numbers, such that the numbers preserve the
lexicographical ordering of the strings. For length-1 prefixes, we can
just use the letters we have, or we can map them if you want to. Let's
do the mapping:

$ => 0
i => 1
m => 2
p => 3
s => 4

then the length-1 prefix are

 0: [2]ississippi$
 1: [1]ssissippi$
 2: [4]sissippi$
 3: [4]issippi$
 4: [1]ssippi$
 5: [4]sippi$
 6: [4]ippi$
 7: [1]ppi$
 8: [3]pi$
 9: [3]i$
10: [1]$
11: [0]

There is nothing new here. However, when we want to sort with length-2
prefixes we can use that the key for suffix i should be the number for i
followed by the number for i+1:

 0: [2]ississippi$ is followed by 1: [1]ssissippi$ so the pair is [2,1]
 1: [1]ssissippi$ is followed by 2: [4]sissippi$ so the pair is [1,4]

and so on:

 0: [2,1]ssissippi$
 1: [1,4]sissippi$
 2: [4,4]issippi$
 3: [4,1]ssippi$
 4: [1,4]sippi$
 5: [4,4]ippi$
 6: [4,1]ppi$
 7: [1,3]pi$
 8: [3,3]i$
 9: [3,1]$
10: [1,0]
11: [0,0]  // put 0 if there isn't a i+1

Again, there isn't much new here, I've just renamed the letters. However,
I have two numbers of magnitude <= n and I can radix sort these in time
O(n).

11: [0,0]
10: [1,0]
 7: [1,3]pi$
 1: [1,4]sissippi$
 4: [1,4]sippi$
 0: [2,1]ssissippi$
 9: [3,1]$
 8: [2,2]i$
 3: [4,1]ssippi$
 6: [4,1]ppi$
 2: [4,4]issippi$
 5: [4,4]ippi$

and once I have the pairs sorted, I can assign new numbers to the pairs
to get a single number for each of the prefixes.

[0,0] => 0
[1,0] => 1
[1,3] => 2
[1,4] => 3
[1,4] => 3
[2,1] => 4
[3,1] => 5
[2,2] => 6
[4,1] => 7
[4,1] => 7
[4,4] => 8
[4,4] => 8

and then you have the 2-prefix encoding

11: [0]
10: [1]
 7: [2]pi$
 1: [3]sissippi$
 4: [3]sippi$
 0: [4]ssissippi$
 9: [5]$
 8: [6]i$
 3: [7]ssippi$
 6: [7]ppi$
 2: [8]issippi$
 5: [8]ippi$

To get pairs of numbers for the length-4 prefixes, you use the same
approach of pairing up numbers, but now you don't want to look at i
and i+1, but at i and i+2.

When you looked at i and i+1 you concatenated length-1 prefixes

     i
....[x] y  z...
.... x [y] z...
       i+1

but if you look at i and i+2 you can concatenate length-2 prefixes

     i
....[xy]  zw ...
.... xy  [zw] ...
         i+2

In general, when you want to concatenate length-k prefixes, you will
look at i and i+k and take the numbers from these.

For the length-2 prefixes, with the numbers they have are

11: [0]
10: [1]
 7: [2]pi$
 1: [3]sissippi$
 4: [3]sippi$
 0: [4]ssissippi$
 9: [5]$
 8: [6]i$
 3: [7]ssippi$
 6: [7]ppi$
 2: [8]issippi$
 5: [8]ippi$

but let us rearrange them so the suffixes come in the order they have
in x, so it is easier to see which suffixes are two letters apart

 0: [4]ssissippi$
 1: [3]sissippi$
 2: [8]issippi$
 3: [7]ssippi$
 4: [3]sippi$
 5: [8]ippi$
 6: [7]ppi$
 7: [2]pi$
 8: [6]i$
 9: [5]$
10: [1]
11: [0]

To get the pair for the first suffix, look at 0 and 2:

 0: [4]ssissippi$
 2: [8]issippi$

which gives you (4,8). For the second suffix, look at 1 and 3:

 1: [3]sissippi$
 3: [7]ssippi$

so it should get the key (3,7). Continue like that to get pairs for
all the suffixes.

 0: [4,8]issippi$
 1: [3,7]ssippi$
 2: [8,3]sippi$
 3: [7,8]ippi$
 4: [3,7]ppi$
 5: [8,2]pi$
 6: [7,6]i$
 7: [2,5]$
 8: [6,1]
 9: [5,0]
10: [1,0]
11: [0,0]

Sort these

11: [0,0]
10: [1,0]
 7: [2,5]$
 1: [3,7]ssippi$
 4: [3,7]ppi$
 0: [4,8]issippi$
 9: [5,0]
 8: [6,1]
 6: [7,6]i$
 3: [7,8]ippi$
 5: [8,2]pi$
 2: [8,3]sippi$

and assign the suffixes new numbers

11: [0]                [0,0] => 0
10: [1]                [1,0] => 1
 7: [2]$               [2,5] => 2
 1: [3]ssippi$         [3,7] => 3
 4: [3]ppi$            [3,7] => 3
 0: [4]issippi$        [4,8] => 4
 9: [5]                [5,0] => 5
 8: [6]                [6,1] => 6
 6: [7]i$              [7,6] => 7
 3: [8]ippi$           [7,8] => 8
 5: [9]pi$             [8,2] => 9
 2: [10]sippi$         [8,3] => 10

and then do the whole thing again (but now with i and i+4 when
building pairs).

Since you can radix sort pairs (i,j) where i,j <= n in O(n), each
sort iteration runs in O(n), and you can trivially assign numbers
to the suffixes in O(n) as well.

The only remaining caveat is that while, true, you can do a radix sort
of pairs of numbers i,j <= n in O(n), you probably don't have the memory
for it (nor the desire to spend the full time on it either, frankly).

If you need bucket tables of size n, and n is in the hundreds of millions,
you might very well run out of memory. But that isn't a problem. When
I said that you could do radix sort, I didn't mean just two iterations of
bucket sort on the two integers. You can use radix sort on the integers
as well. If, say, the integers are 32 bit, we can split them into 4 bytes
each, so sorting the pair becomes 8 iterations of bucket sort. And it ends
up even better than that, because it would be 8 iterations if we had to sort
the pairs from scratch, but we don't have to. After one iteration, the numbers
are already sorted with respect to the first of the two integers in the pair.
When we construct the pairs (i,j), we can exploit that we already have the
i-component in order. Within each i-block, we need to sort the j's, but that
is just one integer, so four bucket-sorts of bytes.

There are different ways of doing this, but they are pretty much all
highly efficient. I haven't done it in the Python code, though, since
Python doesn't have fixed-sized integers, and it is slightly more complicated
(and a bit more ugly) to radix sort Python's kind of integers. If you want
to see the full solution, check out e.g. the C or Go implementations.
"""

from collections import Counter
from dataclasses import dataclass
from alphabet import Alphabet


@dataclass
class Rank:
    """Wrapper around a list that handles out-of-bounds indexing."""

    ranks: list[int]

    def __setitem__(self, i: int, val: int) -> None:
        """Set the rank for a suffix."""
        self.ranks[i] = val

    def __getitem__(self, i: int) -> int:
        """Get the rank of a suffix."""
        return self.ranks[i] if i < len(self.ranks) else 0


def buckets(keys: list[int]) -> list[int]:
    """Compute bucket indices for keys."""
    counts = Counter(keys)
    buckets = [0] * (max(counts.keys()) + 1)
    for i, _ in enumerate(buckets):
        # for i=0 we rely on wrap around and default values
        buckets[i] = buckets[i-1] + counts[i-1]
    return buckets


def sort_with_rank(sa: list[int], k: int, rank: Rank) -> list[int]:
    """
    Return the suffixes in sa, sorted with respect to rank at offset k.

    Returns a new list of indices that are a stable sort of sa according to
    rank[sa[i]+k].
    """
    keys = [rank[i + k] for i in sa]
    out = [0] * len(sa)
    buck = buckets(keys)
    for j, i in enumerate(sa):
        out[buck[keys[j]]] = i
        buck[keys[j]] += 1
    return out


def sort_pairs(sa: list[int], k: int, rank: Rank) -> list[int]:
    """
    Sort sa as pairs taken from rank[sa[i]] and rank[sa[i]+k].

    The sort is a stable radix where we first sort with respect
    to rank[sa[i]+k] and then follow with a sort of rank[sa[i]].

    We return the new sorted array.
    """
    sa = sort_with_rank(sa, k, rank)
    sa = sort_with_rank(sa, 0, rank)
    return sa


def update_rank(sa: list[int], k: int, rank: Rank) -> tuple[int, Rank]:
    """
    Update the rank according to the new ordering.

    First we extract the pairs (rank[i],rank[i+k]) where we take the indices i
    in the order they appear in sa. That way, the pairs are sorted (as long as
    we have sorted up to prefix length k, which we will have done here).

    After that, it is a simple matter of running through the pairs and building
    an alphabet. Notice the indexing (i indices according to the order of pairs
    while j indices according to sa). This is necessary since the order in sa
    is the curren sorted order while rank always has the suffixes in the order
    at which they appear in the string.
    """
    pairs = [(rank[i], rank[i+k]) for i in sa]
    pairs.append(pairs[0])  # removes a special case in the loop
    a, new_rank = 0, [0] * len(pairs)
    for i, j in enumerate(sa):
        a += pairs[i - 1] != pairs[i]
        new_rank[j] = a
    return a + 1, Rank(new_rank)


def prefix_doubling(x: bytearray, asize: int) -> list[int]:
    """
    Compute the suffix array for x using a least-significant digit radix sort.
    """
    rank = Rank(x)
    sa = sort_with_rank(list(range(len(x))), 0, rank)

    k = 1
    while asize < len(sa):
        sa = sort_pairs(sa, k, rank)
        asize, rank = update_rank(sa, k, rank)
        k *= 2

    return sa
