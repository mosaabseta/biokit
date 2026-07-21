"""Tests for approximate matching and the index structures.

The critical tests are the fuzz ones: every fast method must return exactly
the same offsets as brute force.
"""

import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from biokit import (
    naive_approximate_match, approximate_match,
    Index, SubseqIndex, query_index, query_subseq, naive_match,
)


# ------------------------------------------------------------- basic cases
def test_basic_approximate():
    assert naive_approximate_match("ACGT", "ACGTACTTAAAA", 2) == [0, 4]
    assert approximate_match("ACGT", "ACGTACTTAAAA", 2) == [0, 4]


def test_zero_mismatches_equals_exact():
    t, p = "ACGTACGTACGT", "ACGT"
    assert approximate_match(p, t, 0) == naive_match(p, t)
    assert naive_approximate_match(p, t, 0) == naive_match(p, t)


def test_count_ops_shapes():
    hits, a, c = naive_approximate_match("ACG", "ACGACG", 1, count_ops=True)
    assert isinstance(a, int) and isinstance(c, int)
    hits, a, c, ih = approximate_match("ACG", "ACGACG", 1, count_ops=True)
    assert ih >= len(hits)   # index hits always >= occurrences


# ------------------------------------------------------------ index basics
def test_kmer_index():
    idx = Index("GGTATTCGGGA", 3)
    assert idx.query("GGTA") == [0]


def test_subseq_index_reference():
    """Reference values from the course material."""
    ind = SubseqIndex("ATATAT", 2, 2)
    assert ind.index == [("AA", 0), ("AA", 2), ("TT", 1), ("TT", 3)]
    assert ind.query("ATATAT") == [0, 2]
    assert ind.query("TATAT") == [1, 3]


def test_subseq_index_span():
    ind = SubseqIndex("ATATAT", 3, 2)
    assert ind.span == 1 + 2 * (3 - 1)   # 5


def test_subseq_query_text():
    t = "to-morrow and to-morrow and to-morrow creeps in this petty pace"
    p = "to-morrow and to-morrow "
    ind = SubseqIndex(t, 8, 3)
    occurrences, hits = query_subseq(p, t, ind, 2)
    assert occurrences == [0, 14]
    assert hits >= len(occurrences)


# ------------------------------------------------------------- guard rails
def test_query_index_rejects_short_pattern():
    t = "ACGTACGTACGT"
    try:
        query_index("ACGT", t, Index(t, 3), 2)   # needs >= 9 chars
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for too-short pattern")


def test_query_subseq_rejects_small_ival():
    t = "ACGTACGTACGTACGTACGT"
    try:
        query_subseq("ACGTACGTACGT", t, SubseqIndex(t, 3, 2), 2)  # ival <= n
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for ival <= max_mismatches")


# ------------------------------------------------------------------- fuzz
def test_all_methods_agree_with_brute_force():
    """Pigeonhole/BM, k-mer index and subseq index vs naive brute force."""
    random.seed(11)
    for _ in range(1500):
        t = "".join(random.choice("ACGT") for _ in range(random.randint(30, 200)))
        plen = random.randint(12, 24)
        st = random.randint(0, len(t) - plen)
        chars = list(t[st:st + plen])
        for _ in range(random.randint(0, 3)):       # inject substitutions
            chars[random.randrange(plen)] = random.choice("ACGT")
        p = "".join(chars)
        n = 2

        truth = naive_approximate_match(p, t, n)
        assert approximate_match(p, t, n) == truth

        k = plen // (n + 1)
        occ, _ = query_index(p, t, Index(t, k), n)
        assert occ == truth

        occ2, _ = query_subseq(p, t, SubseqIndex(t, k, 3), n)
        assert occ2 == truth


def test_pigeonhole_agrees_across_k():
    random.seed(5)
    for _ in range(2000):
        t = "".join(random.choice("ACGT") for _ in range(random.randint(1, 80)))
        p = "".join(random.choice("ACGT") for _ in range(random.randint(1, 12)))
        for k in (0, 1, 2, 3):
            assert naive_approximate_match(p, t, k) == approximate_match(p, t, k)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print("PASS", fn.__name__)
    print("\nAll %d tests passed." % len(fns))
