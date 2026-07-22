"""Tests for edit-distance variants and de novo assembly."""

import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from biokit import (
    edit_distance, edit_distance_recursive, edit_distance_matrix,
    approximate_edit_distance,
    overlap, overlap_all_pairs, greedy_scs,
    de_bruijn_graph, assemble_de_bruijn,
)


# ------------------------------------------------------------ edit distance
def test_edit_distance_basic():
    assert edit_distance("kitten", "sitting") == 3
    assert edit_distance("", "abc") == 3
    assert edit_distance("abc", "abc") == 0


def test_recursive_matches_dp():
    random.seed(1)
    for _ in range(300):
        a = "".join(random.choice("ACGT") for _ in range(random.randint(0, 7)))
        b = "".join(random.choice("ACGT") for _ in range(random.randint(0, 7)))
        assert edit_distance_recursive(a, b) == edit_distance(a, b)


def test_matrix_corner_is_distance():
    a, b = "GCGTATGC", "TATTGGCTATC"
    assert edit_distance_matrix(a, b)[-1][-1] == edit_distance(a, b)


def test_approximate_edit_distance():
    # exact occurrence -> 0 edits
    assert approximate_edit_distance("GCGTATGC", "TATTGGCTATCGGCGTATGCAAAA") == 0
    # one substitution in the occurrence -> 1 edit
    assert approximate_edit_distance("GCGTATGC", "TATTGGCTATCGGCGTATGGAAAA") == 1
    # pattern absent -> small but nonzero
    assert approximate_edit_distance("TTTT", "ACGACGACG") > 0


def test_approximate_edit_never_exceeds_global():
    """Approx (free start in text) is never worse than the full global distance."""
    random.seed(2)
    for _ in range(200):
        p = "".join(random.choice("AC") for _ in range(random.randint(1, 6)))
        t = "".join(random.choice("AC") for _ in range(random.randint(1, 15)))
        assert approximate_edit_distance(p, t) <= edit_distance(p, t)


# ---------------------------------------------------------------- overlap
def test_overlap():
    assert overlap("TTACGT", "CGTGTGC") == 3
    assert overlap("TTACGT", "GTGTGC") == 0
    assert overlap("TTACGT", "TTACGT") == 6  # whole string


def test_overlap_all_pairs():
    reads = ["ACGGATGATC", "GATCAAGT", "TTCACGGA"]
    edges = set(overlap_all_pairs(reads, 3))
    assert ("ACGGATGATC", "GATCAAGT", 4) in edges
    assert ("TTCACGGA", "ACGGATGATC", 5) in edges
    assert len(edges) == 2


def test_overlap_all_pairs_matches_bruteforce():
    random.seed(3)
    for _ in range(200):
        reads = ["".join(random.choice("ACGT") for _ in range(8))
                 for _ in range(6)]
        k = 3
        fast = set(overlap_all_pairs(reads, k))
        brute = set()
        for a in reads:
            for b in reads:
                if a != b:
                    o = overlap(a, b, k)
                    if o > 0:
                        brute.add((a, b, o))
        assert fast == brute


# --------------------------------------------------------------- assembly
def test_greedy_scs_contains_reads():
    seq = "ACGTACGTACGTTTAA"
    reads = [seq[i:i + 6] for i in range(0, len(seq) - 5, 2)]
    sup = greedy_scs(reads, k=3)
    assert all(r in sup for r in reads)


def test_de_bruijn_reconstructs():
    """A perfect k-mer spectrum reconstructs a sequence with that spectrum."""
    random.seed(7)
    for _ in range(500):
        seq = "".join(random.choice("ACGT") for _ in range(random.randint(30, 80)))
        k = random.randint(6, 12)
        kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
        res = assemble_de_bruijn(kmers, k)
        assert res is not None
        # may differ from seq at repeats, but must share the exact k-mer set
        got = sorted(res[i:i + k] for i in range(len(res) - k + 1))
        assert got == sorted(kmers)


def test_de_bruijn_graph_shape():
    nodes, edges = de_bruijn_graph(["ACGT"], 3)
    assert nodes == {"AC", "CG", "GT"}
    assert edges == [("AC", "CG"), ("CG", "GT")]


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print("PASS", fn.__name__)
    print("\nAll %d tests passed." % len(fns))
