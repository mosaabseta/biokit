"""Tests for edit-distance variants and de novo assembly."""

import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from biokit import (
    scs, scs_all, assemble_greedy_contigs, assemble_greedy_indexed,
    assemble_greedy_contigs_indexed, greedy_scs_indexed, de_bruijn_contigs,
    assemble_from_reads,
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


def test_scs_exact():
    assert scs(["ACGGTACGAGC", "GAGCTTCGGA", "GACACGG"]) == "GACACGGTACGAGCTTCGGA"


def test_scs_all_ties():
    length, sols = scs_all(["ABC", "BCA", "CAB"])
    assert length == 5
    assert sols == ["ABCAB", "BCABC", "CABCA"]
    assert scs(["ABC", "BCA", "CAB"]) in sols   # scs returns one of the ties


def test_scs_never_shorter_than_greedy():
    """Exact SCS is optimal, so it is never longer than the greedy result."""
    random.seed(21)
    for _ in range(100):
        seq = "".join(random.choice("ACGT") for _ in range(random.randint(10, 16)))
        reads = [seq[i:i + 5] for i in range(0, len(seq) - 4, 3)][:5]
        assert len(scs(reads)) <= len(greedy_scs(reads, k=1))


def test_assemble_greedy_contigs():
    seq = "ACGTTGCATTGCAAGGCTA"
    reads = [seq[i:i + 8] for i in range(len(seq) - 7)]
    assert assemble_greedy_contigs(reads, 4) == [seq]


def test_assemble_greedy_contigs_gap_gives_two():
    """A coverage gap must yield separate contigs, not one fused sequence."""
    left = "ACGTTGCATTG"
    right = "GGGAAACCCTT"
    reads = ([left[i:i + 6] for i in range(len(left) - 5)] +
             [right[i:i + 6] for i in range(len(right) - 5)])
    contigs = assemble_greedy_contigs(reads, 4)
    assert len(contigs) == 2
    assert set(contigs) == {left, right}


def test_assemble_from_reads_de_bruijn():
    seq = "ACGTTGCATTGCAAGGCTA"
    k = 7
    kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    assert assemble_from_reads(kmers, k=k) == seq


def test_assemble_from_reads_greedy():
    seq = "ACGTTGCATTGCAAGGCTA"
    reads = [seq[i:i + 8] for i in range(len(seq) - 7)]
    assert assemble_from_reads(reads, min_overlap=4, method="greedy") == [seq]


def test_assemble_from_reads_bad_method():
    try:
        assemble_from_reads(["ACGT"], method="nope")
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError")


def test_indexed_matches_classic_greedy():
    """The fast assembler must reproduce the classic greedy result."""
    random.seed(17)
    for _ in range(200):
        glen = random.randint(40, 120)
        g = "".join(random.choice("ACGT") for _ in range(glen))
        rl = random.randint(8, 15)
        n = random.randint(5, 20)
        reads = [g[i:i + rl]
                 for i in sorted(random.sample(range(glen - rl),
                                               min(n, glen - rl)))]
        mo = random.randint(3, rl - 2)
        assert assemble_greedy_indexed(reads, mo) == \
               assemble_greedy_contigs(reads, mo)


def test_indexed_rejects_bad_min_overlap():
    try:
        assemble_greedy_indexed(["ACGT"], 0)
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError")


def test_greedy_scs_indexed_matches():
    seq = "ACGTTGCATTGCAAGGCTA"
    reads = [seq[i:i + 8] for i in range(len(seq) - 7)]
    assert greedy_scs_indexed(reads, 4) == seq


def test_alias_is_same_function():
    assert assemble_greedy_indexed is assemble_greedy_contigs_indexed


def test_de_bruijn_contigs_single():
    seq = "ACGTTGCATTGCAAGGCTA"
    reads = [seq[i:i + 8] for i in range(len(seq) - 7)]
    assert de_bruijn_contigs(reads, 7) == [seq]


def test_de_bruijn_contigs_breaks_at_repeat():
    """Too small a k cannot resolve a repeat, so the graph must break up."""
    seq = "ACGTTGCATTGCAAGGCTA"          # 'TTGCA' occurs twice
    reads = [seq[i:i + 8] for i in range(len(seq) - 7)]
    assert len(de_bruijn_contigs(reads, 5)) > 1   # repeat unresolved
    assert de_bruijn_contigs(reads, 7) == [seq]   # k past the repeat


def test_de_bruijn_contigs_preserves_kmers():
    random.seed(23)
    for _ in range(150):
        glen = random.randint(40, 150)
        g = "".join(random.choice("ACGT") for _ in range(glen))
        rl = random.randint(10, 20)
        k = random.randint(5, rl - 2)
        reads = [g[i:i + rl] for i in range(glen - rl + 1)]
        contigs = de_bruijn_contigs(reads, k)
        assert contigs
        want = {r[i:i + k] for r in reads for i in range(len(r) - k + 1)}
        got = {c[i:i + k] for c in contigs for i in range(len(c) - k + 1)}
        assert want <= got


def test_de_bruijn_contigs_never_none():
    """Unlike assemble_de_bruijn, this always returns contigs."""
    reads = ["ACGTAC", "TTTTGG"]          # disconnected, no Eulerian path
    assert assemble_de_bruijn(reads, 4) is None
    assert len(de_bruijn_contigs(reads, 4)) >= 2


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print("PASS", fn.__name__)
    print("\nAll %d tests passed." % len(fns))
