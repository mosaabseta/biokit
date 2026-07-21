"""Tests for biokit. Run with:  python -m pytest  (or python tests/test_biokit.py)"""

import os
import random
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from biokit import (
    reverse_complement, transcribe, translate, gc_content, base_composition,
    find_orfs, longest_orf, longest_orf_in_frame,
    count_kmers, most_frequent_kmer, all_max_frequency_kmers,
    naive_match, boyer_moore_match,
    hamming_distance, edit_distance,
    parse_fasta, write_fasta,
)


def test_sequence_ops():
    assert reverse_complement("ATGC") == "GCAT"
    assert transcribe("ATGC") == "AUGC"
    assert translate("ATGGCCTGA") == "MA*"
    assert translate("ATGGCCTGA", to_stop=True) == "MA"
    assert abs(gc_content("GGCCATAT") - 0.5) < 1e-9
    assert base_composition("AATGC") == {"A": 2, "T": 1, "G": 1, "C": 1}


def test_orf():
    seq = "AAAATGAAATAGCCC"  # ATG AAA TAG starting at index 3
    orfs = find_orfs(seq, both_strands=False)
    assert any(o.sequence == "ATGAAATAG" for o in orfs)
    best = longest_orf("ATGAAATAG", both_strands=False)
    assert best.sequence == "ATGAAATAG"
    assert longest_orf_in_frame("ATGAAATAG", 0) == 9


def test_kmers():
    assert count_kmers("ATATA", 2) == {"AT": 2, "TA": 2}
    km, freq = most_frequent_kmer("ATATA", 2)
    assert freq == 2
    max_freq, winners = all_max_frequency_kmers("ATATA", 2)
    assert max_freq == 2
    assert sorted(winners) == ["AT", "TA"]


def test_distance():
    assert hamming_distance("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT") == 7
    assert edit_distance("kitten", "sitting") == 3
    assert edit_distance("", "abc") == 3
    assert edit_distance("abc", "abc") == 0


def test_naive_basic():
    assert naive_match("AA", "AAAT") == [0, 1]
    assert naive_match("XYZ", "AAAA") == []


def test_boyer_moore_basic():
    assert boyer_moore_match("GGTAGGT", "GGTAGGTGGTAGGT") == [0, 7]


def test_boyer_moore_matches_naive_random():
    """The heart of it: BM must agree with naive on many random inputs."""
    random.seed(42)
    alphabet = "ACGT"
    for _ in range(2000):
        text = "".join(random.choice(alphabet) for _ in range(random.randint(1, 60)))
        m = random.randint(1, 8)
        pattern = "".join(random.choice(alphabet) for _ in range(m))
        assert naive_match(pattern, text) == boyer_moore_match(pattern, text), \
            "mismatch on pattern=%r text=%r" % (pattern, text)


def test_boyer_moore_embedded_patterns():
    """Also check when the pattern is drawn from the text (guaranteed hits)."""
    random.seed(7)
    alphabet = "ACGT"
    for _ in range(1000):
        text = "".join(random.choice(alphabet) for _ in range(random.randint(5, 80)))
        start = random.randint(0, len(text) - 1)
        length = random.randint(1, min(6, len(text) - start))
        pattern = text[start:start + length]
        assert naive_match(pattern, text) == boyer_moore_match(pattern, text)


def test_fasta_roundtrip():
    """Write a FASTA file and read it back unchanged (cross-platform)."""
    tmp_dir = tempfile.mkdtemp()
    path = os.path.join(tmp_dir, "_biokit_test.fasta")
    try:
        records = {"seq1": "ATGC" * 30, "seq2": "GGGCCC"}
        write_fasta(records, path, line_width=70)
        assert parse_fasta(path) == records
        write_fasta(records, path, line_width=0)   # unwrapped must also work
        assert parse_fasta(path) == records
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print("PASS", fn.__name__)
    print("\nAll %d test functions passed." % len(fns))
