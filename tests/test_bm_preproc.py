"""Tests for the Boyer-Moore preprocessing tables.

The expected values are the standard published outputs of the Gusfield
Z/N/L'/L/l' algorithms, as used in Ben Langmead's ``bm_preproc.py`` course
material (see ATTRIBUTION in the README). They make an excellent regression
net: an off-by-one anywhere in the table construction shows up immediately.

Run with:  python -m pytest tests/test_bm_preproc.py
       or: python tests/test_bm_preproc.py
"""

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from biokit.matching import (
    z_array, n_array, big_l_prime_array, big_l_array, small_l_prime_array,
    good_suffix_table, good_suffix_mismatch, good_suffix_match,
    dense_bad_char_tab, BoyerMoore, boyer_moore_match, naive_match,
)


# ---------------------------------------------------------------- Z array
def test_z_abb():
    assert z_array("abb") == [3, 0, 0]


def test_z_abababab():
    assert z_array("abababab") == [8, 0, 6, 0, 4, 0, 2, 0]


# ---------------------------------------------------------------- N array
def test_n_abb():
    assert n_array("abb") == [0, 1, 3]


def test_n_abracadabra():
    assert n_array("abracadabra") == [1, 0, 0, 4, 0, 1, 0, 1, 0, 0, 11]


def test_n_abababab():
    assert n_array("abababab") == [0, 2, 0, 4, 0, 6, 0, 8]


# ---------------------------------------------------------------- L' array
def test_big_l_prime_abb():
    assert big_l_prime_array("abb", n_array("abb")) == [0, 0, 2]


def test_big_l_prime_abracadabra():
    s = "abracadabra"
    assert big_l_prime_array(s, n_array(s)) == [0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 8]


# ---------------------------------------------------------------- l' array
def test_small_l_prime_abracadabra():
    expected = [11, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1]
    assert small_l_prime_array(n_array("abracadabra")) == expected


# ------------------------------------------------------- good suffix table
def test_good_suffix_table_abb():
    lp, l, slp = good_suffix_table("abb")
    assert lp == [0, 0, 2]
    assert l == [0, 0, 2]
    assert slp == [3, 0, 0]


def test_good_suffix_table_abracadabra():
    lp, l, slp = good_suffix_table("abracadabra")
    assert lp == [0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 8]
    assert l == [0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 8]
    assert slp == [11, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1]


def test_good_suffix_mismatch_ggtaggt():
    """Full shift table for GGTAGGT, strong (L') vs weak (L) rule."""
    p = "GGTAGGT"
    lp, l, slp = good_suffix_table(p)
    assert lp == [0, 0, 0, 0, 3, 0, 0]
    assert l == [0, 0, 0, 0, 3, 3, 3]
    assert slp == [7, 3, 3, 3, 3, 0, 0]

    # mismatch at the last offset yields no shift from this rule
    assert good_suffix_mismatch(6, lp, slp) == 0
    assert good_suffix_mismatch(6, l, slp) == 0

    # strong rule shifts further than the weak rule at offsets 5 and 4
    assert good_suffix_mismatch(5, lp, slp) == 7
    assert good_suffix_mismatch(5, l, slp) == 4
    assert good_suffix_mismatch(4, lp, slp) == 7
    assert good_suffix_mismatch(4, l, slp) == 4

    # the two rules agree for offsets 3 down to 0
    for i in range(4):
        assert good_suffix_mismatch(i, lp, slp) == 4
        assert good_suffix_mismatch(i, l, slp) == 4


def test_good_suffix_match():
    _, _, slp = good_suffix_table("GGTAGGT")
    assert good_suffix_match(slp) == len(slp) - slp[1]


# ------------------------------------------------------- bad character tab
def test_dense_bad_char_tab():
    amap = {c: i for i, c in enumerate("ACGT")}
    tab = dense_bad_char_tab("ACT", amap)
    assert tab[0] == [0, 0, 0, 0]      # nothing seen yet
    assert tab[1] == [1, 0, 0, 0]      # 'A' seen at index 0
    assert tab[2] == [1, 2, 0, 0]      # 'A' at 0, 'C' at 1


# ------------------------------------------------- class-level integration
def test_boyer_moore_class_tables():
    bm = BoyerMoore("GGTAGGT", alphabet="ACGTN")
    assert bm.big_l_prime == [0, 0, 0, 0, 3, 0, 0]
    assert bm.big_l == [0, 0, 0, 0, 3, 3, 3]
    assert bm.small_l_prime == [7, 3, 3, 3, 3, 0, 0]
    # strong rule shifts at least as far as the weak one
    for i in range(len("GGTAGGT")):
        assert bm.good_suffix_rule(i, strong=True) >= \
               bm.good_suffix_rule(i, strong=False)


def test_bad_character_rule_class():
    bm = BoyerMoore("GGTAGGT", alphabet="ACGTN")
    # mismatching an 'A' at offset 6: last 'A' in P[:6] is at index 3
    assert bm.bad_character_rule(6, "A") == 3


def test_n_alphabet_supported():
    """The default ACGTN alphabet must handle N-containing text."""
    text = "ACGTNNNACGTACGT"
    assert boyer_moore_match("ACGT", text) == naive_match("ACGT", text)


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print("PASS", fn.__name__)
    print("\nAll %d preprocessing tests passed." % len(fns))
