"""Exact pattern matching algorithms.

Two matchers are provided, both returning the 0-based offsets where
``pattern`` occurs in ``text``:

* :func:`naive_match` -- the brute-force baseline.
* :class:`BoyerMoore` + :func:`boyer_moore_match` -- Boyer-Moore using
  both the *bad character* and *strong good suffix* rules. The good-suffix
  tables are built from the Z / N arrays (Gusfield, *Algorithms on Strings,
  Trees and Sequences*).

Both matchers can optionally report how many alignments and character
comparisons they performed, so you can compare their work on the same input.
"""


# ---------------------------------------------------------------------------
# Naive exact matching
# ---------------------------------------------------------------------------
def naive_match(pattern, text, count_ops=False):
    """Find all occurrences of ``pattern`` in ``text`` by brute force.

    Slides ``pattern`` across every alignment and compares left to right,
    bailing out on the first mismatch.

    Parameters
    ----------
    pattern, text : str
    count_ops : bool
        If True, also return ``(alignments_tried, character_comparisons)``.

    Returns
    -------
    list[int]                       (count_ops=False)
    tuple[list[int], int, int]      (count_ops=True)

    >>> naive_match("AA", "AAAT")
    [0, 1]
    """
    occurrences = []
    n, m = len(text), len(pattern)
    alignments = 0
    comparisons = 0
    if m == 0 or m > n:
        return (occurrences, 0, 0) if count_ops else occurrences
    for i in range(n - m + 1):
        alignments += 1
        match = True
        for j in range(m):
            comparisons += 1
            if text[i + j] != pattern[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    if count_ops:
        return occurrences, alignments, comparisons
    return occurrences


# ---------------------------------------------------------------------------
# Z / N array preprocessing (Gusfield)
# ---------------------------------------------------------------------------
def _z_array(s):
    """Return the Z-array of ``s`` (Gusfield 1.4.1)."""
    assert len(s) > 1, "Z-array is undefined for strings shorter than 2"
    z = [len(s)] + [0] * (len(s) - 1)

    # initial comparison of s[1:] against the prefix
    for i in range(1, len(s)):
        if s[i] == s[i - 1]:
            z[1] += 1
        else:
            break

    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1

    for k in range(2, len(s)):
        if k > r:  # Case 1: outside any Z-box
            n = 0
            while k + n < len(s) and s[n] == s[k + n]:
                n += 1
            z[k] = n
            if n > 0:
                r, l = k + n - 1, k
        else:  # Case 2: inside a Z-box
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:  # Case 2a
                z[k] = zkp
            else:  # Case 2b: extend past r
                nmatch = 0
                while r + 1 + nmatch < len(s) and \
                        s[r + 1 + nmatch] == s[nbeta + nmatch]:
                    nmatch += 1
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z


def _n_array(s):
    """Return the N-array: Z-array of the reversed string, reversed."""
    return _z_array(s[::-1])[::-1]


def _big_l_prime_array(p, n):
    """Compile the L' array (Gusfield theorem 2.2.2) from ``p`` and ``N``.

    ``L'[i]`` is the largest index j < n such that ``N[j] == |P[i:]|``.
    """
    lp = [0] * len(p)
    for j in range(len(p) - 1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp


def _big_l_array(p, lp):
    """Compile the L array (Gusfield theorem 2.2.2) from ``p`` and ``L'``.

    ``L[i]`` is the largest index j < n such that ``N[j] >= |P[i:]|``.
    """
    big_l = [0] * len(p)
    big_l[1] = lp[1]
    for i in range(2, len(p)):
        big_l[i] = max(big_l[i - 1], lp[i])
    return big_l


def _small_l_prime_array(n):
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i + 1:  # prefix of length i+1 matches a suffix
            small_lp[len(n) - i - 1] = i + 1
    for i in range(len(n) - 2, -1, -1):  # smear values leftward
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i + 1]
    return small_lp


# ---------------------------------------------------------------------------
# Public preprocessing API
#
# These expose the individual Gusfield tables so you can inspect them,
# teach from them, or build your own matcher. The naming follows the
# convention used in Ben Langmead's course material (see ATTRIBUTION in the
# README) so existing code and exercises drop in unchanged.
# ---------------------------------------------------------------------------
def z_array(s):
    """Z-array of ``s`` (Gusfield theorem 1.4.1).

    ``Z[i]`` is the length of the longest substring starting at ``i`` that
    is also a prefix of ``s``.

    >>> z_array("abababab")
    [8, 0, 6, 0, 4, 0, 2, 0]
    """
    return _z_array(s)


def n_array(s):
    """N-array (Gusfield theorem 2.2.2): the Z-array of the reversed string.

    ``N[j]`` is the length of the longest suffix of ``s[:j+1]`` that is also
    a suffix of ``s``.

    >>> n_array("abb")
    [0, 1, 3]
    """
    return _n_array(s)


def big_l_prime_array(p, n):
    """Compile the L' array from pattern ``p`` and its N array.

    >>> big_l_prime_array("abb", n_array("abb"))
    [0, 0, 2]
    """
    return _big_l_prime_array(p, n)


def big_l_array(p, lp):
    """Compile the L array from pattern ``p`` and its L' array."""
    return _big_l_array(p, lp)


def small_l_prime_array(n):
    """Compile the l' array (Gusfield theorem 2.2.4) from the N array."""
    return _small_l_prime_array(n)


def good_suffix_table(p):
    """Return ``(L', L, l')`` -- every table needed for the good-suffix rule.

    >>> lp, l, slp = good_suffix_table("GGTAGGT")
    >>> lp
    [0, 0, 0, 0, 3, 0, 0]
    >>> l
    [0, 0, 0, 0, 3, 3, 3]
    >>> slp
    [7, 3, 3, 3, 3, 0, 0]
    """
    n = _n_array(p)
    lp = _big_l_prime_array(p, n)
    return lp, _big_l_array(p, lp), _small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    """Shift amount from the good-suffix rule given a mismatch at offset ``i``.

    Pass the ``L'`` table for the *strong* rule, or ``L`` for the weak rule.
    """
    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1  # i now points to the leftmost matching position of P
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    """Shift amount from the good-suffix rule after a full match of P in T."""
    return len(small_l_prime) - small_l_prime[1]


def dense_bad_char_tab(p, amap):
    """Build a dense bad-character table, indexed by offset then character.

    ``amap`` maps each alphabet character to its integer column index.
    Characters of ``p`` outside ``amap`` are skipped rather than raising.
    """
    tab = []
    nxt = [0] * len(amap)
    for i in range(len(p)):
        tab.append(nxt[:])
        if p[i] in amap:
            nxt[amap[p[i]]] = i + 1
    return tab


# ---------------------------------------------------------------------------
# Boyer-Moore
# ---------------------------------------------------------------------------
class BoyerMoore(object):
    """Preprocess a pattern for Boyer-Moore matching.

    Builds the bad-character table and the strong good-suffix tables once,
    so the same pattern can be matched against many texts cheaply.

    Parameters
    ----------
    pattern : str
    alphabet : str
        Characters that may appear. Defaults to the DNA alphabet ``ACGT``;
        pass e.g. ``"ACGTN"`` or the amino acids for other data. Characters
        outside the alphabet are handled gracefully (treated as absent).
    """

    def __init__(self, pattern, alphabet="ACGTN"):
        self.pattern = pattern
        self.p = pattern  # alias, for compatibility with course-style code
        self.alphabet = alphabet
        self.amap = {c: i for i, c in enumerate(alphabet)}
        # bad character table: tab[i][c] holds 1 + (last index of c in p[:i])
        self.bad_char = dense_bad_char_tab(pattern, self.amap)
        # good suffix tables
        if len(pattern) > 1:
            self.big_l_prime, self.big_l, self.small_l_prime = \
                good_suffix_table(pattern)
        else:
            self.big_l_prime = [0] * len(pattern)
            self.big_l = [0] * len(pattern)
            self.small_l_prime = [0] * len(pattern)

    def bad_character_rule(self, i, mismatch_char):
        """Skip count from the bad-character rule at pattern offset ``i``."""
        ci = self.amap.get(mismatch_char, -1)
        if ci == -1:  # character not in alphabet: shift past it
            return i + 1
        return i - (self.bad_char[i][ci] - 1)

    def good_suffix_rule(self, i, strong=True):
        """Skip count from the good-suffix rule at mismatch offset ``i``.

        ``strong=True`` (the default) uses the L' table and gives the larger,
        strictly correct shifts. ``strong=False`` uses the weak L table, which
        is what the classic course implementation applies.
        """
        table = self.big_l_prime if strong else self.big_l
        return good_suffix_mismatch(i, table, self.small_l_prime)

    def match_skip(self):
        """Skip count from the good-suffix rule after a full match."""
        if len(self.pattern) < 2:
            return 1
        return good_suffix_match(self.small_l_prime)


def boyer_moore_match(pattern, text, alphabet="ACGTN", count_ops=False):
    """Find all occurrences of ``pattern`` in ``text`` using Boyer-Moore.

    Parameters
    ----------
    pattern, text : str
    alphabet : str
        Symbols the sequences may contain (default DNA ``ACGTN``).
    count_ops : bool
        If True also return ``(alignments_tried, character_comparisons)``.

    Returns
    -------
    list[int]                       (count_ops=False)
    tuple[list[int], int, int]      (count_ops=True)

    >>> boyer_moore_match("GGTAGGT", "GGTAGGTGGTAGGT")
    [0, 7]
    """
    if pattern == "" or len(pattern) > len(text):
        return ([], 0, 0) if count_ops else []

    bm = BoyerMoore(pattern, alphabet)
    occurrences = []
    n, m = len(text), len(pattern)
    i = 0
    alignments = 0
    comparisons = 0
    while i <= n - m:
        alignments += 1
        shift = 1
        mismatched = False
        for j in range(m - 1, -1, -1):  # compare right to left
            comparisons += 1
            if pattern[j] != text[i + j]:
                skip_bc = bm.bad_character_rule(j, text[i + j])
                skip_gs = bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            shift = max(shift, bm.match_skip())
        i += shift
    if count_ops:
        return occurrences, alignments, comparisons
    return occurrences


# ---------------------------------------------------------------------------
# Strand-aware matching
# ---------------------------------------------------------------------------
def match_both_strands(pattern, text, matcher=None, alphabet="ACGTN"):
    """Find ``pattern`` on both strands of a double-stranded DNA ``text``.

    The plain matchers (:func:`naive_match`, :func:`boyer_moore_match`) are
    string algorithms with no notion of strands -- they only search the
    sequence exactly as given. In DNA the pattern may instead occur on the
    complementary strand, which this function covers by additionally
    searching for the reverse complement of the pattern.

    All offsets are reported in coordinates of the **given** ``text``, so no
    conversion is needed. Searching for ``revcomp(pattern)`` in ``text`` is
    equivalent to searching for ``pattern`` in ``revcomp(text)``, but keeps
    the coordinates straightforward.

    Parameters
    ----------
    pattern, text : str
    matcher : callable, optional
        Exact matcher to use; defaults to :func:`boyer_moore_match`. Any
        function with the signature ``f(pattern, text)`` returning offsets
        works, e.g. :func:`naive_match`.
    alphabet : str
        Passed through to Boyer-Moore when it is the matcher.

    Returns
    -------
    list[tuple[int, int]]
        Sorted ``(offset, strand)`` pairs, where strand is ``+1`` for the
        given sequence and ``-1`` for the complementary strand. A
        palindromic pattern (one equal to its own reverse complement) is
        reported once per site, on the ``+1`` strand only.

    >>> match_both_strands("AGCATG", "TTTTCATGCTTTT")
    [(4, -1)]
    >>> match_both_strands("GGTACC", "AAAAGGTACCAAAA")
    [(4, 1)]
    """
    from .sequence import reverse_complement

    if matcher is None:
        def matcher(p, t):
            return boyer_moore_match(p, t, alphabet=alphabet)

    forward = matcher(pattern, text)
    results = [(offset, 1) for offset in forward]

    rc = reverse_complement(pattern)
    if rc != pattern:  # palindrome: identical hits, do not report twice
        seen = set(forward)
        for offset in matcher(rc, text):
            if offset not in seen:
                results.append((offset, -1))

    return sorted(results)
