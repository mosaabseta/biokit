"""Substring and subsequence indexes for fast candidate lookup.

Two index structures, both built once over a text and then queried many
times. Each is a sorted list of ``(key, offset)`` pairs searched with binary
search, so construction is O(n log n) and a query is O(log n + hits).

* :class:`Index` -- indexes contiguous k-mers.
* :class:`SubseqIndex` -- indexes *spaced* subsequences: k characters taken
  every ``ival`` positions. Spacing makes keys more specific than contiguous
  k-mers of the same length, so fewer false candidates survive verification.

The paired ``query_index`` / ``query_subseq`` functions use these for
approximate matching under the pigeonhole principle.
"""

import bisect


class Index(object):
    """An index of every contiguous k-mer in a text.

    >>> idx = Index("GGTATTCGGGA", 3)
    >>> idx.query("GGTA")
    [0]
    """

    def __init__(self, t, k):
        """Build an index of all length-``k`` substrings of ``t``."""
        self.k = k
        self.index = []
        for i in range(len(t) - k + 1):
            self.index.append((t[i:i + k], i))
        self.index.sort()

    def query(self, p):
        """Return offsets where the first k-mer of ``p`` occurs in the text.

        These are *candidates* only -- the caller must verify the rest of the
        pattern at each offset.
        """
        kmer = p[:self.k]
        i = bisect.bisect_left(self.index, (kmer, -1))
        hits = []
        while i < len(self.index):
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits


class SubseqIndex(object):
    """An index of spaced subsequences of a text.

    Holds every subsequence of ``k`` characters taken ``ival`` positions
    apart. ``SubseqIndex("ATAT", 2, 2)`` extracts ``("AA", 0)`` and
    ``("TT", 1)``.

    Attributes
    ----------
    span : int
        Width of text a single subsequence covers: ``1 + ival * (k - 1)``.

    >>> ind = SubseqIndex("ATATAT", 2, 2)
    >>> ind.index
    [('AA', 0), ('AA', 2), ('TT', 1), ('TT', 3)]
    >>> ind.query("ATATAT")
    [0, 2]
    """

    def __init__(self, t, k, ival):
        """Create index from all subsequences of ``k`` characters spaced
        ``ival`` positions apart."""
        self.k = k          # number of characters per subsequence
        self.ival = ival    # spacing; 1 = adjacent, 2 = every other, etc.
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):
            self.index.append((t[i:i + self.span:ival], i))
        self.index.sort()   # alphabetize by subsequence

    def query(self, p):
        """Return index hits for the first subsequence of ``p``.

        To search the other offsets, query with ``p[1:]``, ``p[2:]`` and so
        on -- that is what :func:`query_subseq` does.
        """
        subseq = p[:self.span:self.ival]
        i = bisect.bisect_left(self.index, (subseq, -1))
        hits = []
        while i < len(self.index):
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits


# ---------------------------------------------------------------------------
# Index-driven approximate matching
# ---------------------------------------------------------------------------
def query_index(p, t, index, max_mismatches=2):
    """Approximate matching using a contiguous k-mer :class:`Index`.

    Partitions ``p`` into ``max_mismatches + 1`` disjoint k-mers and looks
    each one up, then verifies the whole pattern at every candidate.

    Returns
    -------
    (offsets, index_hits) : tuple[list[int], int]
        ``offsets`` are the distinct positions matching with at most
        ``max_mismatches`` substitutions; ``index_hits`` is the number of
        candidates the index produced before verification.

    Notes
    -----
    The index must have been built over the same ``t``. Requires
    ``len(p) >= (max_mismatches + 1) * index.k``, otherwise the segments
    cannot be disjoint and the pigeonhole guarantee fails -- a ValueError is
    raised rather than silently returning wrong answers.
    """
    k = index.k
    segments = max_mismatches + 1
    if len(p) < segments * k:
        raise ValueError(
            "pattern of length %d is too short for %d disjoint %d-mers; "
            "use a smaller k or fewer mismatches" % (len(p), segments, k))

    all_matches = set()
    index_hits = 0
    for i in range(segments):
        start = i * k
        for hit in index.query(p[start:]):
            index_hits += 1
            offset = hit - start
            if offset < 0 or offset + len(p) > len(t):
                continue
            mismatches = 0
            for j in range(len(p)):
                if p[j] != t[offset + j]:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        break
            if mismatches <= max_mismatches:
                all_matches.add(offset)
    return sorted(all_matches), index_hits


def query_subseq(p, t, subseq_index, max_mismatches=2):
    """Approximate matching using a :class:`SubseqIndex`.

    Queries the index once per offset ``0 .. max_mismatches``, which selects
    disjoint spaced subsequences of ``p``. Because the subsequences are
    disjoint, at most ``max_mismatches`` of them can be spoiled, so at least
    one must match exactly.

    Returns
    -------
    (offsets, index_hits) : tuple[list[int], int]

    Notes
    -----
    Requires ``subseq_index.ival > max_mismatches`` so that the queried
    subsequences really are disjoint; otherwise a ValueError is raised.

    >>> t = "to-morrow and to-morrow and to-morrow creeps in this petty pace"
    >>> p = "to-morrow and to-morrow "
    >>> ind = SubseqIndex(t, 8, 3)
    >>> occurrences, hits = query_subseq(p, t, ind, 2)
    >>> occurrences
    [0, 14]
    """
    ival = subseq_index.ival
    if ival <= max_mismatches:
        raise ValueError(
            "ival (%d) must exceed max_mismatches (%d) for the subsequences "
            "to be disjoint" % (ival, max_mismatches))

    all_matches = set()
    index_hits = 0
    for i in range(max_mismatches + 1):
        if len(p) - i < subseq_index.span:
            continue  # not enough pattern left to form a subsequence
        for hit in subseq_index.query(p[i:]):
            index_hits += 1
            offset = hit - i
            if offset < 0 or offset + len(p) > len(t):
                continue
            mismatches = 0
            for j in range(len(p)):
                if p[j] != t[offset + j]:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        break
            if mismatches <= max_mismatches:
                all_matches.add(offset)
    return sorted(all_matches), index_hits
