"""Sequence distance and similarity metrics."""


def hamming_distance(a, b):
    """Number of positions at which two equal-length strings differ.

    >>> hamming_distance("GAGCCTACTAACGGGAT", "CATCGTAATGACGGCCT")
    7
    """
    if len(a) != len(b):
        raise ValueError("hamming_distance requires equal-length sequences")
    return sum(1 for x, y in zip(a, b) if x != y)


def edit_distance(a, b):
    """Levenshtein edit distance (min insertions/deletions/substitutions).

    Computed with the classic O(len(a) * len(b)) dynamic-programming table.

    >>> edit_distance("kitten", "sitting")
    3
    """
    m, n = len(a), len(b)
    if m == 0:
        return n
    if n == 0:
        return m
    prev = list(range(n + 1))
    for i in range(1, m + 1):
        curr = [i] + [0] * n
        for j in range(1, n + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            curr[j] = min(
                prev[j] + 1,        # deletion
                curr[j - 1] + 1,    # insertion
                prev[j - 1] + cost  # substitution / match
            )
        prev = curr
    return prev[n]


def edit_distance_recursive(a, b, _memo=None):
    """Edit distance computed with the naive recursion (plus memoization).

    This is the textbook recursive definition: compare the last characters,
    then recurse on the three ways to reconcile a mismatch (delete, insert,
    substitute). Without memoization the recursion is exponential -- the same
    subproblems are recomputed over and over -- which is exactly the argument
    for the dynamic-programming version in :func:`edit_distance`.

    A ``_memo`` dict is threaded through to cache subproblem results, making
    it usable on real inputs; pass nothing and it is created for you. For long
    sequences prefer :func:`edit_distance`, which is iterative and uses O(n)
    memory.

    >>> edit_distance_recursive("kitten", "sitting")
    3
    >>> edit_distance_recursive("GCGTATGC", "TATTGGCTATC")
    7
    """
    if _memo is None:
        _memo = {}
    key = (len(a), len(b))
    if key in _memo:
        return _memo[key]

    if len(a) == 0:
        result = len(b)
    elif len(b) == 0:
        result = len(a)
    else:
        delt = edit_distance_recursive(a[:-1], b, _memo) + 1
        ins = edit_distance_recursive(a, b[:-1], _memo) + 1
        cost = 0 if a[-1] == b[-1] else 1
        sub = edit_distance_recursive(a[:-1], b[:-1], _memo) + cost
        result = min(delt, ins, sub)

    _memo[key] = result
    return result


def edit_distance_matrix(a, b):
    """Return the full edit-distance DP matrix ``D`` for ``a`` (rows) and ``b``.

    ``D[i][j]`` is the edit distance between ``a[:i]`` and ``b[:j]``. The edit
    distance itself is ``D[-1][-1]``. Useful for teaching, for tracing an
    alignment back, or for the approximate-match variant below.
    """
    m, n = len(a), len(b)
    d = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        d[i][0] = i
    for j in range(n + 1):
        d[0][j] = j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost = 0 if a[i - 1] == b[j - 1] else 1
            d[i][j] = min(d[i - 1][j] + 1,        # deletion
                          d[i][j - 1] + 1,        # insertion
                          d[i - 1][j - 1] + cost)  # substitution / match
    return d


def approximate_edit_distance(pattern, text):
    """Fewest edits needed to match ``pattern`` against its best spot in ``text``.

    This adapts the edit-distance DP to approximate matching. The pattern
    labels the rows and the text the columns, and crucially the **first row
    is all zeros** -- a match may begin at any position in the text for free.
    The answer is the minimum value in the bottom row: the edit distance of
    the closest occurrence of ``pattern`` within ``text``.

    Substitutions, insertions and deletions are all allowed (unlike
    :func:`biokit.approximate.approximate_match`, which permits substitutions
    only). Returns just the distance; the location is not tracked.

    >>> approximate_edit_distance("GCGTATGC", "TATTGGCTATCGGCGTATGCAAAA")
    0
    >>> approximate_edit_distance("GCGTATGC", "TATTGGCTATCGGCGTATGGAAAA")
    1
    """
    p, t = pattern, text
    m, n = len(p), len(t)
    d = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        d[i][0] = i          # first column: cost of consuming the pattern
    # first row stays all zeros: a match may start anywhere in the text
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost = 0 if p[i - 1] == t[j - 1] else 1
            d[i][j] = min(d[i - 1][j] + 1,
                          d[i][j - 1] + 1,
                          d[i - 1][j - 1] + cost)
    return min(d[m])         # best match ends somewhere along the bottom row
