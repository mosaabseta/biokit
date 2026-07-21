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
