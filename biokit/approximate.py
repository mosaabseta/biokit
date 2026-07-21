"""Approximate pattern matching allowing a bounded number of substitutions.

"Approximate" here means Hamming distance: substitutions only, no insertions
or deletions. For indels use :func:`biokit.distance.edit_distance`.
"""

from .matching import boyer_moore_match


def naive_approximate_match(pattern, text, max_mismatches=2, count_ops=False):
    """Find every offset where ``pattern`` occurs with <= ``max_mismatches``.

    Brute force: tries every alignment, counting mismatches and bailing out
    as soon as the budget is blown.

    >>> naive_approximate_match("ACGT", "ACGTACTTAAAA", max_mismatches=2)
    [0, 4]
    """
    occurrences = []
    n, m = len(text), len(pattern)
    alignments = 0
    comparisons = 0
    if m == 0 or m > n:
        return (occurrences, 0, 0) if count_ops else occurrences

    for i in range(n - m + 1):
        alignments += 1
        mismatches = 0
        for j in range(m):
            comparisons += 1
            if text[i + j] != pattern[j]:
                mismatches += 1
                if mismatches > max_mismatches:
                    break
        if mismatches <= max_mismatches:
            occurrences.append(i)

    if count_ops:
        return occurrences, alignments, comparisons
    return occurrences


def approximate_match(pattern, text, max_mismatches=2,
                      alphabet="ACGTN", count_ops=False):
    """Approximate matching via the pigeonhole principle + Boyer-Moore.

    Splits ``pattern`` into ``max_mismatches + 1`` disjoint segments. If the
    pattern occurs with at most ``max_mismatches`` substitutions then at
    least one segment must match *exactly*, since the mismatches can spoil at
    most ``max_mismatches`` of the segments. Exact hits for each segment are
    found with Boyer-Moore, then the full pattern is verified around each hit.

    Parameters
    ----------
    count_ops : bool
        If True, returns ``(offsets, alignments, comparisons, index_hits)``.
        ``index_hits`` counts every exact segment hit that was verified --
        note this is a 4-tuple, unlike the 3-tuple from the exact matchers.

    >>> approximate_match("ACGT", "ACGTACTTAAAA", max_mismatches=2)
    [0, 4]
    """
    m = len(pattern)
    n = len(text)
    if m == 0 or m > n:
        return ([], 0, 0, 0) if count_ops else []

    segment_length = m // (max_mismatches + 1)

    # Too short to split: pigeonhole gives no leverage, so fall back.
    if segment_length < 1:
        result = naive_approximate_match(pattern, text, max_mismatches,
                                         count_ops=count_ops)
        if count_ops:
            hits, a, c = result
            return hits, a, c, 0
        return result

    all_matches = set()
    alignments = 0
    comparisons = 0
    index_hits = 0

    for i in range(max_mismatches + 1):
        start = i * segment_length
        end = min((i + 1) * segment_length, m)
        segment = pattern[start:end]

        hits, a, c = boyer_moore_match(segment, text, alphabet=alphabet,
                                       count_ops=True)
        alignments += a
        comparisons += c

        for hit in hits:
            index_hits += 1
            offset = hit - start          # where the full pattern would start
            if offset < 0 or offset + m > n:
                continue                  # would hang off either end

            mismatches = 0
            for j in range(0, start):     # verify left of the segment
                comparisons += 1
                if pattern[j] != text[offset + j]:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        break
            if mismatches <= max_mismatches:
                for j in range(end, m):   # verify right of the segment
                    comparisons += 1
                    if pattern[j] != text[offset + j]:
                        mismatches += 1
                        if mismatches > max_mismatches:
                            break

            if mismatches <= max_mismatches:
                all_matches.add(offset)

    occurrences = sorted(all_matches)
    if count_ops:
        return occurrences, alignments, comparisons, index_hits
    return occurrences
