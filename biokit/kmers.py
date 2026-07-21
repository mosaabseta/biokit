"""k-mer counting and repeat/motif discovery.

Counts substrings of a fixed length ``k`` across one or many sequences.
"""

from collections import Counter

from .sequence import clean


def count_kmers(seq, k):
    """Return a ``Counter`` of all length-``k`` substrings in ``seq``.

    >>> count_kmers("ATATA", 2) == {"AT": 2, "TA": 2}
    True
    """
    seq = clean(seq)
    counts = Counter()
    for i in range(len(seq) - k + 1):
        counts[seq[i:i + k]] += 1
    return counts


def count_kmers_multi(sequences, k):
    """Count length-``k`` substrings across an iterable of sequences.

    ``sequences`` may be a list of strings or the ``.values()`` of a
    FASTA dict. Counts are pooled across all sequences.
    """
    counts = Counter()
    for seq in sequences:
        seq = clean(seq)
        for i in range(len(seq) - k + 1):
            counts[seq[i:i + k]] += 1
    return counts


def most_frequent_kmer(sequences, k):
    """Return ``(kmer, frequency)`` for the single most common length-``k`` word.

    Accepts either one sequence (str) or an iterable of sequences.
    """
    if isinstance(sequences, str):
        sequences = [sequences]
    counts = count_kmers_multi(sequences, k)
    if not counts:
        return None, 0
    return counts.most_common(1)[0]


def all_max_frequency_kmers(sequences, k):
    """Return ``(max_freq, [kmers])`` for every word tied at the top frequency.

    Handy for the "how many distinct k-mers occur the maximum number of
    times" style of question.
    """
    if isinstance(sequences, str):
        sequences = [sequences]
    counts = count_kmers_multi(sequences, k)
    if not counts:
        return 0, []
    max_freq = max(counts.values())
    winners = sorted(km for km, freq in counts.items() if freq == max_freq)
    return max_freq, winners
