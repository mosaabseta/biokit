"""Open Reading Frame (ORF) detection.

An ORF here is defined the standard way: a run beginning at a start codon
(ATG) and ending at the first in-frame stop codon (TAA/TAG/TGA), inclusive
of both. Positions are reported as 0-based indices into the sequence that
was searched.
"""

from collections import namedtuple

from .codon_table import START_CODONS, STOP_CODONS
from .sequence import clean, reverse_complement

# strand: +1 (given sequence) or -1 (reverse complement)
# frame:  0, 1 or 2, the offset the reading started at
# start/end: 0-based indices INTO THE SEARCHED STRAND
ORF = namedtuple("ORF", ["strand", "frame", "start", "end", "length", "sequence"])


def _find_in_strand(seq, strand):
    """Yield ORFs found on a single strand across its three frames."""
    for frame in range(3):
        pos = frame
        n = len(seq)
        while pos < n - 2:
            if seq[pos:pos + 3] in START_CODONS:
                # walk forward to the first in-frame stop codon
                for end in range(pos, n - 2, 3):
                    if seq[end:end + 3] in STOP_CODONS:
                        orf_seq = seq[pos:end + 3]
                        yield ORF(strand, frame, pos, end + 2,
                                  len(orf_seq), orf_seq)
                        # continue scanning after this ORF's stop codon
                        pos = end + 3
                        break
                else:
                    # no stop codon before the end of the sequence
                    pos = n
                    break
            else:
                pos += 3


def find_orfs(seq, both_strands=True):
    """Return every ORF in ``seq`` as a list of :class:`ORF` records.

    Parameters
    ----------
    seq : str
        DNA sequence.
    both_strands : bool
        If True, also search the reverse-complement strand.

    Notes
    -----
    For reverse-strand ORFs, ``start``/``end`` are indices into the reverse
    complement, and ``sequence`` is the ORF as read 5'->3' on that strand.
    """
    seq = clean(seq)
    orfs = list(_find_in_strand(seq, +1))
    if both_strands:
        orfs.extend(_find_in_strand(reverse_complement(seq), -1))
    return orfs


def longest_orf(seq, both_strands=True):
    """Return the single longest ORF in ``seq``, or ``None`` if there is none.

    >>> longest_orf("ATGAAATAG", both_strands=False).sequence
    'ATGAAATAG'
    """
    orfs = find_orfs(seq, both_strands=both_strands)
    if not orfs:
        return None
    return max(orfs, key=lambda o: o.length)


def longest_orf_in_frame(seq, frame):
    """Length of the longest ORF found only in a single forward frame.

    ``frame`` is 0, 1 or 2. Useful for the classic coursework question of
    "the longest ORF in reading frame N".
    """
    seq = clean(seq)
    best = 0
    pos = frame
    n = len(seq)
    while pos < n - 2:
        if seq[pos:pos + 3] in START_CODONS:
            for end in range(pos, n - 2, 3):
                if seq[end:end + 3] in STOP_CODONS:
                    best = max(best, end + 3 - pos)
                    pos = end + 3
                    break
            else:
                break
        else:
            pos += 3
    return best
