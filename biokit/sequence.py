"""Core sequence operations on DNA/RNA strings.

Everything here works on plain Python ``str`` objects, so there are no
external dependencies. Sequences are treated case-insensitively where it
makes sense, but functions return upper-case results.
"""

from collections import Counter

from .codon_table import CODON_TABLE, COMPLEMENT


def clean(seq):
    """Return ``seq`` upper-cased with surrounding whitespace removed.

    >>> clean("  atgc\\n")
    'ATGC'
    """
    return seq.strip().upper()


def complement(seq):
    """Return the base-by-base complement of a DNA sequence.

    >>> complement("ATGC")
    'TACG'
    """
    return "".join(COMPLEMENT.get(base, "N") for base in seq)


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence.

    >>> reverse_complement("ATGC")
    'GCAT'
    """
    return "".join(COMPLEMENT.get(base, "N") for base in reversed(seq))


def transcribe(seq):
    """Transcribe DNA to RNA (replace T with U).

    >>> transcribe("ATGC")
    'AUGC'
    """
    return clean(seq).replace("T", "U")


def reverse_transcribe(seq):
    """Reverse-transcribe RNA back to DNA (replace U with T).

    >>> reverse_transcribe("AUGC")
    'ATGC'
    """
    return clean(seq).replace("U", "T")


def translate(seq, to_stop=False):
    """Translate a DNA sequence into a protein string.

    The sequence is read in non-overlapping triplets from position 0.
    Stop codons are rendered as ``*``. Trailing bases that do not form a
    full codon are ignored.

    Parameters
    ----------
    seq : str
        DNA sequence (T, not U). RNA is accepted and converted.
    to_stop : bool
        If True, translation halts at (and excludes) the first stop codon.

    >>> translate("ATGGCCTGA")
    'MA*'
    >>> translate("ATGGCCTGA", to_stop=True)
    'MA'
    """
    seq = clean(seq).replace("U", "T")
    protein = []
    for i in range(0, len(seq) - 2, 3):
        amino = CODON_TABLE.get(seq[i:i + 3], "X")
        if to_stop and amino == "*":
            break
        protein.append(amino)
    return "".join(protein)


def gc_content(seq):
    """Return the GC content of a sequence as a fraction in [0, 1].

    Returns 0.0 for an empty sequence.

    >>> round(gc_content("GGCCATAT"), 3)
    0.5
    """
    seq = clean(seq)
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return gc / len(seq)


def base_composition(seq):
    """Return a dict of base -> count for the sequence.

    >>> base_composition("AATGC") == {"A": 2, "T": 1, "G": 1, "C": 1}
    True
    """
    return dict(Counter(clean(seq)))


def melting_temp(seq):
    """Rough melting temperature (Wallace rule) in degrees Celsius.

    Uses Tm = 2*(A+T) + 4*(G+C). A quick estimate best suited to short
    primers (< ~14 nt); use nearest-neighbour models for anything serious.

    >>> melting_temp("ATGC")
    12
    """
    seq = clean(seq)
    at = seq.count("A") + seq.count("T")
    gc = seq.count("G") + seq.count("C")
    return 2 * at + 4 * gc
