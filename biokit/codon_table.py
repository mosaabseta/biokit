"""Standard genetic code (NCBI translation table 1).

Provides the DNA codon -> amino acid mapping used by
:func:`biokit.sequence.translate`, plus the sets of start and stop codons.
"""

# Standard genetic code. Keys are DNA codons (T, not U). "*" marks a stop.
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

START_CODONS = frozenset({"ATG"})
STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})

# Complement lookup for reverse-complement (upper and lower case).
COMPLEMENT = {
    "A": "T", "T": "A", "G": "C", "C": "G", "N": "N",
    "a": "t", "t": "a", "g": "c", "c": "g", "n": "n",
}
