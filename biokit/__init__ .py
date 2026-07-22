"""biokit -- a small, dependency-free bioinformatics toolkit.

Import the functions you need directly from the top level, e.g.::

    from biokit import parse_fasta, translate, longest_orf, boyer_moore_match

Modules
-------
io          FASTA / FASTQ reading and writing
sequence    reverse complement, transcribe, translate, GC content, ...
orf         open reading frame detection (all six frames)
kmers       k-mer counting and repeat/motif discovery
matching    exact matching: naive and Boyer-Moore
approximate approximate matching with k substitutions (pigeonhole)
indexing    k-mer Index and spaced SubseqIndex for fast candidate lookup
distance    Hamming and edit (Levenshtein) distance
stats       FASTA summary statistics
blast       optional NCBI BLAST wrapper (needs Biopython)
"""

from .io import (
    parse_fasta,
    iter_fasta,
    write_fasta,
    parse_fastq,
    phred_scores,
)
from .sequence import (
    clean,
    complement,
    reverse_complement,
    transcribe,
    reverse_transcribe,
    translate,
    gc_content,
    base_composition,
    melting_temp,
)
from .orf import (
    ORF,
    find_orfs,
    longest_orf,
    longest_orf_in_frame,
)
from .kmers import (
    count_kmers,
    count_kmers_multi,
    most_frequent_kmer,
    all_max_frequency_kmers,
)
from .matching import (
    naive_match,
    BoyerMoore,
    boyer_moore_match,
    match_both_strands,
    # Boyer-Moore preprocessing internals, exposed for inspection/teaching
    z_array,
    n_array,
    big_l_prime_array,
    big_l_array,
    small_l_prime_array,
    good_suffix_table,
    good_suffix_mismatch,
    good_suffix_match,
    dense_bad_char_tab,
)
from .approximate import (
    naive_approximate_match,
    approximate_match,
)
from .indexing import (
    Index,
    SubseqIndex,
    query_index,
    query_subseq,
)
from .distance import (
    hamming_distance,
    edit_distance,
)
from .stats import (
    summarize_fasta,
    length_distribution,
)

__version__ = "0.2.0"

__all__ = [
    # io
    "parse_fasta", "iter_fasta", "write_fasta", "parse_fastq", "phred_scores",
    # sequence
    "clean", "complement", "reverse_complement", "transcribe",
    "reverse_transcribe", "translate", "gc_content", "base_composition",
    "melting_temp",
    # orf
    "ORF", "find_orfs", "longest_orf", "longest_orf_in_frame",
    # kmers
    "count_kmers", "count_kmers_multi", "most_frequent_kmer",
    "all_max_frequency_kmers",
    # matching
    "naive_match", "BoyerMoore", "boyer_moore_match", "match_both_strands",
    "z_array", "n_array", "big_l_prime_array", "big_l_array",
    "small_l_prime_array", "good_suffix_table", "good_suffix_mismatch",
    "good_suffix_match", "dense_bad_char_tab",
    # approximate matching
    "naive_approximate_match", "approximate_match",
    # indexing
    "Index", "SubseqIndex", "query_index", "query_subseq",
    # distance
    "hamming_distance", "edit_distance",
    # stats
    "summarize_fasta", "length_distribution",
]
