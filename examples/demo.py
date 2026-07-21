"""End-to-end demo of biokit.

Run from the repository root with:

    python examples/demo.py examples/sample.fasta

If no path is given it falls back to the bundled sample file.
"""

import os
import sys

# Make the package importable when running the script directly from a clone.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import biokit as bk


def main(path):
    print("=" * 60)
    print("FASTA summary")
    print("=" * 60)
    summary = bk.summarize_fasta(path)
    for key, value in summary.items():
        if isinstance(value, float):
            value = round(value, 4)
        print("  {:<16} {}".format(key, value))

    records = bk.parse_fasta(path)

    print("\n" + "=" * 60)
    print("Longest ORF across all records (both strands)")
    print("=" * 60)
    best_overall = None
    for rid, seq in records.items():
        orf = bk.longest_orf(seq)
        if orf and (best_overall is None or orf.length > best_overall[1].length):
            best_overall = (rid, orf)
    if best_overall:
        rid, orf = best_overall
        print("  record : {}".format(rid))
        print("  strand : {:+d}   frame: {}".format(orf.strand, orf.frame))
        print("  span   : {}-{}   length: {}".format(orf.start, orf.end, orf.length))
        print("  protein: {}".format(bk.translate(orf.sequence)))

    print("\n" + "=" * 60)
    print("Most frequent 6-mer (repeat detection)")
    print("=" * 60)
    kmer, freq = bk.most_frequent_kmer(list(records.values()), 6)
    print("  {} occurs {} times".format(kmer, freq))
    max_freq, winners = bk.all_max_frequency_kmers(list(records.values()), 6)
    print("  {} distinct 6-mers tie at the max frequency of {}".format(
        len(winners), max_freq))

    print("\n" + "=" * 60)
    print("Exact matching: naive vs Boyer-Moore")
    print("=" * 60)
    text = next(iter(records.values()))
    pattern = "ATG"
    naive_hits, na, nc = bk.naive_match(pattern, text, count_ops=True)
    bm_hits, ba, bc = bk.boyer_moore_match(pattern, text, count_ops=True)
    print("  pattern {!r} in first record".format(pattern))
    print("  naive       -> {} hits, {} alignments, {} comparisons".format(
        len(naive_hits), na, nc))
    print("  boyer-moore -> {} hits, {} alignments, {} comparisons".format(
        len(bm_hits), ba, bc))
    print("  identical results:", naive_hits == bm_hits)

    print("\n" + "=" * 60)
    print("Boyer-Moore preprocessing tables for 'GGTAGGT'")
    print("=" * 60)
    lp, l, slp = bk.good_suffix_table("GGTAGGT")
    print("  L' (strong) : {}".format(lp))
    print("  L  (weak)   : {}".format(l))
    print("  l'          : {}".format(slp))
    print("  shift on mismatch at offset 5 -> strong {}, weak {}".format(
        bk.good_suffix_mismatch(5, lp, slp),
        bk.good_suffix_mismatch(5, l, slp)))


if __name__ == "__main__":
    default = os.path.join(os.path.dirname(__file__), "sample.fasta")
    main(sys.argv[1] if len(sys.argv) > 1 else default)
