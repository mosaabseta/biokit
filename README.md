# biokit

A small, **dependency-free** bioinformatics toolkit in pure Python. It bundles
the everyday building blocks — FASTA/FASTQ I/O, sequence manipulation, ORF
finding, k-mer counting, sequence distances — together with two exact
pattern-matching algorithms implemented from scratch: **naive** and
**Boyer–Moore** (bad-character + strong good-suffix rules).

The core library imports nothing outside the Python standard library, so it is
easy to read, easy to vendor, and easy to teach from. The only optional
dependency is Biopython, used solely by the `blast` module.

---

## Contents

| Module            | What it gives you                                                        |
|-------------------|--------------------------------------------------------------------------|
| `biokit.io`       | `parse_fasta`, `iter_fasta`, `write_fasta`, `parse_fastq`, `phred_scores`|
| `biokit.sequence` | `reverse_complement`, `transcribe`, `translate`, `gc_content`, …          |
| `biokit.orf`      | `find_orfs`, `longest_orf`, `longest_orf_in_frame` (all six frames)       |
| `biokit.kmers`    | `count_kmers`, `most_frequent_kmer`, `all_max_frequency_kmers`            |
| `biokit.matching` | `naive_match`, `BoyerMoore`, `boyer_moore_match`, plus the Gusfield preprocessing tables (`z_array`, `n_array`, `good_suffix_table`, …) |
| `biokit.distance` | `hamming_distance`, `edit_distance`                                       |
| `biokit.stats`    | `summarize_fasta`, `length_distribution`                                  |
| `biokit.blast`    | `blast_sequence` — optional NCBI BLAST (needs Biopython + network)        |

---

## Installation

### Option A — install as a package (recommended)

```bash
git clone https://github.com/mosaabseta/biokit.git
cd biokit
pip install .
```

Editable install for development:

```bash
pip install -e ".[dev]"      # includes pytest
```

Add the optional BLAST support:

```bash
pip install ".[blast]"       # pulls in biopython
```

### Option B — just copy the folder

The `biokit/` folder is self-contained. Drop it next to your script and
`import biokit`. No installation needed.

**Requirements:** Python 3.7+. No third-party packages for the core library.

---

## Quick start

```python
import biokit as bk

# 1. Read a FASTA file into {id: sequence}
records = bk.parse_fasta("sequences.fasta")

# 2. High-level summary
print(bk.summarize_fasta("sequences.fasta"))
# {'n_records': 18, 'longest_id': '...', 'longest_length': 4805, ...}

# 3. Sequence operations
seq = records["seq1"]
print(bk.reverse_complement(seq))
print(bk.translate(seq, to_stop=True))     # protein up to first stop
print(bk.gc_content(seq))                   # 0.0–1.0

# 4. Longest ORF (searches all 6 reading frames)
orf = bk.longest_orf(seq)
print(orf.strand, orf.frame, orf.start, orf.end, orf.length)
print(bk.translate(orf.sequence))

# 5. Most frequent k-mer (repeat / motif detection)
kmer, freq = bk.most_frequent_kmer(list(records.values()), k=12)

# 6. Exact matching
print(bk.naive_match("ATG", seq))           # [0, 15, 43, ...]
print(bk.boyer_moore_match("ATG", seq))     # same result, fewer comparisons
```

Run the bundled demo end-to-end:

```bash
python examples/demo.py examples/sample.fasta
```

---

## Exact matching: naive vs Boyer–Moore

Both matchers return a list of **0-based** offsets where the pattern occurs,
and both give identical results — Boyer–Moore just gets there with far fewer
character comparisons.

```python
from biokit import naive_match, boyer_moore_match

text = "GGTAGGTGGTAGGT"
naive_match("GGTAGGT", text)                 # [0, 7]
boyer_moore_match("GGTAGGT", text)           # [0, 7]
```

### Counting the work each one does

Pass `count_ops=True` to get `(offsets, alignments_tried, comparisons)`:

```python
offsets, aligns, comps = boyer_moore_match("GCTAGCTA", genome, count_ops=True)
```

On a 400 kb toy genome, Boyer–Moore typically performs **~20% of the character
comparisons** the naive scan does for the same query.

### The `BoyerMoore` class (preprocess once, search many)

If you match the same pattern against many texts, build the tables once:

```python
from biokit import BoyerMoore, boyer_moore_match

# Default alphabet is DNA "ACGT". For proteins or other data, pass your own:
boyer_moore_match(pattern, protein_text, alphabet="ACDEFGHIKLMNPQRSTVWY")
```

**Implementation note.** Boyer–Moore combines two shift rules and takes the
larger shift at each mismatch:

- **Bad-character rule** — on a mismatch, shift so the offending text
  character lines up with its rightmost occurrence in the pattern (or past it).
- **Strong good-suffix rule** — reuse the suffix that already matched. The
  tables (`L'`, `L`, `l'`) are built from the pattern's Z/N arrays following
  Gusfield's *Algorithms on Strings, Trees and Sequences*.

Correctness is enforced by tests that cross-check Boyer–Moore against the naive
matcher on tens of thousands of random inputs (see below).

### Inspecting the preprocessing tables

All the Gusfield tables are public, so you can inspect them, teach from them,
or build your own matcher on top:

```python
from biokit import (
    z_array, n_array, big_l_prime_array, big_l_array, small_l_prime_array,
    good_suffix_table, good_suffix_mismatch, good_suffix_match,
    dense_bad_char_tab,
)

z_array("abababab")            # [8, 0, 6, 0, 4, 0, 2, 0]
n_array("abracadabra")         # [1, 0, 0, 4, 0, 1, 0, 1, 0, 0, 11]

lp, l, slp = good_suffix_table("GGTAGGT")
# lp  (L')  [0, 0, 0, 0, 3, 0, 0]
# l   (L)   [0, 0, 0, 0, 3, 3, 3]
# slp (l')  [7, 3, 3, 3, 3, 0, 0]

good_suffix_mismatch(5, lp, slp)   # 7  -- strong rule
good_suffix_mismatch(5, l,  slp)   # 4  -- weak rule
```

The `BoyerMoore` object carries all three tables as attributes
(`.big_l_prime`, `.big_l`, `.small_l_prime`, `.bad_char`) and its
`good_suffix_rule(i, strong=True)` lets you switch between the strong (`L'`)
and weak (`L`) rule:

```python
bm = BoyerMoore("GGTAGGT")          # default alphabet "ACGTN"
bm.bad_character_rule(6, "A")       # 3
bm.good_suffix_rule(5, strong=True) # 7
bm.good_suffix_rule(5, strong=False)# 4
bm.match_skip()
```

The matcher uses the **strong** rule by default, which shifts at least as far
as the weak rule at every offset and never skips a valid occurrence.

---

## ORF finding

```python
from biokit import find_orfs, longest_orf, longest_orf_in_frame

orfs = find_orfs(seq)                 # every ORF, both strands
orfs = find_orfs(seq, both_strands=False)   # forward strand only

best = longest_orf(seq)               # an ORF namedtuple, or None
# ORF(strand=+1|-1, frame=0|1|2, start, end, length, sequence)

# Longest ORF restricted to one forward reading frame:
n = longest_orf_in_frame(seq, frame=1)
```

An ORF runs from a start codon (`ATG`) to the first in-frame stop codon
(`TAA`/`TAG`/`TGA`), inclusive. For reverse-strand hits, positions index into
the reverse complement and `sequence` is read 5'→3' on that strand.

---

## Running the tests

```bash
pip install pytest
pytest                       # or: python tests/test_biokit.py
```

The suite includes randomized fuzz tests asserting that `boyer_moore_match`
returns exactly the same offsets as `naive_match` across thousands of random
sequences and alphabets — the safety net that keeps the fast matcher honest.

---

## Project layout

```
biokit/
├── biokit/            # the importable package (pure standard library)
│   ├── __init__.py    # public API re-exports
│   ├── io.py
│   ├── sequence.py
│   ├── orf.py
│   ├── kmers.py
│   ├── matching.py    # naive + Boyer–Moore
│   ├── distance.py
│   ├── stats.py
│   ├── blast.py       # optional (Biopython)
│   └── codon_table.py
├── examples/
│   ├── demo.py
│   └── sample.fasta
├── tests/
│   ├── test_biokit.py       # core library tests + BM-vs-naive fuzzing
│   └── test_bm_preproc.py   # Boyer–Moore table reference vectors
├── pyproject.toml
├── requirements.txt
├── LICENSE
└── README.md
```

---

## Attribution

The Boyer–Moore preprocessing in `biokit.matching` implements the Z / N / L' /
L / l' table constructions from:

> Dan Gusfield, *Algorithms on Strings, Trees and Sequences: Computer Science
> and Computational Biology*, Cambridge University Press, 1997
> (theorems 1.4.1, 2.2.2 and 2.2.4).

The public function names (`z_array`, `n_array`, `big_l_prime_array`,
`good_suffix_table`, `dense_bad_char_tab`, …) intentionally mirror the
`bm_preproc.py` teaching module by **Ben Langmead**, distributed with the
Johns Hopkins *Algorithms for DNA Sequencing* course, so that code and
exercises written against that module drop in unchanged. The expected values
in `tests/test_bm_preproc.py` are the same published reference vectors used
there. The implementations in this repository were written independently and
are MIT-licensed as part of this project; credit for the interface design and
the excellent teaching material belongs to Ben Langmead.

If you redistribute this project, please keep this section intact.

---

## License

MIT — see [LICENSE](LICENSE). Use it, fork it, teach with it.
