"""Optional NCBI BLAST helper.

This is the only module with an external dependency (Biopython) and it
requires a network connection to reach NCBI. It is imported lazily so the
rest of :mod:`biokit` works with zero dependencies.

Install the extra with::

    pip install biokit-seq[blast]
"""


def blast_sequence(sequence, program="blastn", database="nt", top_n=5):
    """Run a remote NCBI BLAST search and return the best hits.

    Parameters
    ----------
    sequence : str
        Query nucleotide or protein sequence.
    program : str
        BLAST program, e.g. ``"blastn"``, ``"blastp"``, ``"blastx"``.
    database : str
        Target database, e.g. ``"nt"`` or ``"nr"``.
    top_n : int
        Number of hits to return, sorted by ascending E-value.

    Returns
    -------
    list[dict]
        Each dict has ``title``, ``length``, ``e_value``, ``score`` and a
        short ``alignment`` preview.

    Raises
    ------
    ImportError
        If Biopython is not installed.
    """
    try:
        from Bio.Blast import NCBIWWW, NCBIXML
    except ImportError as exc:  # pragma: no cover - depends on optional dep
        raise ImportError(
            "blast_sequence requires Biopython. Install it with "
            "'pip install biokit-seq[blast]' or 'pip install biopython'."
        ) from exc

    handle = NCBIWWW.qblast(program, database, sequence)
    hits = []
    for record in NCBIXML.parse(handle):
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                hits.append({
                    "title": alignment.title,
                    "length": alignment.length,
                    "e_value": hsp.expect,
                    "score": hsp.score,
                    "alignment": hsp.query[:75] + "...",
                })
    hits.sort(key=lambda h: h["e_value"])
    return hits[:top_n]
