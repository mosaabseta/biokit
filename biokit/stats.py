"""Summary statistics over a set of FASTA records."""

from .io import parse_fasta
from .sequence import gc_content


def summarize_fasta(path):
    """Return a summary dict for a FASTA file.

    Keys: ``n_records``, ``longest_id``, ``longest_length``,
    ``shortest_id``, ``shortest_length``, ``mean_length``,
    ``mean_gc`` (average GC fraction across records).
    """
    records = parse_fasta(path)
    if not records:
        return {"n_records": 0}

    lengths = {rid: len(seq) for rid, seq in records.items()}
    longest_id = max(lengths, key=lengths.get)
    shortest_id = min(lengths, key=lengths.get)
    total_len = sum(lengths.values())
    mean_gc = sum(gc_content(seq) for seq in records.values()) / len(records)

    return {
        "n_records": len(records),
        "longest_id": longest_id,
        "longest_length": lengths[longest_id],
        "shortest_id": shortest_id,
        "shortest_length": lengths[shortest_id],
        "mean_length": total_len / len(records),
        "mean_gc": mean_gc,
    }


def length_distribution(path):
    """Return an ``{id: length}`` dict for every record in a FASTA file."""
    return {rid: len(seq) for rid, seq in parse_fasta(path).items()}
