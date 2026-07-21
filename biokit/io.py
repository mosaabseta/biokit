"""Reading and writing FASTA / FASTQ files, dependency-free.

The FASTA reader mirrors the common pattern of building an ``id -> sequence``
dictionary, but also offers a streaming iterator for large files.
"""


def parse_fasta(path):
    """Read a FASTA file into an ``{id: sequence}`` dict.

    The record id is the first whitespace-delimited token after ``>``
    (the ``>`` itself is stripped). Multi-line sequences are concatenated.

    Parameters
    ----------
    path : str
        Path to a FASTA file.

    Returns
    -------
    dict
        Ordered mapping of record id to its (upper-case) sequence.
    """
    records = {}
    current = None
    with open(path, "r") as handle:
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            elif current is not None:
                records[current].append(line.strip().upper())
    return {rid: "".join(chunks) for rid, chunks in records.items()}


def iter_fasta(path):
    """Yield ``(id, sequence)`` pairs from a FASTA file one at a time.

    Use this instead of :func:`parse_fasta` when the file is too large to
    hold in memory.
    """
    current = None
    chunks = []
    with open(path, "r") as handle:
        for line in handle:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    yield current, "".join(chunks)
                current = line[1:].split()[0]
                chunks = []
            elif current is not None:
                chunks.append(line.strip().upper())
    if current is not None:
        yield current, "".join(chunks)


def write_fasta(records, path, line_width=70):
    """Write an ``{id: sequence}`` mapping to a FASTA file.

    Parameters
    ----------
    records : dict
        Mapping of record id to sequence.
    path : str
        Output path.
    line_width : int
        Wrap sequence lines at this many characters (0 disables wrapping).
    """
    with open(path, "w") as handle:
        for rid, seq in records.items():
            handle.write(">{}\n".format(rid))
            if line_width and line_width > 0:
                for i in range(0, len(seq), line_width):
                    handle.write(seq[i:i + line_width] + "\n")
            else:
                handle.write(seq + "\n")


def parse_fastq(path):
    """Yield ``(id, sequence, quality)`` tuples from a FASTQ file.

    Quality is returned as the raw ASCII-encoded string. Use
    :func:`phred_scores` to convert it to numeric Phred values.
    """
    with open(path, "r") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().rstrip()
            handle.readline()  # the "+" separator line
            qual = handle.readline().rstrip()
            rid = header[1:].split()[0] if header.startswith("@") else header.strip()
            yield rid, seq.upper(), qual


def phred_scores(quality, offset=33):
    """Convert a FASTQ quality string to a list of integer Phred scores.

    ``offset`` is 33 for Sanger / Illumina 1.8+ (the modern default).
    """
    return [ord(ch) - offset for ch in quality]
