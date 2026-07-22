"""De novo assembly primitives: overlaps, overlap graphs, and assembly.

Shotgun sequencing breaks the genome into many short reads with no positional
information. Reassembling them relies on one idea: reads that came from
overlapping stretches of the genome share suffix/prefix overlaps. This module
provides the building blocks:

* :func:`overlap` -- length of the longest suffix/prefix match between two
  reads.
* :func:`overlap_all_pairs` -- the overlap graph edges for a set of reads,
  found efficiently with a k-mer index.
* :func:`greedy_scs` -- greedy shortest-common-superstring assembly.
* :func:`de_bruijn_graph` / :func:`assemble_de_bruijn` -- k-mer (de Bruijn)
  graph assembly, the approach real assemblers are built on.

These are teaching-grade implementations: correct and readable, not tuned for
millions of reads.
"""

from collections import defaultdict


# ---------------------------------------------------------------------------
# Pairwise overlap
# ---------------------------------------------------------------------------
def overlap(a, b, min_length=3):
    """Length of the longest suffix of ``a`` that exactly matches a prefix of ``b``.

    Returns 0 if no overlap of at least ``min_length`` exists. Only proper
    suffix/prefix overlaps are considered (the whole of ``a`` need not match).

    >>> overlap("TTACGT", "CGTGTGC")
    3
    >>> overlap("TTACGT", "GTGTGC")
    0
    """
    start = 0
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:
            return 0
        if b.startswith(a[start:]):            # verify full suffix/prefix match
            return len(a) - start
        start += 1


# ---------------------------------------------------------------------------
# Overlap graph (all pairs) via a k-mer index
# ---------------------------------------------------------------------------
def overlap_all_pairs(reads, k):
    """Return overlap-graph edges for ``reads`` using minimum overlap ``k``.

    Naively this is O(N^2) in the number of reads. Instead, every read is
    indexed by the k-mers it contains; then for each read ``a`` only the
    reads sharing ``a``'s length-``k`` suffix are tested as overlap partners.

    Parameters
    ----------
    reads : sequence[str]
    k : int
        Minimum overlap length (also the k-mer size for the index).

    Returns
    -------
    list[tuple[str, str, int]]
        ``(a, b, overlap_length)`` for every ordered pair ``a != b`` whose
        suffix/prefix overlap is at least ``k``.

    Notes
    -----
    If the same read string appears more than once this compares by value,
    which can produce self-referential-looking edges between the duplicates;
    de-duplicate reads first if that matters for your use.
    """
    index = defaultdict(set)
    for read in reads:
        for i in range(len(read) - k + 1):
            index[read[i:i + k]].add(read)

    edges = []
    for a in reads:
        suffix = a[-k:]
        for b in index[suffix]:
            if a == b:
                continue
            olen = overlap(a, b, min_length=k)
            if olen > 0:
                edges.append((a, b, olen))
    return edges


# ---------------------------------------------------------------------------
# Greedy shortest common superstring
# ---------------------------------------------------------------------------
def _pick_max_overlap(reads, k):
    best_a, best_b, best_olen = None, None, 0
    for i in range(len(reads)):
        for j in range(len(reads)):
            if i == j:
                continue
            olen = overlap(reads[i], reads[j], min_length=k)
            if olen > best_olen:
                best_a, best_b, best_olen = reads[i], reads[j], olen
    return best_a, best_b, best_olen


def greedy_scs(reads, k=1):
    """Assemble reads by greedily merging the pair with the largest overlap.

    Repeatedly finds the two reads with the longest overlap (at least ``k``)
    and merges them, until no overlap of length >= ``k`` remains; anything
    left is concatenated. Returns a single superstring.

    This is a heuristic: it does **not** guarantee the *shortest* common
    superstring (that problem is NP-hard) and the result can depend on how
    ties are broken. It is the standard illustration of greedy assembly.

    >>> greedy_scs(["ABC", "BCA", "CAB"], k=2) in (
    ...     "ABCAB", "BCABC", "CABCA")
    True
    """
    reads = list(reads)
    a, b, olen = _pick_max_overlap(reads, k)
    while olen > 0:
        reads.remove(a)
        reads.remove(b)
        reads.append(a + b[olen:])
        a, b, olen = _pick_max_overlap(reads, k)
    return "".join(reads)


# ---------------------------------------------------------------------------
# De Bruijn graph assembly
# ---------------------------------------------------------------------------
def de_bruijn_graph(reads, k):
    """Build a de Bruijn graph from ``reads`` for k-mer size ``k``.

    Each length-``k`` substring becomes an edge from its (k-1)-prefix node to
    its (k-1)-suffix node.

    Returns
    -------
    (nodes, edges) : tuple[set[str], list[tuple[str, str]]]
        ``nodes`` are the distinct (k-1)-mers; ``edges`` is a list of
        ``(left, right)`` node pairs, one per k-mer occurrence.
    """
    nodes = set()
    edges = []
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i + k]
            left, right = kmer[:-1], kmer[1:]
            nodes.add(left)
            nodes.add(right)
            edges.append((left, right))
    return nodes, edges


def assemble_de_bruijn(reads, k):
    """Reconstruct a sequence from ``reads`` via an Eulerian path of its
    de Bruijn graph.

    Builds the de Bruijn graph for k-mer size ``k`` and walks an Eulerian path
    (every edge used exactly once), stitching the (k-1)-mer nodes into a single
    string. Works when such a path exists -- i.e. the k-mer spectrum is
    complete and the graph is connected with at most the usual two odd-degree
    endpoints. Returns ``None`` if no Eulerian path exists.

    This is the clean, error-free case that motivates why real assemblers use
    de Bruijn graphs; genuine data needs error correction and repeat handling
    on top.

    >>> seq = "ACGTTGCA"
    >>> reads = [seq[i:i+4] for i in range(len(seq)-3)]
    >>> assemble_de_bruijn(reads, 4) == seq
    True
    """
    graph = defaultdict(list)
    outdeg = defaultdict(int)
    indeg = defaultdict(int)

    _, edges = de_bruijn_graph(reads, k)
    for left, right in edges:
        graph[left].append(right)
        outdeg[left] += 1
        indeg[right] += 1

    if not edges:
        return None

    # Find the Eulerian-path start: a node with outdeg - indeg == 1, else any.
    nodes = set(outdeg) | set(indeg)
    start = None
    starts = ends = 0
    for node in nodes:
        d = outdeg[node] - indeg[node]
        if d == 1:
            start = node
            starts += 1
        elif d == -1:
            ends += 1
        elif d != 0:
            return None  # a node is too unbalanced: no Eulerian path
    # A single Eulerian path needs exactly one (start, end) odd pair, or none
    # (an Eulerian circuit). More than one odd pair => ambiguous, give up.
    if not ((starts == 1 and ends == 1) or (starts == 0 and ends == 0)):
        return None
    if start is None:
        start = edges[0][0]

    # Hierholzer's algorithm for the Eulerian path.
    adj = {node: list(dests) for node, dests in graph.items()}
    stack = [start]
    path = []
    while stack:
        node = stack[-1]
        if adj.get(node):
            stack.append(adj[node].pop())
        else:
            path.append(stack.pop())
    path.reverse()

    if len(path) != len(edges) + 1:
        return None  # graph was disconnected; no single Eulerian path

    # Stitch (k-1)-mer nodes: first node in full, then last char of each next.
    result = path[0]
    for node in path[1:]:
        result += node[-1]
    return result
