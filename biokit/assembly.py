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

    ``overlap(s, s)`` returns ``len(s)`` -- a string trivially overlaps itself
    fully. That is correct for a raw string operation; when building an overlap
    graph you almost always want to skip identical reads, which
    :func:`overlap_all_pairs` does. Note also that a *periodic* read (e.g.
    ``"ATATAT"``) can produce short overlaps with other reads that reflect its
    internal repeat rather than a genuine suffix/prefix join; resolving those
    needs coverage or paired-end information beyond this function.

    >>> overlap("TTACGT", "CGTGTGC")
    3
    >>> overlap("TTACGT", "GTGTGC")
    0
    """
    if min_length < 1:
        raise ValueError("min_length must be at least 1")
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
        ``(a, b, overlap_length)`` for every ordered pair of *distinct* read
        strings whose suffix/prefix overlap is at least ``k``. Each ordered
        pair appears at most once, even if a read string is supplied more
        than once in ``reads``.
    """
    # De-duplicate read strings: repeats would otherwise produce spurious
    # self-edges and duplicate tuples. dict.fromkeys preserves first-seen order.
    reads = list(dict.fromkeys(reads))

    index = defaultdict(set)
    for read in reads:
        for i in range(len(read) - k + 1):
            index[read[i:i + k]].add(read)

    edges = []
    seen = set()
    for a in reads:
        suffix = a[-k:]
        for b in index[suffix]:
            if a == b:
                continue
            if (a, b) in seen:
                continue
            olen = overlap(a, b, min_length=k)
            if olen > 0:
                edges.append((a, b, olen))
                seen.add((a, b))
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


# ---------------------------------------------------------------------------
# Exact shortest common superstring
# ---------------------------------------------------------------------------
def scs(reads):
    """Shortest common superstring of ``reads``, found by brute force.

    Tries every ordering of the reads, merges each with maximum overlap, and
    keeps the shortest result. This is guaranteed optimal but runs in O(n!)
    time -- finding the SCS is NP-hard -- so it is only usable for a handful
    of reads. :func:`greedy_scs` is the practical heuristic.

    Returns ``None`` for an empty input.

    >>> scs(["ACGGTACGAGC", "GAGCTTCGGA", "GACACGG"])
    'GACACGGTACGAGCTTCGGA'
    """
    import itertools

    reads = list(reads)
    if not reads:
        return None

    shortest = None
    for perm in itertools.permutations(reads):
        superstring = perm[0]
        for i in range(len(reads) - 1):
            olen = overlap(perm[i], perm[i + 1], min_length=1)
            superstring += perm[i + 1][olen:]
        if shortest is None or len(superstring) < len(shortest):
            shortest = superstring
    return shortest


def scs_all(reads):
    """Return *every* shortest common superstring, not just one.

    The SCS is frequently not unique: several distinct orderings can tie at
    the minimum length. :func:`scs` returns an arbitrary one of them, which
    is why two correct implementations can disagree. This returns the full
    sorted list of tied solutions.

    Returns
    -------
    (length, superstrings) : tuple[int, list[str]]

    >>> length, sols = scs_all(["ABC", "BCA", "CAB"])
    >>> length
    5
    >>> sols
    ['ABCAB', 'BCABC', 'CABCA']
    """
    import itertools

    reads = list(reads)
    if not reads:
        return 0, []

    found = {}
    for perm in itertools.permutations(reads):
        superstring = perm[0]
        for i in range(len(reads) - 1):
            olen = overlap(perm[i], perm[i + 1], min_length=1)
            superstring += perm[i + 1][olen:]
        found.setdefault(len(superstring), set()).add(superstring)

    best = min(found)
    return best, sorted(found[best])


# ---------------------------------------------------------------------------
# Assembling a genome from reads
# ---------------------------------------------------------------------------
def assemble_greedy_contigs(reads, min_overlap=3):
    """Assemble ``reads`` into one or more contigs by greedy merging.

    Like :func:`greedy_scs`, but instead of concatenating whatever cannot be
    merged, the unmergeable pieces are returned separately. That is what real
    assembly produces: a set of **contigs** (contiguous assembled stretches),
    not a single sequence, because coverage gaps and repeats break the graph.

    Returns
    -------
    list[str]
        Contigs, longest first.

    >>> seq = "ACGTTGCATTGCAAGGCTA"
    >>> reads = [seq[i:i+8] for i in range(len(seq)-7)]
    >>> assemble_greedy_contigs(reads, 4) == [seq]
    True
    """
    reads = list(dict.fromkeys(reads))  # duplicates add nothing here
    a, b, olen = _pick_max_overlap(reads, min_overlap)
    while olen > 0:
        reads.remove(a)
        reads.remove(b)
        reads.append(a + b[olen:])
        a, b, olen = _pick_max_overlap(reads, min_overlap)
    return sorted(reads, key=len, reverse=True)


def assemble_from_reads(reads, k=None, min_overlap=3, method="de_bruijn"):
    """Assemble a sequence from ``reads``, the practical entry point.

    Parameters
    ----------
    reads : sequence[str]
    k : int, optional
        k-mer size for the de Bruijn method. Defaults to
        ``min(len(read)) - 1``, capped at 31.
    min_overlap : int
        Minimum overlap for the greedy method.
    method : {"de_bruijn", "greedy"}
        ``"de_bruijn"`` builds a k-mer graph and walks an Eulerian path --
        fast, and the basis of real short-read assemblers, but it needs a
        complete k-mer spectrum and returns ``None`` if none exists.
        ``"greedy"`` merges by largest overlap and always returns contigs,
        but is a heuristic and slower.

    Returns
    -------
    str or list[str] or None
        A single sequence for ``"de_bruijn"`` (or ``None`` if no Eulerian
        path exists), a list of contigs for ``"greedy"``.

    Notes
    -----
    Both methods are strand-blind: a real sequencing run yields reads from
    both strands, so reverse-complement the reads (or their partners) before
    assembling genuine data.
    """
    reads = [r for r in reads if r]
    if not reads:
        return None if method == "de_bruijn" else []

    if method == "greedy":
        return assemble_greedy_contigs_indexed(reads, min_overlap)

    if method != "de_bruijn":
        raise ValueError("method must be 'de_bruijn' or 'greedy'")

    if k is None:
        k = min(min(len(r) for r in reads) - 1, 31)
    if k < 2:
        raise ValueError("reads are too short to assemble (need k >= 2)")
    return assemble_de_bruijn(reads, k)


# ---------------------------------------------------------------------------
# Index-accelerated greedy assembly
# ---------------------------------------------------------------------------
def assemble_greedy_contigs_indexed(reads, min_overlap=30):
    """Greedy assembly accelerated with a k-mer index and a priority queue.

    Produces the same contigs as :func:`assemble_greedy_contigs` but without
    the cubic rescan. :func:`_pick_max_overlap` re-compares every pair after
    every merge, costing O(n^3) overlap calls; here the overlap graph is built
    once with a k-mer index, edges are kept in a max-heap, and after each merge
    only the edges touching the new node are recomputed.

    Parameters
    ----------
    reads : sequence[str]
    min_overlap : int
        Minimum overlap to merge on; also the k-mer size for the index. For
        100 bp reads, 30-40 is a sensible range -- low values (< 10) admit
        chance overlaps between unrelated reads and corrupt the assembly.

    Returns
    -------
    list[str]
        Contigs, longest first.

    Notes
    -----
    Merging is greedy and therefore a heuristic: it does not guarantee the
    shortest common superstring, and results can vary with tie-breaking.

    >>> seq = "ACGTTGCATTGCAAGGCTA"
    >>> reads = [seq[i:i+8] for i in range(len(seq)-7)]
    >>> assemble_greedy_contigs_indexed(reads, 4) == [seq]
    True
    """
    import heapq

    k = min_overlap
    if k < 1:
        raise ValueError("min_overlap must be at least 1")

    nodes = {}          # id -> sequence
    for i, r in enumerate(dict.fromkeys(reads)):
        if len(r) >= k:
            nodes[i] = r
    if not nodes:
        return sorted(dict.fromkeys(reads), key=len, reverse=True)

    index = defaultdict(set)          # k-mer -> ids of nodes containing it

    def add_to_index(node_id):
        seq = nodes[node_id]
        for i in range(len(seq) - k + 1):
            index[seq[i:i + k]].add(node_id)

    def drop_from_index(node_id, seq):
        for i in range(len(seq) - k + 1):
            index[seq[i:i + k]].discard(node_id)

    for node_id in nodes:
        add_to_index(node_id)

    heap = []

    def push_edges_for(node_id):
        """Queue edges out of and into ``node_id`` against live partners."""
        seq = nodes[node_id]
        # outgoing: partners whose prefix matches this node's suffix
        for other in list(index.get(seq[-k:], ())):
            if other == node_id or other not in nodes:
                continue
            olen = overlap(seq, nodes[other], min_length=k)
            if olen > 0:
                heapq.heappush(heap, (-olen, node_id, other))
        # incoming: partners whose suffix matches this node's prefix
        for other in list(index.get(seq[:k], ())):
            if other == node_id or other not in nodes:
                continue
            olen = overlap(nodes[other], seq, min_length=k)
            if olen > 0:
                heapq.heappush(heap, (-olen, other, node_id))

    for node_id in list(nodes):
        push_edges_for(node_id)

    next_id = max(nodes) + 1 if nodes else 0

    while heap:
        neg_olen, a_id, b_id = heapq.heappop(heap)
        # lazy deletion: skip edges whose endpoints were already merged away
        if a_id not in nodes or b_id not in nodes or a_id == b_id:
            continue
        a_seq, b_seq = nodes[a_id], nodes[b_id]
        # the stored length may be stale; re-verify before merging
        olen = overlap(a_seq, b_seq, min_length=k)
        if olen != -neg_olen or olen == 0:
            if olen > 0:
                heapq.heappush(heap, (-olen, a_id, b_id))
            continue

        drop_from_index(a_id, a_seq)
        drop_from_index(b_id, b_seq)
        del nodes[a_id]
        del nodes[b_id]

        merged_id = next_id
        next_id += 1
        nodes[merged_id] = a_seq + b_seq[olen:]
        add_to_index(merged_id)
        push_edges_for(merged_id)

    leftovers = [r for r in dict.fromkeys(reads) if len(r) < k]
    return sorted(list(nodes.values()) + leftovers, key=len, reverse=True)


def greedy_scs_indexed(reads, min_overlap=30):
    """Index-accelerated counterpart of :func:`greedy_scs`.

    Same greedy merging strategy and same single-superstring output, but
    driven by the k-mer index and heap used in
    :func:`assemble_greedy_contigs_indexed` instead of rescanning every pair
    after each merge. Anything that cannot be merged is concatenated, exactly
    as :func:`greedy_scs` does.

    Note that ``min_overlap`` doubles as the index k-mer size, so unlike
    :func:`greedy_scs` (whose ``k`` defaults to 1) it cannot be set below 1
    and very small values are impractical.

    >>> seq = "ACGTTGCATTGCAAGGCTA"
    >>> reads = [seq[i:i+8] for i in range(len(seq)-7)]
    >>> greedy_scs_indexed(reads, 4) == seq
    True
    """
    return "".join(assemble_greedy_contigs_indexed(reads, min_overlap))


def de_bruijn_contigs(reads, k):
    """Assemble ``reads`` into contigs via maximal non-branching paths.

    The robust de Bruijn counterpart to :func:`assemble_de_bruijn`. Rather
    than demanding a single Eulerian path through the whole graph -- which
    real data rarely admits, so that function returns ``None`` -- this walks
    every **maximal non-branching path** (a run of nodes each with in-degree
    and out-degree exactly 1). Those paths are the *unitigs* that real
    assemblers emit: stretches the graph proves are unambiguous, broken
    wherever a repeat or coverage gap creates a branch.

    Always returns contigs, never ``None``.

    Parameters
    ----------
    reads : sequence[str]
    k : int
        k-mer size. For 100 bp reads, 25-31 is a good range.

    Returns
    -------
    list[str]
        Contigs, longest first.

    >>> seq = "ACGTTGCATTGCAAGGCTA"
    >>> reads = [seq[i:i+8] for i in range(len(seq)-7)]
    >>> de_bruijn_contigs(reads, 7) == [seq]
    True
    """
    # Build from DISTINCT k-mers. de_bruijn_graph emits one edge per k-mer
    # *occurrence*, so overlapping reads create parallel edges that would
    # inflate every degree and make nothing look non-branching. Coverage
    # multiplicity belongs in a separate depth statistic, not in the topology.
    kmers = set()
    for read in reads:
        for i in range(len(read) - k + 1):
            kmers.add(read[i:i + k])
    edges = [(kmer[:-1], kmer[1:]) for kmer in sorted(kmers)]
    if not edges:
        return []

    adj = defaultdict(list)
    indeg = defaultdict(int)
    outdeg = defaultdict(int)
    for left, right in edges:
        adj[left].append(right)
        outdeg[left] += 1
        indeg[right] += 1

    all_nodes = set(indeg) | set(outdeg)

    def one_in_one_out(node):
        return indeg[node] == 1 and outdeg[node] == 1

    paths = []
    # Every maximal non-branching path starts at a node that is not 1-in-1-out
    for node in all_nodes:
        if one_in_one_out(node):
            continue
        for nxt in adj[node]:
            path = [node, nxt]
            while one_in_one_out(nxt):
                nxt = adj[nxt][0]
                path.append(nxt)
            paths.append(path)

    # Whatever remains untouched is an isolated cycle of 1-in-1-out nodes
    covered = set()
    for path in paths:
        covered.update(path)
    for node in all_nodes:
        if node in covered or not one_in_one_out(node):
            continue
        cycle = [node]
        covered.add(node)
        nxt = adj[node][0]
        while nxt != node:
            cycle.append(nxt)
            covered.add(nxt)
            nxt = adj[nxt][0]
        cycle.append(node)
        paths.append(cycle)

    contigs = []
    for path in paths:
        seq = path[0]
        for node in path[1:]:
            seq += node[-1]
        contigs.append(seq)
    return sorted(contigs, key=len, reverse=True)


#: Backwards-compatible alias for :func:`assemble_greedy_contigs_indexed`.
assemble_greedy_indexed = assemble_greedy_contigs_indexed
