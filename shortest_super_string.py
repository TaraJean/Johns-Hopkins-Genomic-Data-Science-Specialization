
import itertools

# Function finding the length of the overlap between 2 strings
def overlap(a,b, min_length):
    start = 0

    while True:
        start = a.find(b[:min_length],start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

# Function finding the shortest common superstring of substrings by brute force
def scs_brute(subs):
    shortest_sup = None
    for subs_perm in itertools.permutations(subs):
        sup = subs_perm[0]
        for i in range(len(subs)-1):
            olen = overlap(subs_perm[i], subs_perm[i+1], min_length=1)
            sup += subs_perm[i+1][olen:]
        if shortest_sup is None or len(sup) < len(shortest_sup):
            shortest_sup = sup
    return shortest_sup

# Function finding the 2 reads with the greatest overlap from a set of reads
def pick_max_overlap(reads,k):
    reada, readb = None, None
    best_olen = 0
    for a,b in itertools.permutations(reads, 2):
        olen = overlap(a, b, k)
        if olen > best_olen:
            reada, readb = a,b
            best_olen = olen
    return reada, readb, best_olen

# Function finding the shortest common superstring of substrings using a greedy algorithm
def scs_greedy(reads, k):
    read_a, read_b, olen = pick_max_overlap(reads,k)
    while olen > 0:
        reads.remove(read_a)
        reads.remove(read_b)
        reads.append(read_a + read_b[olen:])
        read_a, read_b, olen = pick_max_overlap(reads,k)
    return ''.join(reads)

# Function converting to a De Bruijn graph
def de_bruijn_ize(st, k):
    edges = []
    nodes = set()
    for i in range(len(st)-k+1):
        edges.append((st[i:i+k-1], st[i+1:i+k]))
        nodes.add(st[i:i+k-1])
        nodes.add(st[i+1:i+k])
    return nodes, edges

    





