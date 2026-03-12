
# Homework - Implement and Adjust k_mer Indexing algorithm for sequence alignment, using the pigeonhole principle

import bisect

def read_genome(file):
    genome = ''
    for line in file:
        if not line[0] == ">":
            genome += line.rstrip()
    return genome

class Index(object):
    def __init__(self, t, k):
        # Create index from all substrings of size 'length'
        self.k = k
        self.index = []
        for i in range(len(t) - k + 1):
            self.index.append((t[i:i+k], i))
        self.index.sort()

    def query(self, p):
        # Query the index
        n_parts = 3
        len_parts = int(len(p) / n_parts)
        parts = [(p[:len_parts],0),(p[len_parts:2*len_parts], len_parts), (p[2*len_parts:3*len_parts], 2*len_parts)]
        hits = []
        num_hits = 0
        for part,offset in parts:
            kmer = part[:self.k]
            j = bisect.bisect_left(self.index, (kmer, -1))
            while j < len(self.index):
                if self.index[j][0] != kmer:
                    break
                hits.append(self.index[j][1]-offset)
                num_hits += 1
                j += 1
        return set(hits),num_hits

    def queryIndex(self, p, t):
        # Validate Index queries
        k = self.k
        offsets = []
        hits, num_hits = self.query(p)
        for i in hits:
            if i < 0 or i + len(p) > len(t):
                continue
            
            mismatches = 0
            for j in range(len(p)):
                if p[j] != t[i+j]:
                    mismatches += 1
                    if mismatches > 2:
                        break
            if mismatches <= 2:
                offsets.append(i)
        return offsets, num_hits

t = read_genome(open("C:\\Users\\taraj\\Downloads\\chr1.GRCh38.excerpt.fasta"))
p = "GGCGCGGTGGCTCACGCCTGTAAT"

index = Index(t, 8)
offsets, num_hits = index.queryIndex(p, t)
print(num_hits)
