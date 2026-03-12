
from itertools import permutations

# This function computes the overlap length between 2 strings
def overlap(a,b, min_length):
    start = 0

    while True:
        start = a.find(b[:min_length],start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

# This function creates a map of overlaps between reads in a given list
def naive_overlap_map(reads, k):
    olaps = {}
    for a,b in permutations(reads, 2):
        print("a=",a,",b=",b)
        olen = overlap(a,b,k)
        if olen > 0:
            olaps[(a,b)] = olen
    return olaps

reads = ['ACGGTAGATC', 'GATCAAGT', 'TTCACGGA']
print(naive_overlap_map(reads, 3))
