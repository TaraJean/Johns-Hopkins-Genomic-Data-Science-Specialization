
# Homework Assignment: Implementing and Adjusting Boyer_Moore Algorithm, Comparing to Naive Matching

import sys
sys.path.append("C:\\Windows")
import bm_preproc

def read_genome(file):
    genome = ''
    for line in file:
        if not line[0] == ">":
            genome += line.rstrip()
    return genome

def boyer_moore(p, p_bm, t):
    i = 0
    characters_compared = 0
    alignments_tried = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        alignments_tried += 1
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            characters_compared += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences, characters_compared, alignments_tried

def naive(p,t):
    characters_compared = 0
    alignments_tried = 0
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        alignments_tried += 1
        match = True
        for j in range(len(p)):
            characters_compared += 1
            if t[i+j] != p[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    return occurrences, characters_compared, alignments_tried

t = read_genome(open("C:\\Users\\taraj\\Downloads\\chr1.GRCh38.excerpt.fasta"))

p = "GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG"

_,characters, alignments = naive(p,t)
print(characters)

p_bm = bm_preproc.BoyerMoore(p)

_,characters, alignments = boyer_moore(p, p_bm, t)
print(alignments)








