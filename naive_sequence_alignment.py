
import matplotlib.pyplot as plt

# Naive Exact Matching Algorithm for Sequence Alignment
def naive(p,t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        match = True
        for j in range(len(p)):
            if t[i+j] != p[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    return occurrences

# Reverse Complement Function
def reverse_complement(seq):
    complements = {'A':'T','C':'G','G':'C','T':'A', 'N':'N'}
    rev_comp = ''
    for base in seq:
        rev_comp = complements[base] + rev_comp
    return rev_comp

# Naive Exact Matching Algorithm, Including Reverse Complements
def naive_rc(p,t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        match = True
        for j in range(len(p)):
            if t[i+j] != p[j]:
                match = False
                break
        if match:
            occurrences.append(i)

    rev_p = reverse_complement(p)
    if not rev_p == p:
        p = rev_p
        for i in range(len(t) - len(p) + 1):
            match = True
            for j in range(len(p)):
                if t[i+j] != p[j]:
                    match = False
                    break
            if match:
                occurrences.append(i)
    
    return occurrences

# Naive Matching Algorithm, Allowing up to 2 Mismatches
def naive_2mm(p,t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        match = True
        mismatches = 0
        for j in range(len(p)):
            if t[i+j] != p[j]:
                match = False
                if mismatches == 3:
                    break
                else:
                    mismatches += 1
        if match or mismatches <= 2:
            occurrences.append(i)
    return occurrences







