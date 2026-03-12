
import collections
import random
import matplotlib.pyplot as plt

#This function takes two strings, s1 and s2, and returns the longest common prefix between them
def longest_common_prefix(s1,s2):
    i = 0
    while i < len(s1) and i < len(s2) and s1[i] == s2[i]:
        i+= 1
    return s1[:i]

#This function returns the reverse complement of a given DNA sequence
def reverse_complement(seq):
    complements = {'A':'T','C':'G','G':'C','T':'A', 'N':'N'}
    rev_comp = ''
    for base in seq:
        rev_comp = complements[base] + rev_comp
    return rev_comp

#This function converts a .fasta file to a Python string
def read_genome(file):
    genome = ''
    for line in file:
        if not line[0] == ">":
            genome += line.rstrip()
    return genome

#The below two functions convert quality scores to phred33 and vice-versa
def QtoPhred33(Q):
    return chr(Q+33)

def Phred33toQ(qual):
    return ord(qual)-33

#This function converts a .fastq file to two Python lists of strings, sequences and qualities
def readFastq(file):
    sequences = []
    qualities = []
    while True:
        file.readline()
        seq = file.readline().rstrip()
        file.readline()
        qual = file.readline().rstrip()
        if len(seq) == 0:
            break
        else:
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities

# This function creates a histogram of quality scores given in Phred33 format
def createHist(qualities):
    hist = [0] * 50
    for qual in qualities:
        for phred in qual:
            q = Phred33toQ(phred)
            hist[q] += 1
    return hist

# This function returns the GC percent of each base position
def findGCbyPos(reads):
    gc = [0] * len(reads[0])
    totals = [0] * len(reads[0])
    for read in reads:
        for i in range(len(read)):
            if read[i] == 'C' or read[i] == 'G':
                gc[i] += 1
            totals[i] += 1
    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])
    return gc








        


        




