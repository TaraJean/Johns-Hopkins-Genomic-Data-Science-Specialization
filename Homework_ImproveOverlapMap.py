
# Homework - Improve time complexity of previous overlap map function and explore results

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

def overlap(a,b, min_length):
    start = 0

    while True:
        start = a.find(b[:min_length],start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def overlap_map(reads, k):
    k_mer_sets = {}
    for read in reads:
        for i in range(len(read)-k+1):
            k_mer = read[i:i+k]
            if k_mer in k_mer_sets:
                k_mer_sets[k_mer].add(read)
            else:
                k_mer_sets[k_mer] = {read}

    olaps = {}
    for read in reads:
        suffix = read[len(read)-k:]
        if suffix in k_mer_sets:
            for match in k_mer_sets[suffix]:
                if read != match:
                    olen = overlap(read,match,k)
                    if olen > 0:
                        olaps[(read,match)] = olen

    return olaps

seqs, _ = readFastq(open("C:\\Users\\taraj\\Downloads\\ERR266411_1.for_asm.fastq"))

print(len(overlap_map(seqs, 30)))
olaps = overlap_map(seqs, 30)

nodes_with_outgoing = set()

for read1, read2 in olaps.keys():
    nodes_with_outgoing.add(read1)

print(len(nodes_with_outgoing))
