
# Homework - Implementing local edit distance calculation function

def read_genome(file):
    genome = ''
    for line in file:
        if not line[0] == ">":
            genome += line.rstrip()
    return genome

def editDistance(x, y):
    # Create distance matrix
    D = []
    for i in range(len(x)+1):
        D.append([0]*(len(y)+1))
    # Initialize first row and column of matrix
    for i in range(len(x)+1):
        D[i][0] = i
    for i in range(len(y)+1):
        D[0][i] = i
    # Fill in the rest of the matrix
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if x[i-1] == y[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    # Edit distance is the value in the bottom right corner of the matrix
    return D[-1][-1]

# Local edit distance
def editDistanceLocal(p,t):
    D = []
    for i in range(len(p)+1):
        D.append([0]*(len(t)+1))
    for i in range(len(p)+1):
        D[i][0] = i
    for i in range(1, len(p)+1):
        for j in range(1, len(t)+1):
            distHor = D[i][j-1] + 1
            distVer = D[i-1][j] + 1
            if p[i-1] == t[j-1]:
                distDiag = D[i-1][j-1]
            else:
                distDiag = D[i-1][j-1] + 1
            D[i][j] = min(distHor, distVer, distDiag)
    return min(D[-1])

genome = read_genome(open("C:\\Users\\taraj\\Downloads\\chr1.GRCh38.excerpt (1).fasta"))

print(editDistanceLocal("GATTTACCAGATTGAG", genome))



