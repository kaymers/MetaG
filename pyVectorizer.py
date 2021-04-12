from Bio import SeqIO
from collections import defaultdict
from itertools import product

all_kmers_cache = {}


def reverse_compliment(s):
        reverse_letter = {
            'A': 'T',
            'C': 'G',
            'G': 'C',
            'T': 'A'
        }
        return ''.join([reverse_letter[i] for i in s[::-1]])

def all_kmers(k):
    if k in all_kmers_cache: return all_kmers_cache[k]
    all_strings = [''.join(i) for i in product("ACGT", repeat=k)]
    proper_strings = [min(i,reverse_compliment(i)) for i in all_strings]
    unique_kmers = list(set(proper_strings))
    unique_kmers.sort()
    all_kmers_cache[k] = unique_kmers
    return unique_kmers

def vectorize_seq(seq, k=4):
    d = defaultdict(lambda: 0)
    for c in range(len(seq)-k+1):
        d[seq[c:c+k]]+=1
    return [d[kmer] if kmer==reverse_compliment(kmer) else d[kmer]+d[reverse_compliment(kmer)]for kmer in all_kmers(k)]

def vectorize_file(path, k):
    file=SeqIO.parse(path, "fasta")
    current_seq=0
    vectors = []
    for i in file:
        vec = vectorize_seq(str(i.seq), k)
        vectors.append(vec)
        if not(current_seq%1000): print("Processing sequence "+str(current_seq))
        current_seq+=1
    return vectors

def count_acgt_file(path):
    file=SeqIO.parse(path, "fasta")
    current_seq=0
    acgt = [0,0,0,0]
    for i in file:
        seq = str(i.seq)
        acgt[0] += (seq.count("A")+seq.count("a"))
        acgt[1] += (seq.count("C")+seq.count("c"))
        acgt[2] += (seq.count("G")+seq.count("g"))
        acgt[3] += (seq.count("T")+seq.count("t"))

        if not(current_seq%1000): print("Processing sequence "+str(current_seq))
        current_seq+=1
    return acgt

#vectors = vectorize_file("./data/AAGA01.1.fsa_nt",4)