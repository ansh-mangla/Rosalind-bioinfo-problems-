from function import *
from Bio import SeqIO

# count nucleotide
# seq = open("data/rosalind_dna.txt", 'r').read()
# nucleotide_counter(seq)

# replace t by u
# seq = open("data/rosalind_rna.txt", 'r').read()
# replacing_tu(seq)

# get reverse compiment
# seq = open("data/rosalind_revc.txt", 'r').read()
# reverse_compliment(seq)

# recurrent relation
# values = open("data/rosalind_fib.txt", 'r').read()
# n, k = map(int, values.split(" "))
# find_rabbit_pop(n, k)

# gc contant
# sequences = read_fasta("data/rosalind_gc.txt")
# get_max_gc(sequences)

# pint mutation
# s, t = open("data/rosalind_hamm.txt", "r").readlines()
# point_mutation(s, t)

# mendle's fist law
with open("data/rosalind_iprb.txt", 'r') as file:
    k, m, n = map(int, file.read().strip().split(" "))
dom_pheno(k, m, n)

# dom_pheno(k, m, n)
