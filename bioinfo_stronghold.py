import os
from Bio import SeqIO
import numpy as np
import pandas as pd

# directory to save answaers to 

# for creating file


def save_answer(name, text: str, folder="answer_file"):
    """
    Saves the output (text) in a file named (name) in the path (folder)
    """
    file_path = f"{folder}/ans_rosalind_{name}.txt"
    os.makedirs(folder, exist_ok=True)
    with open(file_path, 'w') as file:
        file.write(text)
    print(f"{file_path} has been saved. ")

# read FASTA file


def read_fasta(filepath, only_seqs=False):
    sequences = []
    if only_seqs == True:
        for record in SeqIO.parse(filepath, 'fasta'):
            sequences.append(str(record.seq))
        return sequences
    
    sequences = {}
    for record in SeqIO.parse(filepath, 'fasta'):
        sequences[record.id] = str(record.seq)
    return sequences


#############################################
# count number of each nucleotide in a DNA sequence


def nucleotide_counter(seq: str):
    # counter
    dict = {
        "A": seq.count("A"),
        "C": seq.count("C"),
        "G": seq.count("G"),
        "T": seq.count("T"),
    }

    # answer format string
    ans = ""
    for x in dict:
        ans = ans + str(dict[x]) + ' '

    print(dict)  # what is returned
    print(ans)  # the answer
    # save_answer("dna", str(ans))
    return dict

#############################################
# replce t by b


def replacing_tu(seq):
    replaced = seq.replace("T", "U")
    # save_answer("rna", str(replaced))
    return replaced

#############################################
# reverse compliment


def reverse_compliment(seq: str):
    compliment = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }
    rev_comp = ""
    for n in seq:
        if n in compliment:
            rev_comp += str(compliment[n])
    rev_comp = rev_comp[::-1]
    # save_answer("revc", str(rev_comp))
    return rev_comp

######################################################
# recurrence relation


def find_rabbit_pop(n: int, k: int):
    rabbit = [1, 1]
    for m in range(1, n+1):
        if m > 2:
            r = rabbit[-1] + k*rabbit[-2]
            rabbit.append(r)

    # save_answer("fib", str(rabbit[-1]))
    return rabbit[-1]


########################################################
# calculate GC contant


def GC_contant(seq):
    gc_contant = (seq.count("G") + seq.count("C"))/len(seq) * 100
    return gc_contant


def get_max_gc(sequences: list):
    dict = {}
    for seq in sequences:
        dict[seq[0]] = GC_contant(seq[1])

    max_id = max(dict, key=lambda k: dict[k])
    max_gc = max(dict.values())
    string = max_id + "\n" + str(max_gc)
    # save_answer("gc", string)
    return (max_id, max_gc)


####################################################
# counting point mutations


def point_mutation(s: str, t: str):
    mutations = 0
    for n in range(0, len(s)):
        if s[n] != t[n]:
            mutations += 1
    print(mutations)
    # save_answer("hamm", str(mutations))
    return mutations


#####################################################
# mendle's first law


def dom_pheno(k, m, n):
    t = k + m + n
    total_crosses = ((t*(t-1))/2)

    # intercrosses
    ic_AA = (k*(k-1))/2
    ic_Aa = (m*(m-1))/2

    # outccrossesm
    ou_AA_Aa = k*m
    ou_AA_aa = k*n
    ou_Aa_aa = m*n

    p = (ic_AA + (3/4)*(ic_Aa) + ou_AA_Aa +
         ou_AA_aa + (1/2)*(ou_Aa_aa))/total_crosses

    # save_answer("iprb", str(p))
    return p

#####################################################
# mRNA to amino acid


def translte(seq):
    # print(codon_aa)
    df_c2a = pd.read_csv("resource/codon_amino_acid.csv")
    c2a_dict = dict(zip(df_c2a["Codon"], df_c2a["aa"]))
    aa = ''
    for n in range(0, len(seq), 3):
        codon = seq[n:n+3]
        if codon in c2a_dict:
            aa += c2a_dict[codon]
    aa = aa.replace("Stop", "")
    print(aa)
    # save_answer("prot", aa)
    return aa

##################################################
# finding mortis


def find_motifs(seq, mot):
    """
    Finds the motif in a sequence 


    Args:
    seq (str): The sequsence.
    mot (str): the motif.  

    Returns: 
    str: the positions of the motif in the sequence. 

    """
    pos = []
    for p in range(0, len(seq)-len(mot) + 1):
        if seq[p:p+len(mot)] == mot:
            pos.append(p+1)

    ans = " ".join(map(str, pos))
    # save_answer("subs", ans)
    print(ans)
    return pos


##############################################
# finding consensus sequence and profile matrix

def find_coensus_seq(seqs: list):
    len_seq = len(seqs[0])
    c_seq = ""

    # converting sequences to a matrix where each row one seq
    s = np.array([[n for n in seq] for seq in seqs])

    # creating a list of dictionaries contining counts of ATGC in each position on all the strings
    list_of_dicts = []
    for n in range(len_seq):
        # get the unique elemtns in each column and find counts
        nucs, counts = np.unique(s[:, n], return_counts=True)
        # create dicts and add them to a list
        list_of_dicts.append(dict(zip(nucs, counts)))

    # list to df but rows index = ACGT
    profile_m = pd.DataFrame(list_of_dicts).T
    profile_m = profile_m.loc[["A", "C", "G", "T"]]  # reoder
    profile_m = profile_m.fillna(0).astype(int)  # NaN = 0

    # getting conensus seq based on profile mtrix
    for i in range(len_seq):
        m = profile_m.index[profile_m[i] == profile_m[i].max()][0]
        c_seq += str(m)

    # saving answer

    ans_string = profile_m.to_string(header=False)
    for base in "ACGT":
        ans_string = ans_string.replace(base, base+":")

    ans = c_seq + "\n" + ans_string
    print(ans)

    # save_answer("cons", ans)

    return {"profile_matrix": profile_m, "c_seq": c_seq}
