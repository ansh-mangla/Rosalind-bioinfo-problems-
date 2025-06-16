from Bio import SeqIO
import os

# for creating file


def save_answer(name, text: str):
    file_path = f"answer_file/ans_rosalind_{name}.txt"
    os.makedirs("answer_file", exist_ok=True)
    with open(file_path, 'w') as file:
        file.write(text)
    print(f"{file_path} has been saved. ")

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
    save_answer("dna", str(ans))
    return dict

#############################################
# replce t by b


def replacing_tu(seq):
    replaced = seq.replace("T", "U")
    save_answer("rna", str(replaced))
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
    save_answer("revc", str(rev_comp))
    return rev_comp

######################################################
# recurrence relation


def find_rabbit_pop(n: int, k: int):
    rabbit = [1, 1]
    for m in range(1, n+1):
        if m > 2:
            r = rabbit[-1] + k*rabbit[-2]
            rabbit.append(r)

    save_answer("fib", str(rabbit[-1]))
    return rabbit[-1]

# read FASTA file


def read_fasta(filepath):
    sequences = []
    for record in SeqIO.parse(filepath, 'fasta'):
        sequences.append((record.id, str(record.seq)))
    return sequences


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
    save_answer("gc", string)
    return (max_id, max_gc)


####################################################
# counting point mutations


def point_mutation(s: str, t: str):
    mutations = 0
    for n in range(0, len(s)):
        if s[n] != t[n]:
            mutations += 1
    print(mutations)
    save_answer("hamm", str(mutations))
    return mutations


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

    save_answer("iprb", str(p))
    return p
