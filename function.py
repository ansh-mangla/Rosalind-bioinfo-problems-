from Bio import SeqIO


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
    open("answer_file/ans_rosalind_dna.txt", 'w').write(ans)
    return dict

#############################################
# replce t by b


def replacing_tu(seq):
    replaced = seq.replace("T", "U")
    open("answer_file/ans_rosalind_rna.txt", 'w').write(replaced)
    return replaced

#############################################
# reverse compliment


def reverse_compliment(seq: str):
    None
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

    open("answer_file/ans_rosalind_revc.txt", 'w').write(rev_comp)
    return rev_comp

######################################################
# recurrence relation


def find_rabbit_pop(n: int, k: int):
    rabbit = [1, 1]
    for m in range(1, n+1):
        if m > 2:
            r = rabbit[-1] + k*rabbit[-2]
            rabbit.append(r)

    open("answer_file/ans_rosalind_fib.txt", "w").write(str(rabbit[-1]))
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
    open("answer_file/ans_rosalind_gc.txt", 'w').write(string)
    return (max_id, max_gc)
