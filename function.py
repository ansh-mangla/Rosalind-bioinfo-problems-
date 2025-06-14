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
