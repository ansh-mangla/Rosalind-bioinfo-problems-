import os
from Bio import SeqIO, Entrez
import numpy as np
import pandas as pd
from bioinfo_stronghold import * # answers from bioinformatics stronghold
from variable import ANSH_EMAIL
# folder to save answers to 
FOLDER = "answer_bioinfo_army"

# Introduction to the Bioinformatics Armory
# nucleotide_counter()

# Gen BAnk Intoduction 
def publication_count(email, genus, start_dt, end_dt, ):
    """This function returs the number of published articles between two dates related to a particular Genus.

    Args:
        email (str): Email for Entrez
        genus (str): Genus for which to search GenBank
        start_dt (str): YYYY/MM/DD start date
        end_dt (str): YYYY/MM/DD end date

    Returns:
        int: number of published articles 
    """
    Entrez.email = email
    handle = Entrez.esearch(db = "nucleotide", term=f"{genus}[organism]",
                        datetype = "pdat", 
                        mindate = start_dt, maxdate=end_dt)
    record = Entrez.read(handle)
    handle.close()
    count = record['Count']
    # print(count)
    save_answer("gbk", str(count), folder=FOLDER)
    return int(count)

# Data Formats

def fetch_fasta(ids):
    """This functions takes in accestion ids and fetch the associeated GenBanak record. Output is both retured as a fasta and a dictionary. 

    Args:
        ids (str): eg "FJ817486 JX069768 JX469983"

    Returns:
        list: [fasta, dict]
        dict -> {id: Seq.Record.obj}
    """
    Entrez.email = ANSH_EMAIL
    handle = Entrez.efetch(db = "nucleotide", id=ids,rettype="fasta", retmode="text")
    fasta_parsed = {record.id : record for record in  SeqIO.parse(handle, 'fasta')}
    handle.close()
    
    fasta_seqs =[]
    for record in fasta_parsed.values():
        fasta_seqs.append(record.format("fasta"))

    fasta_file = "\n".join(fasta_seqs)
    return [fasta_file, fasta_parsed]

def shortest_seq(ids):
    """Takes in a string of GenBank ids, featch the sequences and output the shortest sequence as a fasta. 

    Args:
        ids (str): eg "FJ817486 JX069768 JX469983"

    Returns:
        str: shortest sequence 
    """
    x , seq_dict = fetch_fasta(ids)
    len_seq = {id: len(record.seq) for id, record in seq_dict.items()}
    short_seq_id = min(len_seq, key=len_seq.get)
    seq = seq_dict[short_seq_id].format("fasta")

    save_answer("frmt", seq.strip(), FOLDER)
    return seq

