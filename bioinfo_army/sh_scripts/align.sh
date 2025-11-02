echo "To run this script make sure you execute it in the same folder where it is at" 
# extract the GenBank IDs
id1=$(cut -d " " -f1 ../data_bioinfro_army/rosalind_need.txt)
id2=$(cut -d " " -f2 ../data_bioinfro_army/rosalind_need.txt)

# Fethc the sequence
efetch -db nucleotide -id $id1 -format fasta >> seq1.fasta
efetch -db nucleotide -id $id2 -format fasta >> seq2.fasta

# alig
needle -asequence seq1.fasta -bsequence seq2.fasta \
       -gapopen 10 -gapextend 1 \
       -endweight -endopen 10 -endextend 1 \
       -outfile output.pair

# get the score and save it at the right location
score=$(grep "Score" output.pair | cut -d " " -f3 | head -n 1)
int=${score%.*}
echo "$int" >> ans_rosalind_need.txt
mv ans_rosalind_need.txt ../answer_bioinfo_army

# remove files
rm seq1.fasta seq2.fasta output.pair