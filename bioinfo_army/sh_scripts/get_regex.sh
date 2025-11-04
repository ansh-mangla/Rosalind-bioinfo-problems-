# run the meme
meme ../data_bioinfro_army/rosalind_meme.txt

# get the regular expression 
python get_regex.py ./meme_out/meme.xml > ans_rosalind_meme.txt

# move the file to the coreect location
mv ans_rosalind_meme.txt ../answer_bioinfo_army

# remove the unwanted files
rm -r meme_out
