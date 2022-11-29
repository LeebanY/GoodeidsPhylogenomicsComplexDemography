conda activate vcf2fasta_env

for all in AminoAcids/*;do
mafft $all > AminoAcids_Maaft/$all
done

for all in AminoAcids_Maaft/*;do
trimal 
