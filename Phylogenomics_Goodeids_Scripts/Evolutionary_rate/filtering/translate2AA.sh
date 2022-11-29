#!/bin/bash
#SBATCH --job-name=masce2translate    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=macse2.log   # Standard output and error log
pwd; hostname; date



conda activate vcf2fasta_env

for all in FASTA_ALIGNMENTS_1_CDS/exported_fasta/*;do
fbname=$(basename "$all" .fas)
java -jar macse_v2.05.jar -prog translateNT2AA -seq $all -out_AA $fbname.AA.fas
sed -i s'/!/-/g' $fbname.AA.fas
pal2nal.pl $fbname.AA.fas $all -nogap -nomismatch -output fasta > $fbname.aligned.pal2nal.fasta
done 
