#!/bin/bash
#SBATCH --job-name=masce2translate    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=40gb                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=macse2.log   # Standard output and error log
pwd; hostname; date

conda activate vcf2fasta_env

for all in FASTA_ALIGNMENTS_1_CDS/2/*;do
fbname=$(basename "$all" .fas)
java -jar macse_v2.05.jar -prog refineAlignment -optim 2 -local_realign_init 0.001 -local_realign_dec 0.001 -align $all -out_NT FASTA_ALIGNMENTS_1_CDS/2/$fbname.NT.fas
java -jar macse_v2.05.jar -prog exportAlignment -align FASTA_ALIGNMENTS_1_CDS/2/$fbname.NT.fas -codonForInternalStop NNN -codonForInternalFS --- -charForRemainingFS - -out_NT FASTA_ALIGNMENTS_1_CDS/2/$fbname.aligned.exported.fas -out_AA FASTA_ALIGNMENTS_1_CDS/2/$fbname.AA.fas
done

rm FASTA_ALIGNMENTS_1_CDS/2/*.AA.fas
rm FASTA_ALIGNMENTS_1_CDS/2/*.NT.fas
