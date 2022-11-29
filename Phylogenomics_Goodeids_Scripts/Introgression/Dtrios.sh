#!/bin/bash
#SBATCH --job-name=dsuiteinvest    # Job name
#SBATCH --ntasks=4                    # Run on a single CPU
#SBATCH --partition=long
#SBATCH --output=dsuiteinvest.log   # Standard output and error log
pwd; hostname; date

VCF=../Goodeids.GATKHardFilt.vcf.gz
sets=SETS.txt

conda activate vcftools_env

./Build/Dsuite Dinvestigate -w 250,100 ../Goodeids.GATKHardFilt.NoMissingGenotypes.vcf.gz SETS_Dinvestigate.txt GoodeidTestTrio.txt 
