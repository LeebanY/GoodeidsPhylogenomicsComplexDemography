#!/bin/bash
#SBATCH --job-name=plink    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=2gb                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=plink.log   # Standard output and error log
pwd; hostname; date

VCF=../Goodeids.GATKHardFilt.NoMissingGenotypes.BiAllelicOnly.PrunedEvery100bp.PASSED.header.vcf.gz
#vcf_test=../Goodeids.GATKHardFilt.NoMissingGenotypes.BiAllelicOnly.PrunedEvery100bp.vcf.gz
conda activate plink

plink --vcf $VCF --threads 16 --double-id --allow-extra-chr --set-missing-var-ids @:# --vcf-half-call m --pca --out Goodeids

#while read line;do
#plink --vcf $VCF --chr $line --threads 16 --double-id --allow-extra-chr --set-missing-var-ids @:# --vcf-half-call m --pca --out Goodeid.$line.PCA
#done < Guppychr.txt
