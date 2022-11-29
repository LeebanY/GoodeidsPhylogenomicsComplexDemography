#!/bin/bash
#SBATCH --job-name=variantfilt    # Job name
#SBATCH --ntasks=8                    # Run on a single CPU
#SBATCH --mem=5gb                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=variantfilt.log   # Standard output and error log
pwd; hostname; date


conda activate vcftools_env

#remove sites with no alternative alleles or only alternative alleles. Also remove non-biallelic positions. Also removing positions that don't have complete genotype information (no missing data).
bcftools view -e 'AC==0 || AC==AN || F_MISSING > 0.2' -m2 -M2 -O z -o Goodeids.GATKHardFilt.NoMissingGenotypes.BiAllelicOnly.vcf.gz Goodeids.GATKHardFilt.NoMissingGenotypes.vcf.gz
bcftools +prune -w 100bp -n 1 -o Goodeids.GATKHardFilt.NoMissingGenotypes.BiAllelicOnly.PrunedEvery100bp.vcf.gz Goodeids.GATKHardFilt.NoMissingGenotypes.BiAllelicOnly.vcf.gz

mv Goodeids.GATKHardFilt.NoMissingGenotypes.BiAllelicOnly.vcf.gz Species_Tree/
