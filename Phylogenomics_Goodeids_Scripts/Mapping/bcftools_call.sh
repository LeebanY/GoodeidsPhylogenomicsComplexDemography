#!/bin/bash
#SBATCH --job-name=bcftools   # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=16gb                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=bcftools.log   # Standard output and error log
pwd; hostname; date


ref=~/scratch/Gmultiradiatus_genome/G_multiradiatus_assembly_and_annotation/GMassembly.fasta
genes=GM_Overlap_With_CB_sorted.gene.merged.bed
out=Goodeids.bcftools.vcf.gz
out2=Goodeids.bcftools.filtered.vcf.gz
bamlist=bamlist.txt

conda activate vcftools_env

#bcftools mpileup -f $ref -b $bamlist -R $genes | bcftools call -m -Oz -f GQ -o Goodeids.bcftools.vcf.gz

vcftools --gzvcf $out --remove-indels --max-missing 0.8 --min-meanDP 20 --recode --stdout | gzip -c > $out2
