#!/bin/bash
#SBATCH --job-name=freebayes    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=16gb                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=freebayes_1.log   # Standard output and error log
pwd; hostname; date

source activate freebayesENV
ref=~/scratch/Gmultiradiatus_genome/G_multiradiatus_assembly_and_annotation/GMassembly.fasta

freebayes-parallel <(fasta_generate_regions.py $ref 1000000) 16 -f $ref --report-genotype-likelihood-max --no-population-priors --use-best-n-alleles 4 --hwe-priors-off --use-mapping-quality --ploidy 2 --theta 0.02 --haplotype-length -1 --genotype-qualities --bam final_bams/CB.final.bam --bam final_bams/AS.final.bam --bam final_bams/AT.final.bam --bam final_bams/CL.final.bam --bam final_bams/GA.final.bam --bam final_bams/IF.final.bam --bam final_bams/XC.final.bam --bam final_bams/XR.final.bam --bam ../GM_GENOMIC_READ_DATA/fGirMul1_S1_merged.final.bam > Goodeids.vcf
