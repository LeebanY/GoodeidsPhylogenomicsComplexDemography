#!/bin/bash
#SBATCH --job-name=getbams    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=16gb                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=indexbams.log   # Standard output and error log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ly36@st-andrews.ac.uk
pwd; hostname; date

REF=~/scratch/Gmultiradiatus_genome/G_multiradiatus_assembly_and_annotation/GMassembly.fasta 
loc=/home/lyusuf/scratch/Biogeography_Goodeids/FASTQC

conda activate samtoolsENV
#conda activate freebayesENV

bwa index $REF

#bwa mem $REF $loc/CB_1.fq.gz $loc/CB_2.fq.gz > CB.sam
#bwa mem $REF $loc/Ga_1.fq.gz $loc/Ga_2.fq.gz > GA.sam

#Borealis 75
samtools view -bS CB.sam > CB.bam
samtools view -bS GA.sam > GA.bam


samtools sort CB.bam -o CB.srt.bam
samtools sort GA.bam -o GA.srt.bam


samtools rmdup -s CB.srt.bam CB.srt.rmdup.bam
samtools rmdup -s GA.srt.bam GA.srt.rmdup.bam

samtools view -b -F 4 CB.srt.rmdup.bam > CB.mapped.bam
samtools view -b -F 4 GA.srt.rmdup.bam > GA.mapped.bam

mv *.mapped.bam mapped_bams/
#rm *.bam
#rm *.sam


picard AddOrReplaceReadGroups \
I=mapped_bams/CB.mapped.bam  \
O=CB.final.bam  \
RGID=CB \
RGLB=goodeids \
RGPL=illumina \
RGPU=unit1 \
RGSM=CB

picard AddOrReplaceReadGroups \
I=mapped_bams/GA.mapped.bam \
O=GA.final.bam \
RGID=GA \
RGLB=goodeids \
RGPL=illumina \
RGPU=unit1 \
RGSM=GA  

samtools index mapped_bams/CB.final.bam
#samtools index mapped_bams/AS.final.bam
#samtools index mapped_bams/AT.final.bam
#samtools index mapped_bams/CL.final.bam
samtools index mapped_bams/GA.final.bam
#samtools index mapped_bams/IF.final.bam
#samtools index mapped_bams/XC.final.bam
#samtools index mapped_bams/XR.final.bam



#freebayes-parallel <(fasta_generate_regions.py $ref 1000000) 16 -f $ref --report-genotype-likelihood-max --no-population-priors --use-best-n-alleles 4 --hwe-priors-off --use-mapping-quality --ploidy 2 --theta 0.02 --haplotype-length -1 --genotype-qualities --bam mapped_bams/CB.final.bam --bam mapped_bams/AS.final.bam --bam mapped_bams/AT.final.bam --bam mapped_bams/CL.final.bam --bam mapped_bams/GA.final.bam --bam mapped_bams/IF.final.bam --bam mapped_bams/XC.final.bam --bam mapped_bams/XR.final.bam  > Goodeids.vcf
