#!/bin/bash
#SBATCH --job-name=getbams4    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=16gb                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=indexbams4.log   # Standard output and error log
pwd; hostname; date

REF=~/scratch/Gmultiradiatus_genome/G_multiradiatus_assembly_and_annotation/GMassembly.fasta 
loc=/home/lyusuf/scratch/Biogeography_Goodeids/FASTQC

conda activate samtoolsENV
#conda activate freebayesENV

#bwa index $REF

bwa mem -t 16 $REF $loc/IF_R1.fq.gz $loc/IF_R2.fq.gz > IF.sam
#bwa mem $REF $loc/XC_1.fq.gz $loc/XC_2.fq.gz > XC.sam
#bwa mem $REF $loc/XR_1.fq.gz $loc/XR_2.fq.gz > XR.sam

#Borealis 75
samtools view -bS IF.sam > IF.bam

samtools sort IF.bam -o IF.srt.bam

samtools rmdup -s IF.srt.bam IF.srt.rmdup.bam

samtools view -b -F 4 IF.srt.rmdup.bam > IF.mapped.bam

mv *.mapped.bam mapped_bams/
#rm *.bam
#rm *.sam


picard AddOrReplaceReadGroups \
I=mapped_bams/IF.mapped.bam \
O=IF.final.bam \
RGID=IF \
RGLB=goodeids \
RGPL=illumina \
RGPU=unit1 \
RGSM=IF

samtools index mapped_bams/IF.final.bam


#freebayes-parallel <(fasta_generate_regions.py $ref 1000000) 16 -f $ref --report-genotype-likelihood-max --no-population-priors --use-best-n-alleles 4 --hwe-priors-off --use-mapping-quality --ploidy 2 --theta 0.02 --haplotype-length -1 --genotype-qualities --bam mapped_bams/CB.final.bam --bam mapped_bams/AS.final.bam --bam mapped_bams/AT.final.bam --bam mapped_bams/CL.final.bam --bam mapped_bams/GA.final.bam --bam mapped_bams/IF.final.bam --bam mapped_bams/XC.final.bam --bam mapped_bams/XR.final.bam  > Goodeids.vcf
