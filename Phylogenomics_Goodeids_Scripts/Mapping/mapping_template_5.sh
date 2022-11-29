#!/bin/bash
#SBATCH --job-name=getbams5    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=16gb                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=indexbams5.log   # Standard output and error log
pwd; hostname; date

REF=~/scratch/Gmultiradiatus_genome/G_multiradiatus_assembly_and_annotation/GMassembly.fasta 
loc=/home/lyusuf/scratch/Biogeography_Goodeids/FASTQC

conda activate samtoolsENV
#conda activate freebayesENV

#bwa index $REF

#bwa mem $REF $loc/AS_R1.fq.gz $loc/AS_R2.fq.gz > AS.sam
bwa mem -t 16 $REF $loc/AT_R1.fq.gz $loc/AT_R2.fq.gz > AT_2.sam
bwa mem -t 16 $REF $loc/CL_R1.fq.gz $loc/CL_R2.fq.gz > CL.sam

#Borealis 75
samtools view -bS AS.sam > AS.bam
samtools view -bS AT_2.sam > AT.bam
samtools view -bS CL.sam > CL.bam


samtools sort AS.bam -o AS.srt.bam
samtools sort AT.bam -o AT.srt.bam
samtools sort CL.bam -o CL.srt.bam


samtools rmdup -s AS.srt.bam AS.srt.rmdup.bam
samtools rmdup -s AT.srt.bam AT.srt.rmdup.bam
samtools rmdup -s CL.srt.bam CL.srt.rmdup.bam

samtools view -b -F 4 AS.srt.rmdup.bam > AS.mapped.bam
samtools view -b -F 4 AT.srt.rmdup.bam > AT.mapped.bam
samtools view -b -F 4 CL.srt.rmdup.bam > CL.mapped.bam

mv *.mapped.bam mapped_bams/
#rm *.bam
#rm *.sam


picard AddOrReplaceReadGroups \
I=mapped_bams/AS.mapped.bam \
O=AS.final.bam \
RGID=AS \
RGLB=goodeids \
RGPL=illumina \
RGPU=unit1 \
RGSM=AS

picard AddOrReplaceReadGroups \
I=mapped_bams/AT.mapped.bam \
O=AT.final.bam \
RGID=AT \
RGLB=goodeids \
RGPL=illumina \
RGPU=unit1 \
RGSM=AT  

picard AddOrReplaceReadGroups \
I=mapped_bams/CL.mapped.bam \
O=CL.final.bam \
RGID=CL \
RGLB=goodeids \
RGPL=illumina \
RGPU=unit1 \
RGSM=CL

samtools index mapped_bams/AS.final.bam
samtools index mapped_bams/AT.final.bam
samtools index mapped_bams/CL.final.bam



