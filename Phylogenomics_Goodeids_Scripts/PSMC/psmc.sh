#!/bin/bash
#SBATCH --job-name=psmc    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=40G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=psmc.log   # Standard output and error log
pwd; hostname; date

Sample_List="AS GM GA CB IF AT CL XC XR"
vcf=Goodeids.GATKHardFilt.vcf.gz
ref=~/scratch/Gmultiradiatus_genome/G_multiradiatus_assembly_and_annotation/GMassembly.fasta
#conda activate vcftools_env
conda activate PSMC 
#conda install bcftools=1.3.1

#for i in $Sample_List;do
#samtools mpileup -Q 30 -q 30 -u -v -f $ref final_bams/$i.final.bam | bcftools call -c | vcfutils.pl vcf2fq -d 10 -D 100 -Q 30 > $i.fq.gz
#PSMC/psmc/utils/fq2psmcfa -q20 $i.fq.gz > $i.psmcfa
#PSMC/psmc/psmc -p "4+25*2+4+6" -o $i.psmc $i.psmcfa
#PSMC/psmc/utils/psmc_plot.pl -u 4.89e-8 -g 1 -p $i.plot $i.psmc
#done 

#Generation time = 1 based on Yoli's observations of Goodeid generation time in the wild.
#Mutation rate from Poecilia reticulata mutation rate estimate:originally estimated in eznick DN, Shaw FH, Rodd FH, Shaw RG. Evaluation of the Rate of Evolution in Natural Populations of Guppies (Poecilia reticulata). Science. 1997;275(5308):1934–7
#But used in The Genome of the Trinidadian Guppy, Poecilia reticulata, and Variation in the Guanapo Population (Kunster et al 2016.) for PSMC analysis. 

PSMC/psmc/psmc -p "4+25*2+4+6" -o IF.psmc IF.psmcfa

samtools mpileup -Q 30 -q 30 -u -v -f $ref ../GM_GENOMIC_READ_DATA/fGirMul1_S1_merged.final.bam | bcftools call -c | vcfutils.pl vcf2fq -d 10 -D 100 -Q 30 > GM.fq.gz
PSMC/psmc/utils/fq2psmcfa -q20 GM.fq.gz > GM.psmcfa
PSMC/psmc/psmc -p "4+25*2+4+6" -o GM.psmc GM.psmcfa
