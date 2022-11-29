#!/bin/bash
#SBATCH --job-name=vcf2ali    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=16gb                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=vcf2ali.log   # Standard output and error log
pwd; hostname; date

ref=~/scratch/Gmultiradiatus_genome/G_multiradiatus_assembly_and_annotation/GMassembly.fasta
vcf=../Goodeids.GATKHardFilt.vcf.gz 
gff1=bedfiles/GM_Overlap_With_CB_sorted.bed
gff2a=bedfiles/GM_Overlap_SplitFileAtfaultygene.bed
gff2b=bedfiles/GM_Overlap_SplitFileAtfaultygene_2.bed
gff3=bedfiles/GM_Overlap_SplitFileAtfaultygene_3.bed
gff4=bedfiles/GM_Overlap_SplitFileAtfaultygene_4.bed
gff5=bedfiles/GM_Overlap_SplitFileAtfaultygene_5.bed



conda activate vcf2fasta_env

python vcf2fasta/vcf2fasta.py --fasta $ref --vcf $vcf --gff $gff1 --feat CDS --blend --skip --out FASTA_ALIGNMENTS_1
python vcf2fasta/vcf2fasta.py --fasta $ref --vcf $vcf --gff $gff2a --feat CDS --blend --skip --out FASTA_ALIGNMENTS_2a
python vcf2fasta/vcf2fasta.py --fasta $ref --vcf $vcf --gff $gff2b --feat CDS --blend --skip --out FASTA_ALIGNMENTS_2b
python vcf2fasta/vcf2fasta.py --fasta $ref --vcf $vcf --gff $gff3 --feat CDS --blend --skip --out FASTA_ALIGNMENTS_3
python vcf2fasta/vcf2fasta.py --fasta $ref --vcf $vcf --gff $gff4 --feat CDS --blend --skip --out FASTA_ALIGNMENTS_4
python vcf2fasta/vcf2fasta.py --fasta $ref --vcf $vcf --gff $gff5 --feat CDS --blend --skip --out FASTA_ALIGNMENTS_5
