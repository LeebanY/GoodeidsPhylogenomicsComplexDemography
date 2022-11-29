#!/bin/bash
#SBATCH --job-name=FitMG94    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=16gb                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=FitMG94.log   # Standard output and error log
pwd; hostname; date

conda activate HYPHY

location=/home/lyusuf/scratch/Biogeography_Goodeids/Mapping/Evolutionary_Rate/FINAL_FASTA
treefile=unannotated_tree.txt
comp=/home/lyusuf/scratch/FISH_WGA/Genomes/WGA_CYPRINODON/MafFilter/HYPHY/hyphy-develop/res/TemplateBatchFiles/FitMG94.bf

#This script will estimate dn/ds for each branch in the tree for all 7000 genes that have been successfully aligned.
for fas in $location/*.fasta;do
alignment=$(basename "$fas" .fasta) 
hyphy $comp --alignment $fas --tree $treefile --type global --output FitMG94_global/$alignment.fitmg94.txt
done
