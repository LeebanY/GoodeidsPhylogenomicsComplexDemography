#!/bin/bash
#SBATCH --job-name=meme    # Job name
#SBATCH --ntasks=8                    # Run on a single CPU
#SBATCH --mem=4gb                     # Job memory request
#SBATCH --cpus-per-task=1
#SBATCH --partition=long
#SBATCH --output=meme.log   # Standard output and error log
pwd; hostname; date

conda activate HYPHY

location=/home/lyusuf/scratch/Biogeography_Goodeids/Mapping/Evolutionary_Rate/FINAL_FASTA
absrel_tree=/home/lyusuf/scratch/Biogeography_Goodeids/Mapping/Evolutionary_Rate/absrel_spp_tree.txt

for fas in $location/*.fasta;do 
genename=$(basename "$fas" .fasta)
hyphy meme --alignment $fas --tree $absrel_tree --branches Internal --CPU 8 --output meme/$genename.absrel.JSON.txt --pvalue 0.05
done
