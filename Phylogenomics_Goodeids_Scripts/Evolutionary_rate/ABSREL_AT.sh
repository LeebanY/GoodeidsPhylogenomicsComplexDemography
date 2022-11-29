#!/bin/bash
#SBATCH --job-name=absrel    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=20gb                     # Job memory request
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=long
#SBATCH --output=absrel.log   # Standard output and error log
pwd; hostname; date

conda activate HYPHY

location=/home/lyusuf/scratch/Biogeography_Goodeids/Mapping/Evolutionary_Rate/FINAL_FASTA
absrel_tree=/home/lyusuf/scratch/Biogeography_Goodeids/Mapping/Evolutionary_Rate/at_absrel_tree.txt

for fas in $location/*.fasta;do 
genename=$(basename "$fas" .fasta)
hyphy absrel --alignment $fas --tree $absrel_tree --CPU 16 --branches Foreground --output ATOWERI_TEST_absrel/$genename.absrel.JSON.txt
done
