#!/bin/bash
#SBATCH --job-name=relax    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=20gb                     # Job memory request
#SBATCH --cpus-per-task=1
#SBATCH --partition=long
#SBATCH --output=relax.log   # Standard output and error log
pwd; hostname; date

conda activate HYPHY


location=/home/lyusuf/scratch/Biogeography_Goodeids/Mapping/Evolutionary_Rate/FINAL_FASTA
absrel_tree=/home/lyusuf/scratch/Biogeography_Goodeids/Mapping/Evolutionary_Rate/absrel_tree.txt

for fas in $location/*.fasta;do 
genename=$(basename "$fas" .fasta)
hyphy RELAX.bf --alignment $fas --tree $absrel_tree > RELAX/$genename.txt
less RELAX/$genename.txt | grep -A20 "### Tests of individual" | sed 's/|//g' | sed 's/-//g' | sed 's/://g' | sed -e 's/    /\t/g' | sed 's/ //' | grep -v "##RELAX" | grep -v "significant" | sed -e 's/\t/ /g' > RELAX/$genename.adjusted.txt
done
