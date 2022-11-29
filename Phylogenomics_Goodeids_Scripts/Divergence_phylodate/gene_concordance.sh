#!/bin/bash
#SBATCH --job-name=IQTREE_Gcf   # Job name
#SBATCH --ntasks=4                    # Run on a single CPU
#SBATCH --mem=2gb                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=gcf.log   # Standard output and error log
pwd; hostname; date

conda activate iqtree_env

iqtree -t Goodeid_gene_trees.contact.treefile.OUTPUTWITHOUTGROUPSPECIFIED.txt --gcf Goodeid_gene_trees.contact.treefile -p /home/lyusuf/scratch/Biogeography_Goodeids/Mapping/Evolutionary_Rate/FINAL_FASTA --scf 100 --prefix Goodeid.GCF
