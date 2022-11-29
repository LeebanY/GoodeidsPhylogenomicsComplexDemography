#!/bin/bash
#SBATCH --job-name=LIFTOFF    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=40gb                     # Job memory request
#SBATCH --cpus-per-task=1
#SBATCH --partition=long
#SBATCH --output=liftoff.log   # Standard output and error log
pwd; hostname; date


conda activate liftoff

guppy=~/scratch/FISH_WGA/Guppy_genome/data/GCF_000633615.1/chr_combined_guppy.fas
GM=~/scratch/Gmultiradiatus_genome/G_multiradiatus_assembly_and_annotation/GMassembly.fasta
gff=Gmultiradiatus_genome.augustus.gff
outfile=/home/lyusuf/scratch/Biogeography_Goodeids/Mapping/Evolutionary_Rate/liftoff/Gmultiradiatus_to_guppy.annotations.gff
chromfile=/home/lyusuf/scratch/Biogeography_Goodeids/Mapping/Evolutionary_Rate/liftoff/Gmultiradiatus_to_guppy.chroms.txt
unplacedfile=/home/lyusuf/scratch/Biogeography_Goodeids/Mapping/Evolutionary_Rate/liftoff/Gmultiradiatus_to_guppy.unplacedscafs.txt

#Transfer over annotaton to guppy genome coordinates to be able to plot over chromosomes for all species.
liftoff -g $gff -o $outfile -copies -cds -p 16 $guppy $GM
