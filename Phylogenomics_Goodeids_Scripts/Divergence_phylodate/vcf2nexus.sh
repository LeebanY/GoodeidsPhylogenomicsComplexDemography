#!/bin/bash
#SBATCH --job-name=vcf2phylip    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --mem=2G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=vcf2phylip.log   # Standard output and error log
pwd; hostname; date

conda activate vcftools_env 

#python vcf2phylip/vcf2phylip.py -i Goodeids.GATKHardFilt.NoMissingGenotypes.BiAllelicOnly.vcf -o CB -m 9 --output-folder Goodeid_SppTree_Nexus -n -b
python vcf2phylip/vcf2phylip.py -i ../Goodeids.GATKHardFilt.NoMissingGenotypes.BiAllelicOnly.PrunedEvery100bp.PASSED.header.vcf -o CB -m 9 --output-folder Goodeid_SppTree_Nexus_Filtered -n -b
