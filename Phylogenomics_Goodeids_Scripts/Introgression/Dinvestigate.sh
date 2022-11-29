#!/bin/bash
#SBATCH --job-name=dsuiteinvest    # Job name
#SBATCH --ntasks=4                    # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --partition=long
#SBATCH --output=dsuiteinvest.log   # Standard output and error log
pwd; hostname; date

VCF=../Goodeids.GATKHardFilt.NoMissingGenotypes.BiAllelicOnly.PrunedEvery100bp.PASSED.header.vcf.gz
VCF_newcords=../Local_Introgression/Goodeids.Freebayes.GATKFiltered.NomissingGens.Biallelic.CoordinateToNewGM.bcftools.vcf.gz
sets=SETS.txt

conda activate vcftools_env

./Build/Dsuite Dinvestigate -w 250,100 $VCF_newcords SETS_Dinvestigate.txt Trios_Goodeids.txt
