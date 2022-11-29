#!/bin/bash
#SBATCH --job-name=bpp    # Job name
#SBATCH --ntasks=16                    # Run on a single CPU
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G                     # Job memory request
#SBATCH --partition=long
#SBATCH --output=bpp.log   # Standard output and error log
pwd; hostname; date

ctl=Goodeids.ctl.file

#src/bpp --cfile Goodeids.ctl.file
#src/bpp --cfile A00.bpp.ctl
src/bpp --cfile A00.2.bpp.ctl 
