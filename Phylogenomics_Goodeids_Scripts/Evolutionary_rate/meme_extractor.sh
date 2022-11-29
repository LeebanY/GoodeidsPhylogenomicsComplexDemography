#!/bin/bash

new_loc="/home/lyusuf/scratch/Biogeography_Goodeids/Mapping/Evolutionary_Rate/meme_extractor"
outdir="/home/lyusuf/scratch/Biogeography_Goodeids/Mapping/Evolutionary_Rate/meme_extractor/sites_under_selection_meme"

#python extractor_meme.py

for all in $new_loc/*;do
fbname=$(basename "$fas" .tsv)
awk '$8<"0.05"' $all > $new_loc/$fbname.temp.txt
awk -F, 'NR==1{print "file_name",$0;next}{print FILENAME , $0}' $fbname.temp.txt > $outdir/$fbname.sites_under_selection.txt
done





