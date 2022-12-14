#!/usr/bin/bash

# copying folder structure over to UniqueLigands directory
mkdir ~/research/Annulment/UniqueLigands;
cd ExtractedLigand;
find . -type d > ../dirs.txt;
cd ../UniqueLigands;
xargs mkdir -p < ../dirs.txt;

# For each subfolder in extracted ligand directory
for dir in ~/research/Annulment/ExtractedLigand/*/*; do 
    # Run openbabel duplicate checking 
    obabel $dir/*.mol2 -O $dir/U*.mol -m --unique | /dev/null 
    # Delete empty files generated by OBabel
    find $dir -type f -empty -print -delete;
    # move unique files over to corresponding directory in UniqueLigands folder
    file_path=$dir;
    # prefix="/home/brendan/research/Annulment/ExtractedLigand/";
    suffix="${file_path##*/ExtractedLigand/}";
    mv $dir/U*mol /home/brendan/research/Annulment/UniqueLigands/$suffix;
done;







