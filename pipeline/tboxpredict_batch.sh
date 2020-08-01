#!/bin/bash
#Script to run T-box prediction and feature annotation
#Arguments:
#$1: The target directory. Contains .fa inputs and/or preprocessed .csv files
#$2: The output file name
#$3: (Optional: the INFERNAL score cutoff to use)

cm='RF00230.cm' #The Rfam transcriptional T-box covariance model
target="$1" #the target directory
output=$(echo "${2}" | cut -f 1 -d '.') #the output .csv file name

for i in ${target}/*.fa; do
    name=$(echo "${i}" | cut -f 1 -d '.')
    echo "Running cmsearch: ${i}"
    cmsearch --notrunc --notextw $cm ${i} > ${name}_INFERNAL.txt
    echo "Running pipeline: ${i}"
    python3 tbox_pipeline_master.py ${name}_INFERNAL.txt ${name}_PREDICTED.csv ${i} $3 > ${name}_TBOX_LOG.txt
done
#Merge all .csv files in the target directory into one master file
echo "Merging to: ${target}/${output}"
python3 tbox_pipeline_merge.py ${target} ${output}_temp.csv
echo "Postprocessing to: ${output}_temp"
python3 tbox_pipeline_postprocess.py ${target}/${output}_temp.csv tempfiles/${output}_0.csv
#Cleanup intermediate
rm ${target}/${output}_temp.csv
echo "Adding tRNAs"
python3 trna_tree.py tempfiles/${output}_0.csv tempfiles/${output}_1.csv
echo "Refining tRNA predictions"
python3 trna_refinement.py tempfiles/${output}_1.csv tempfiles/${output}_2.csv
echo "Assigning codons"
python3 add_aatrna.py tempfiles/${output}_2.csv tempfiles/${output}_3.csv
echo "Filtering output"
python3 tbox_pipeline_filter.py tempfiles/${output}_3.csv tempfiles/${output}_4.csv
echo "Writing URLs for proteins"
python3 add_protein_url.py tempfiles/${output}_4.csv tempfiles/${output}_5.csv
echo "Adding information on regulation and completeness"
python3 add_regulation.py tempfiles/${output}_5.csv tempfiles/${output}_6.csv
python3 add_completeness.py tempfiles/${output}_6.csv tempfiles/${output}_7.csv
python3 add_stem_lengths.py tempfiles/${output}_7.csv ${output}.csv

echo "Cleaning up"
rm tempfiles/*
rm checkpoint.csv
echo "Done. Predictions saved to ${output}.csv"