#!/bin/bash
#Runs feature predictions on Class II T-boxes
cm='translational_ILE.cm'
name=$(echo "$1" | cut -f 1 -d '.')

cmsearch --notrunc --notextw $cm $1 > ${name}_INFERNAL.txt
python3 tbox_translational.py ${name}_INFERNAL.txt ${name}_PREDICTED.csv $1 $2 > ${name}_TBOX_LOG.txt
echo "Translational predictions complete"
echo "The output: ${name}_PREDICTED.csv"
echo "can now be input into the pipeline."