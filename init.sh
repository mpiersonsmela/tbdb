#!/bin/bash
#Initialize T-box predictions pipeline

echo "Unzipping files"
cd pipeline
tar xzf fasta.tar.gz
rm fasta.tar.gz
tar xzf LUTs.tar.gz
rm LUTs.tar.gz
tar xzf tRNAscan.tar.gz
rm tRNAscan.tar.gz
cd ..
echo "Installing dependencies"
conda init bash
conda env create -f environment.yml
conda activate tbdb

echo "Done. Remember to add your email and NCBI API key to tbox_pipeline_postprocess.py"