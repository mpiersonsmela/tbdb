# tbdb
This repository hosts scripts for T-box feature prediction and annotation, as well as a database of annotated T-boxes.

## Downloading the database
The master file used for the TBDB website is: `Master_tboxes.csv.gz`

Column descriptions are provided in the [accompanying website](http://tbdb.io/about.html).

## How the predictions are generated
First, the query sequences are passed to INFERNAL, which uses a covariance model (`RF00230.cm`, obtained from [Rfam](http://rfam.xfam.org/family/RF00230/cm)) to predict secondary structure.

Second, the secondary structures are searched for motifs including Stem 1, the specifier loop, and the antiterminator. The codon and discriminator are extracted from the position of these motifs.

Thermodynamic calculations on antiterminator and terminator folds are performed using ViennaRNA. The NCBI accession numbers of the input sequences are used to gather various annotations, including taxonomy and donwstream gene ontology. tRNAscan-SE was used to generate a list of tRNAs for each organism (cached in the `pipeline/trnascan` folder), and this is compared to the predicted codon and discriminator to identify matching tRNAs.

For each T-box, the most likely codon within the specifier loop is chosen based on position relative to the loop end, with additional refinement using tRNA discriminator base and downstream gene ontology (where present). Alternative codon-frames, where found, are also presented.

# Running your own predictions
Code for running your own predictions is contained in the `/pipeline` and `/translational` subdirectories.

## First time setup
*NOTE: Before running the pipeline, you MUST add your own email and NCBI API key to the `tbox_pipeline_postprocess.py`! Otherwise it will not work.*

For more information see: https://www.ncbi.nlm.nih.gov/account/

### Dependencies
Dependencies are listed in `environment.yml`.

A script `init.sh` is provided to install them automatically using conda.

Code was designed and tested on OSX 10.14; other Unix operating systems should also work but have not been tested.

## Pipeline usage
Place input .fa files in a target directory (for example, `/fasta`).
You can also include output .csv files (for example, outputs from translational T-box prediction)

Then run: `./tboxpredict_batch.sh fasta output.csv [optional score cutoff]`

where score is the INFERNAL score cutoff to use (see [INFERNAL manual](http://eddylab.org/infernal/Userguide.pdf) for how score is calculated). If no input is given, the cutoff will default to 15 (which is relatively low).

## Translational T-box predictions
With input.fa containing your sequences, run: `./tbox_translational.sh input.fa [optional score cutoff]`
To generate an INFERNAL output from a genome file, run: `cmsearch --notrunc --notextw translational_ILE.cm output.txt`