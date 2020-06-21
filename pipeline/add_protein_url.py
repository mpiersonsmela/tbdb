import pandas as pd
import sys

def add_protein_url(predseq):
    if 'protein_url' not in predseq:
        predseq['protein_url'] = None
        
    for i in range(0, len(predseq)):
        proteinid_string = predseq['downstream_protein_id'].iloc[i]
        if pd.isna(proteinid_string):
            continue
        proteinid = proteinid_string.split('|')[1].split('.')[0]
        if proteinid_string.startswith('gb|'):
            predseq['protein_url'].iloc[i] = "https://www.ncbi.nlm.nih.gov/protein/"+proteinid
        elif proteinid_string.startswith('emb|'):
            predseq['protein_url'].iloc[i] = "https://www.ebi.ac.uk/ena/data/view/"+proteinid
        elif proteinid_string.startswith('dbj|'): 
            predseq['protein_url'].iloc[i] = "http://getentry.ddbj.nig.ac.jp/getentry/dad/"+proteinid
    
    return predseq

tboxes = add_protein_url(pd.read_csv(sys.argv[1]))
tboxes.to_csv(sys.argv[2], index = False)
