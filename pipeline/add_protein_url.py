import pandas as pd
import sys

def add_protein_url(proteinid_string):
    if pd.isna(proteinid_string):
        return None, None
    id_short = proteinid_string.split('|')[1]
    proteinid = proteinid_string.split('|')[1].split('.')[0]
    
    if proteinid_string.startswith('gb|'):
        url = "https://www.ncbi.nlm.nih.gov/protein/"+proteinid
        id_short = "GenBank: "+id_short
    elif proteinid_string.startswith('ref|'):
        url = "https://www.ncbi.nlm.nih.gov/protein/"+proteinid
        id_short = "RefSeq: "+id_short
    elif proteinid_string.startswith('emb|'):
        url = "https://www.ebi.ac.uk/ena/data/view/"+proteinid
        id_short = "ENA: "+id_short
    elif proteinid_string.startswith('dbj|'): 
        url = "http://getentry.ddbj.nig.ac.jp/getentry/dad/"+proteinid
        id_short = "DDBJ: "+id_short
    else:
        url = ""
        print(proteinid_string)
    
    return url, id_short

tboxes = pd.read_csv(sys.argv[1])

tboxes[["protein_url", "protein_id_short"]] = tboxes.apply(lambda x: add_protein_url(x['downstream_protein_id']), axis = 'columns', result_type = 'expand')

tboxes.to_csv(sys.argv[2], index = False)