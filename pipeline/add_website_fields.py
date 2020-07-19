#add_website_fields.py
import pandas as pd
import sys

def add_website_fields(entry, unique_name, accession_url):
    
    if not pd.isna(unique_name):
        tbox_url = '<a href=tboxes/'+unique_name+'.html >'+unique_name+'</a>'
    else:
        tbox_url = None
    
    accession_delimiter=':'

    accession_url_html = ""
    direction_label = ""
    
    genome_accession=entry[0:(entry.find('.',0,-1))]
    accession_url_html = '<a href='+accession_url+'>'+genome_accession+'</a>'
    
    accession_locus=entry[(entry.find(accession_delimiter,0,-1)):len(entry)]

    seq_start=accession_locus[1:(accession_locus.find('-',0,-1))]
    seq_end=accession_locus[(accession_locus.find('-',0,-1))+1:len(accession_locus)]

    if int(seq_end)-int(seq_start) > 0:
        direction_label='T%96Box%3E%3E%3E'
    else:
        direction_label='%3C%3C%3CT%96Box'
        
    return tbox_url, accession_url_html, direction_label

tboxes = pd.read_csv(sys.argv[1])

tboxes[["tbox_url","accession_url_html","direction_label"]] = tboxes.apply(lambda x: add_website_fields(x['Name'], x['unique_name'], x['accession_url']), axis = 'columns', result_type = 'expand')

tboxes.to_csv(sys.argv[2], index = False)