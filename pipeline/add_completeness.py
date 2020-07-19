import pandas as pd
import sys

def add_completeness(codon, a_struct, a_errors, t_struct, t_errors, tRNA):
    complete = ""
    bad_aterm = pd.isna(a_struct) or not(pd.isna(a_errors))
    bad_term = pd.isna(t_struct) or not(pd.isna(t_errors))
    if pd.isna(codon) and bad_aterm and bad_term:
        complete = "None"
    elif not (pd.isna(codon) or bad_aterm or bad_term):
        complete = "Full"
    else:
        complete = "Partial"
    if pd.isna(tRNA):
        tRNA = "False"
    else:
        tRNA = "True"
    return complete, tRNA

tboxes = pd.read_csv(sys.argv[1])

tboxes[["Completeness","tRNA_match"]] = tboxes.apply(lambda x: add_completeness(x['codon'], x['Trimmed_antiterm_struct'], x['vienna_antiterminator_errors'], x['Trimmed_term_struct'], x['new_term_errors'], x['trna_seq_top']), axis = 'columns', result_type = 'expand')

tboxes.to_csv(sys.argv[2], index = False)