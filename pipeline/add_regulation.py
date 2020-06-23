import pandas as pd
import sys

def add_regulation(type, term_errors):
    if pd.isna(type):
        return "Unknown"
    if type == "Transcriptional" and not pd.isna(term_errors):
        return "Unknown"
    return type

tboxes = pd.read_csv(sys.argv[1])

tboxes["Regulation"] = tboxes.apply(lambda x: add_regulation(x['type'], x['new_term_errors']), axis = 'columns', result_type = 'expand')

tboxes.to_csv(sys.argv[2], index = False)