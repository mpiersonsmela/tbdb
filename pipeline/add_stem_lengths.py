import sys
import pandas as pd


def get_lengths(s1_start, s1_end, antiterm_start, antiterm_end, term_start, term_end, other_stems):
    
    s1_length = None
    
    s2_region_start = None
    s2_region_end = None
    s2_region_length = None
    
    s3_start = None
    s3_end = None
    s3_length = None
    
    antiterm_length = None
    term_length = None
    
    if not pd.isna(s1_start) and not pd.isna(s1_end):
        s1_start = int(s1_start)
        s1_end = int(s1_end)
        if s1_start > 1 and s1_end > 1:
            s1_length = s1_end - s1_start + 1
    
    #note: "other_stems" is actually all stems
    if not pd.isna(other_stems) and not pd.isna(s1_start) and not pd.isna(antiterm_start):
        other_stem_list = other_stems.replace('[', '').replace(']', '').replace(' ','').split(',')
        if len(other_stem_list) > 4: #there are more than 2 stems
            s3_start = int(other_stem_list[-4])
            s3_end = int(other_stem_list[-3])
            s3_length = s3_end - s3_start + 1 #the last stem
            
        if len(other_stem_list) > 6: #there are more than 3 stems
            s2_region_start = int(other_stem_list[2])
            s2_region_end = int(other_stem_list[-5])
            s2_region_length = s2_region_end - s2_region_start + 1
    
    if not pd.isna(antiterm_start) and not pd.isna(antiterm_end):
        antiterm_start = int(antiterm_start)
        antiterm_end = int(antiterm_end)
        if antiterm_start > 1 and antiterm_end > 1:
            antiterm_length = antiterm_end - antiterm_start + 1
    
    if not pd.isna(term_start) and not pd.isna(term_end):
        term_start = int(term_start)
        term_end = int(term_end)
        if term_start > 1 and term_end > 1:
            term_length = term_end - term_start + 1
    
    return s1_length, s2_region_start, s2_region_end, s2_region_length, s3_start, s3_end, s3_length, antiterm_length, term_length

tboxes = pd.read_csv(sys.argv[1])

tboxes[['stem1_length', 'stem2_region_start', 'stem2_region_end', 'stem2_region_length', 'stem3_start', 'stem3_end', 'stem3_length','antiterm_length', 'term_length']] = tboxes.apply(lambda x: get_lengths(x['s1_start'], x['s1_end'], x['antiterm_start'], x['antiterm_end'], x['term_start'], x['term_end'], x['other_stems']), axis = 'columns', result_type = 'expand')

tboxes.to_csv(sys.argv[2], index = False)