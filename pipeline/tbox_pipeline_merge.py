#tbox_pipeline_merge.py
#By Jorge Marchand and Merrick Pierson Smela
#Merges output .csv files into a single one for downstream processing

import pandas as pd
import sys
import os
import glob

folder_to_merge = sys.argv[1]

os.chdir(folder_to_merge)
extension = 'csv'
all_filenames = [i for i in glob.glob('*.{}'.format(extension))]
#combine all files in the list
combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames ], sort = False, ignore_index = True)

print('Combined = '+str(len(combined_csv)))
#Drop duplicates
upan = combined_csv.drop_duplicates(subset='Sequence')

print('Unique dropped = '+str(len(upan)))
    
#export to csv
upan.to_csv(sys.argv[2], index=False, encoding='utf-8-sig')