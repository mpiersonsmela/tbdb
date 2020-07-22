#compare_scores.py

import sys
import pandas as pd

def overlap(name1, name2):
    accession1, l1, r1 = extract_from_accession(name1)
    accession2, l2, r2 = extract_from_accession(name2)
    if accession1 != accession2:
        return False
    return (l1 <= r2 and r1 >= l2) 

def extract_from_accession(name):
    accession = name.split(':')[0]
    left_end = int(name.split(':')[1].split('-')[0])
    right_end = int(name.split(':')[1].split('-')[1])
    if left_end > right_end: #need to swap
        temp = left_end
        left_end = right_end
        right_end = temp
    return accession, left_end, right_end


model1 = pd.read_csv(sys.argv[1])
model2 = pd.read_csv(sys.argv[2])
combined = [] #List of rows in output

for i in range(len(model1)):
    name1 = model1["Name"][i]
    include = True
    for j in range(len(model2)):
        name2 = model2["Name"][j]
        if name1 == name2:  #overlap(name1, name2):
            if (model1["Score"][i] < model2["Score"][j]) or (not pd.isna(model1["warnings"][i]) and "TRUNCATED_STEM_1" in model1["warnings"][i] and model2["Score"][j] > 30):
                combined.append(model2.iloc[[j]])
                print("Swapped " + name2)
                include = False
            break
    if include:
        combined.append(model1.iloc[[i]])

output = pd.concat(combined)
print("Combined: " + str(len(output)))
filtered = output.drop_duplicates(subset = "Sequence")
print("Deduplicated: " + str(len(filtered)))
filtered.to_csv(sys.argv[3])
