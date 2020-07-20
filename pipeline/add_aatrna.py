#add_aatrna.py
#By Jorge Marchand and Merrick Pierson Smela
#Adds amino acid and tRNA data to annotate T-boxes
#Requires LUT folder
#Requires tRNA scan folder
#Requires tempfiles folder

import pandas as pd
import sys
import os
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna

def rc(codon):
    scodon=Seq(codon,generic_rna)
    rc_codon=str(scodon.reverse_complement())
    return rc_codon
    
def add_CCA(seq, str):
    seq = seq.strip()
    str = str.strip()
    if str[-5:] == "))))." or str[-5:] == ")))..":
        str = str + '...'
        seq = seq + 'CCA'
    return seq, str

def get_trna(genome_accession, codon):
    trna_file='trnascan/'+genome_accession+'_struc.txt'
    anticodons = "LUTs/rccodonLUT.csv"
    antic=pd.read_csv(anticodons)
    anticodon=rc(codon).replace('U','T')
    aa = antic.loc[antic['AC']==rc(codon),['AA']].values[0][0]
    aa= aa.lower().capitalize()

    g1='grep -A 3 -B 1 "'+aa+'\tAnticodon: '+anticodon+'" trnascan/'+genome_accession+'_struc.txt > tempfiles/temp.txt'
    g2='grep -A 1 -m 1 "Seq" tempfiles/temp.txt > tempfiles/temp_seq.txt'
    g3='grep -A 1 -m 1 "Str" tempfiles/temp.txt > tempfiles/temp_str.txt'
    
    os.system(g1)
    os.system(g2)
    os.system(g3)

    temp_seq=open('tempfiles/temp_seq.txt','r').read().split('Seq: ')[1]
    temp_str=open('tempfiles/temp_str.txt','r').read().split('Str: ')[1].replace('>','(').replace('<',')')
    
    return add_CCA(temp_seq, temp_str)
    
def check_trna_file(genome_accession):
    trna_file='trnascan/'+genome_accession+'_struc.txt'
    return os.path.isfile(trna_file)

def check_trna(genome_accession, codon):
    trna_file='trnascan/'+genome_accession+'_struc.txt'
    anticodon=rc(codon).replace('U','T')

    g1='grep -A 3 -B 1 "Anticodon: '+anticodon+'" trnascan/'+genome_accession+'_struc.txt > tempfiles/temp.txt'
    g2='grep -A 1 -m 1 "Seq" tempfiles/temp.txt > tempfiles/temp_seq.txt'

    os.system(g1)
    os.system(g2)
    
    return (os.stat('tempfiles/temp_seq.txt').st_size > 0)

def add_aatrna(infile,outfile):
    
    anticodons = "LUTs/rccodonLUT.csv"
    antic=pd.read_csv(anticodons)
    predseq=pd.read_csv(infile)

    aa=antic['AA'] #amino acid
    ac=antic['AC']#anti codon
    max_n=len(predseq)

    predseq['amino_acid_top']= None
    predseq['trna_family_top']= None
    predseq['trna_seq_top']= None
    predseq['trna_struc_top']= None

    predseq['amino_acid_alt_1']= None
    predseq['trna_family_alt_1']= None
    predseq['trna_seq_alt_1']= None
    predseq['trna_struc_alt_1']= None

    predseq['amino_acid_alt_2']= None
    predseq['trna_family_alt_2']= None
    predseq['trna_seq_alt_2']= None
    predseq['trna_struc_alt_2']= None


    for i in range(0,max_n):
        trna_seq = ""
        trna_struct = ""

        genome_accession =  predseq['Name'].iloc[i].split('.')[0]
        codon_1 =  predseq['refine_codon_top'].iloc[i]
        codon_2 =  predseq['refine_codon_alt_1'].iloc[i]
        codon_3 =  predseq['refine_codon_alt_2'].iloc[i]

        ####CHECK CODON_1
        if isinstance(codon_1,str) and all(c in 'AUGC' for c in codon_1):
            aa = antic.loc[antic['AC']==rc(codon_1),['AA']].values[0][0]
            trna_family = aa+' ('+rc(codon_1)+')'
            predseq['amino_acid_top'].iloc[i]= aa
            predseq['trna_family_top'].iloc[i]= trna_family

            if check_trna_file(genome_accession):
                if check_trna(genome_accession, codon_1):
                    try:
                        trna_seq, trna_struc = get_trna(genome_accession, codon_1)
                        if isinstance(trna_seq,str):
                            predseq['trna_seq_top'].iloc[i]= trna_seq
                            predseq['trna_struc_top'].iloc[i]= trna_struc
                    except Exception as e:
                        print(e)
                else: 
                    print("tRNA Not found1")
            else: 
                print("File Not found1")

        ####CHECK CODON_2
        if isinstance(codon_2,str) and all(c in 'AUGC' for c in codon_2):
            aa = antic.loc[antic['AC']==rc(codon_2),['AA']].values[0][0]
            trna_family = aa+' ('+rc(codon_2)+')'
            predseq['amino_acid_alt_1'].iloc[i]= aa
            predseq['trna_family_alt_1'].iloc[i]= trna_family

            if check_trna_file(genome_accession):
                if check_trna(genome_accession, codon_2):
                    try:
                        trna_seq, trna_struc = get_trna(genome_accession, codon_2)
                        if isinstance(trna_seq,str):
                            predseq['trna_seq_alt_1'].iloc[i]= trna_seq
                            predseq['trna_struc_alt_1'].iloc[i]= trna_struc
                    except Exception as e:
                        print(e)
                else: 
                    print("tRNA Not found2")
            else: 
                print("File Not found2")

        ####CHECK CODON_3
        if isinstance(codon_3,str) and all(c in 'AUGC' for c in codon_3):
            aa = antic.loc[antic['AC']==rc(codon_2),['AA']].values[0][0]
            trna_family = aa+' ('+rc(codon_3)+')'
            predseq['amino_acid_alt_2'].iloc[i]= aa
            predseq['trna_family_alt_2'].iloc[i]= trna_family

            if check_trna_file(genome_accession):
                if check_trna(genome_accession, codon_3):
                    try:
                        trna_seq, trna_struc = get_trna(genome_accession, codon_3)
                        if isinstance(trna_seq,str):
                            predseq['trna_seq_alt_2'].iloc[i]= trna_seq
                            predseq['trna_struc_alt_2'].iloc[i]= trna_struc
                    except Exception as e:
                        print(e)
                else: 
                    print("tRNA Not found3")
            else: 
                print("File Not found3")

    predseq.to_csv(outfile, index=False, header=True)

infile_path = sys.argv[1]
outfile_path = sys.argv[2]
add_aatrna(infile_path, outfile_path)