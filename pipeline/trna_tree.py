#trna_tree.py
#Code by Jorge Marchand and Merrick Pierson Smela
#Matches T-boxes with tRNAs

import pandas as pd
import sys
import os
import itertools
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna

pd.options.mode.chained_assignment = None

#Initialize LUT
lut = pd.read_csv('LUTs/rccodonLUT.csv', low_memory=False)

####GET REVERSE COMPLEMENT #####
def rc(codon):
    scodon=Seq(codon,generic_rna)
    rc_codon=str(scodon.reverse_complement())
    rc_codon=rc_codon.upper()
    return rc_codon

def wobblepair(base1, base2):
    base1 = base1.upper().replace('T','U')
    base2 = base2.upper().replace('T','U')
    if base1 == rc(base2):
        return True
    if base1 == 'G' and base2 == 'U':
        return True
    if base2 == 'G' and base1 == 'U':
        return True
    return False

def get_discrim(tRNA_seq, tRNA_struct):
    if len(tRNA_seq) != len(tRNA_struct):
        print("Error: structure and sequence length mismatch")
        return None
        
    if tRNA_struct[-5:] == "....." and tRNA_seq[-3:] == 'CCA': 
        return tRNA_seq[-4] #for example, UACCA for Methionine
        
    if tRNA_struct[-3:] == ").." or tRNA_struct[-4:] == ")...":
        return tRNA_seq[-1] #these would get CCA tails added
    
    last_pair = tRNA_struct.rfind(').') #the general case
    if last_pair > len(tRNA_seq) or last_pair == -1:
        print("Error: Discriminator not found")
        return None
    return tRNA_seq[last_pair+1]
    
#####CHECK DS GENE #####
def check_dsgene(pr,codon):
    aa=lut.loc[lut['AC']==rc(codon),['AA']].values[0][0] #Abbreviation
    al=lut.loc[lut['AC']==rc(codon),['ALT']].values[0][0] #Full name

    c1=pr['protein_desc']
    c6=pr['downstream_protein']

    try:
        c12345=c1+c6
    except:
        c12345=c6 #12345 is empty, so only use 6?
    
    try:
        if aa.upper() in c12345.upper() or al.upper() in c12345.upper():
        
            #Check only Alanine/Phenylalanine
            if aa.upper()=='ALA':
                return 'ALA' in c12345.upper().replace('PHENYLALA','.')

            #Check only Leucine/Isoleucine
            elif aa.upper()=='LEU':
                return 'LEU' in c12345.upper().replace('ISOLEU','.')
                    
            elif aa.upper()=='PRO':
                return 'PRO' in c12345.upper().replace('PROTEIN','.')

            else:
                return True
        else:
            return False
    except:
        return False

####CHECK DISCRIMINATOR BASE#####
def check_disc_lut(pr,codon):
    db=pr['discriminator'][3] #Predicted Discriminator Base
    dbc=lut.loc[lut['AC']==rc(codon),['DB']].values[0][0]
    return db in dbc

######GET tRNA #######
def check_disc_trna(genome_accession, codon, pr):
    trna_file='trnascan/'+genome_accession+'_struc.txt'
    db=pr['discriminator'][3] #Predicted Discriminator Base

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
    
    discrim = get_discrim(temp_seq, temp_str)
    print('tRNA Discriminator base = '+discrim)
    print('T-box Anti-Discriminator base = '+db)
    return wobblepair(db, discrim)
    
def check_trna(genome_accession, codon):
    trna_file='trnascan/'+genome_accession+'_struc.txt'
    anticodon=rc(codon).replace('U','T')
    
    g1='grep -A 3 -B 1 "Anticodon: '+anticodon+'" trnascan/'+genome_accession+'_struc.txt > tempfiles/temp.txt'
    g2='grep -A 1 -m 1 "Seq" tempfiles/temp.txt > tempfiles/temp_seq.txt'
    
    os.system(g1)
    os.system(g2)

    return os.stat('tempfiles/temp_seq.txt').st_size > 0

def trna_tree(infile, outfile):
    predseq=pd.read_csv(infile, low_memory=False)
    
    predseq['refine_codon'] = None
    predseq['refine_codon_io'] = None
    predseq['refine_codon_code'] = None

    l=[False,True]
    logic_table = list(itertools.product(l,repeat=6))
    logic_sum= [0]*64

    for i in range(0,len(predseq)):
        try:
            pr=predseq.iloc[i]
            genome_accession=pr['Name'].split('.')[0]
            
            tf_disc=[0,0,0]
            tf_gene=[0,0,0]
            tf_gesc=[0,0,0]
            tf_res=[0,0,0,0,0,0]
            
            print('Currently on: '+str(i))
            
            for rf in range(0,3):
                print("Reading frame: " + str(rf)) 
                codon=pr['codon_region'][rf:rf+3].replace('T','U')
                
                ######################################
                #####Check downstream gene
                dba=check_dsgene(predseq.iloc[i],codon)
                ######################################

                #####Check discriminator base 
                if codon == "CAC" or codon == "CAT":
                    dbb = True #Histidine can have any discriminator
                else:
                    trna_scan_file = 'trnascan/'+genome_accession+'_struc.txt'
                    if os.path.isfile(trna_scan_file):
                        if check_trna(genome_accession, codon):
                            try:
                                dbb=check_disc_trna(genome_accession,codon,pr)
                                #print("tRNA check successful for codon " + codon)
                            except Exception as e:
                                print(e)
                                dbb=check_disc_lut(pr,codon)
                                #print("discriminator error; LUT check successful for codon " + codon)
                        else:
                            dbb=check_disc_lut(pr,codon)
                            #print("tRNA not present; LUT check successful for codon " + codon)
                    else:
                        dbb=check_disc_lut(pr,codon)
                        #print("tRNA file not present; LUT check successful for codon " + codon)
                
                tf_res[rf]=dbb #discriminator
                tf_res[rf+3]=dba #Gene
            
            for y in range(0,len(logic_table)):
                if list(logic_table[y])==tf_res:
                    predseq['refine_codon_code'].iloc[i]=y
                    logic_sum[y]=logic_sum[y]+1
        except Exception as e:
            print(e)
            print('Skipped '+str(i))

    predseq.to_csv(outfile, index=False, header=True)

trna_tree(sys.argv[1],sys.argv[2])