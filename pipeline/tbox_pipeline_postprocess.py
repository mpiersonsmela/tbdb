#tbox_pipeline_postprocess.py
#By Jorge Marchand and Merrick Pierson Smela
#Gets feature data from various databases to annotate T-box predictions

import pandas as pd
import sys
import hashlib
import base64

from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna

import urllib.request

#IMPORTANT: you need to put your email and NCBI API key here in order for this to work
#For more information see: https://www.ncbi.nlm.nih.gov/account/
#You can also run without a key but this will be slower.

Entrez.email = #Your email, required for Entrez
Entrez.api_key =  #Your API key

def add_accession(predseq):
    
    #Code by Jorge Marchand
    
    ###########*****************###########
    #IMPORTANT DELIMITER
    accession_delimiter=':'
    ###########*****************###########
    
    #Accession information for URL
    accession_url = [None] * len(predseq)
    
    #Features for genetic context viewer - Genetic locations
    locus_start = [None] * len(predseq)
    locus_end = [None] * len(predseq)
    locus_view_start = [None] * len(predseq)
    locus_view_end = [None] * len(predseq)
    genome_accession_name = [None] * len(predseq)
    length = [None] * len(predseq)
    
    
    for i in range(0, len(predseq)):
        print('Adding accession for index'+str(i))
        entry = predseq['Name'].iloc[i]
        accession_name=predseq['Name'].iloc[i]
        genome_accession=entry[0:(entry.find('.',0,-1))]
        accession_locus=entry[(entry.find(accession_delimiter,0,-1)):len(entry)]

        genome_accession_name[i]=genome_accession
        
        seq_start=accession_locus[1:(accession_locus.find('-',0,-1))]
        seq_end=accession_locus[(accession_locus.find('-',0,-1))+1:len(accession_locus)]

        #Locus information
        tbox_length=len(predseq['FASTA_sequence'].iloc[i])
        length[i]=str(round(tbox_length)) #Save tbox length as a string


        dif = int(seq_end)-int(seq_start)

        if dif>0:
            real_end=str(int(seq_start)+(tbox_length-1)) #Subtract 1 to get full legnth
            locus_start[i]=seq_start
            locus_end[i]=real_end
            locus_view_start[i]=str(int(seq_start)-500)
            locus_view_end[i]=str(int(real_end)+5000)
            
            accession_url[i] = "https://www.ncbi.nlm.nih.gov/nuccore/"+genome_accession+"?report=genbank&from="+locus_start[i]+"&to="+locus_end[i]
 
        if dif<0:
            real_end=str(int(seq_start)-(tbox_length-1)) #Subtract 1 to get full legnth
            locus_start[i]=seq_start
            locus_end[i]=real_end
            locus_view_end[i]=str(int(seq_start)+500)
            locus_view_start[i]=str(int(real_end)-5000)
    
            accession_url[i] = "https://www.ncbi.nlm.nih.gov/nuccore/"+genome_accession+"?report=genbank&from="+locus_end[i]+"&to="+locus_start[i]+"&strand=2"
            
    
        print('Progress adding accesion ids: '+str(i)+' of '+str(len(predseq)))

    predseq['accession_url']=accession_url
    predseq['accession_name']=genome_accession_name
    predseq['locus_start']=locus_start
    predseq['tbox_length']=length
    predseq['locus_end']=locus_end
    predseq['locus_view_start']=locus_view_start
    predseq['locus_view_end']=locus_view_end
    
    return predseq
    
def add_hash(predseq):
    hashid = [None] * len(predseq)
    unique_name = [None] * len(predseq)
    #tbox_url=[None] * len(predseq)
    for i in range(0, len(predseq)):
    
        print('Hashing : '+str(i))
        if isinstance(predseq['Sequence'][i],str)==1:
            hasher = hashlib.sha256(predseq['Sequence'][i].encode('utf-8'))
            hash =  str(base64.urlsafe_b64encode(hasher.digest()))
            hashid[i] = hash[2:10]
            unique_name[i] = hash[2:10].upper()
            unique_name[i] =  unique_name[i].replace('-','')
            unique_name[i] =  unique_name[i].replace('_','')


    predseq['hash_string']=hashid
    predseq['unique_name']=unique_name

    return predseq
    
def fix_name(predseq):
    predseq.Name.str.replace('c','',regex = False)
    predseq.Name.str.replace('/',':',regex = False)
    return predseq

def add_thermocalc(predseq):
    deltadeltaG = [None] * len(predseq)
    
    for i in range(0, len(predseq)):
        print('Thermocalculation for index - '+str(i))
        try:
            term_e=float(predseq['new_term_energy'].iloc[i])
            anti_e=float(predseq['vienna_antiterminator_energy'].iloc[i])
            ddG=float(term_e)-float(anti_e)
            deltadeltaG[i]=str(ddG)

        except ValueError:
            pass
    
    predseq['deltadelta_g']=deltadeltaG
    
    return predseq

def rc(codon):
    scodon=Seq(codon,generic_rna)
    rc_codon=str(scodon.reverse_complement())
    rc_codon=rc_codon.upper()
    return rc_codon

def wobblepair(base1, base2):
    base1 = base1.upper().replace('T','U')
    base2 = base2.upper().replace('T','U')
    if base1 == 'N' or base2 == 'N':
        return True
    if base1 == rc(base2):
        return True
    if base1 == 'G' and base2 == 'U':
        return True
    if base1 == 'U' and base2 == 'G':
        return True
    return False

def clean(dot_structure, seq):
    if pd.isna(dot_structure) or pd.isna(seq):
        return None
    str_clean = ['.']*len(dot_structure)
    
    str_len = len(dot_structure)
    
    for i in range(str_len - 1):
        if dot_structure[i] == '(':
            l_count = 1 #counting the original '(' that we're looking for
            for j in range(i+1, str_len):
                if dot_structure[j] == '(':
                    l_count += 1
                if dot_structure[j] == ')':
                    l_count -= 1
                if l_count == 0: #We found the pair!
                    if wobblepair(seq[i],seq[j]): #But is it good?
                        str_clean[i] = '('
                        str_clean[j] = ')'
                    break
    return "".join(str_clean) #Convert to string from list


def clean_sequences(predseq):
    predseq["codon_region"]=predseq["codon_region"].str.upper()
    predseq["codon"]=predseq["codon"].str.upper()
    predseq["FASTA_sequence"]=predseq["FASTA_sequence"].str.upper()
    predseq["Name"]=predseq["Name"].str.replace('c','',regex=False)
    
    predseq[["Trimmed_antiterm_struct"]] = predseq.apply(lambda x: clean(x['Trimmed_antiterm_struct'], x['Trimmed_sequence']), axis = 'columns', result_type = 'expand')
    predseq[["Trimmed_term_struct"]] = predseq.apply(lambda x: clean(x['Trimmed_term_struct'], x['Trimmed_sequence']), axis = 'columns', result_type = 'expand')
    
    return predseq

def clean_values(predseq):

    predseq.Score.fillna(0)
    predseq.Bias.fillna(0)
    predseq.CM_accuracy.fillna(0)
    predseq.E_value.fillna(0)
    predseq["E_value"] = predseq["E_value"].map("{:.1e}".format)
            
    for i in range(0, len(predseq)):
        try:
            predseq["deltadelta_g"].iloc[i]=str("{:.2e}".format(predseq["deltadelta_g"].iloc[i]))
        except:
            pass
    
    return predseq
    
def add_organism_and_taxid(predseq):
    
    LUTfile = "LUTs/OrganismLUT.csv"
    LUT = pd.read_csv(LUTfile)
    
    if 'TaxId' not in predseq:
        predseq['TaxId'] = None
    if 'GBSeq_organism' not in predseq:
        predseq['GBSeq_organism'] = None
    if 'phylum' not in predseq:
        predseq['phylum'] = None
    if 'class' not in predseq:
        predseq['class'] = None
    if 'order' not in predseq:
        predseq['order'] = None
    if 'family' not in predseq:
        predseq['family'] = None
    if 'genus' not in predseq:
        predseq['genus'] = None
     
    for i in range(0, len(predseq)):
        taxid = ""
        print('Progress adding organism names: '+str(i)+' of '+str(len(predseq)))
        
        if isinstance(predseq['GBSeq_organism'].iloc[i],str)==0:
            genome_accession=(predseq['Name'].iloc[i])[0:((predseq['Name'].iloc[i]).find('.',0,-1))]
            
            try:
                row =  LUT[LUT['Accession'] == genome_accession]
                predseq.loc[pd.Index([i]),['TaxId','GBSeq_organism','phylum','class','order','family','genus']] = row.values[0,1:8]
            
            except IndexError: #not found in LUT
                try:
                    genome = Entrez.efetch(db="nucleotide", id=genome_accession, rettype="docsum", retmode="xml")
                    target_records = Entrez.read(genome)


                    taxid=target_records[0]['TaxId']


                    genome = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")

                    data = Entrez.read(genome)
                
                #Initialize empty 
                    tax_phylum = 'NA'
                    tax_class = 'NA'
                    tax_order  = 'NA'
                    tax_family = 'NA'
                    tax_genus = 'NA'
                    tax_organism = 'NA'

                    tax_organism = data[0]['ScientificName']
                    tax_phylum = {d['Rank']:d['ScientificName'] for d in data[0]['LineageEx'] if d['Rank'] in ['phylum']}.get('phylum')
                    tax_class = {d['Rank']:d['ScientificName'] for d in data[0]['LineageEx'] if d['Rank'] in ['class']}.get('class')
                    tax_order = {d['Rank']:d['ScientificName'] for d in data[0]['LineageEx'] if d['Rank'] in ['order']}.get('order')
                    tax_family = {d['Rank']:d['ScientificName'] for d in data[0]['LineageEx'] if d['Rank'] in ['family']}.get('family')
                    tax_genus = {d['Rank']:d['ScientificName'] for d in data[0]['LineageEx'] if d['Rank'] in ['genus']}.get('genus')
        

                    print('Finding taxonomy for '+genome_accession+' i='+str(i))
                
                
                    try:
                        predseq['TaxId'].iloc[i]=taxid
                    except:
                        predseq['TaxId'].iloc[i]='NA'
                
                    try:
                        predseq['phylum'].iloc[i]=tax_phylum
                    except:
                        predseq['phylum'].iloc[i]='NA'
            
                    try:
                        predseq['class'].iloc[i]=tax_class
                    except:
                        predseq['class'].iloc[i]='NA'

                    try:
                        predseq['order'].iloc[i]=tax_order
                    except:
                        predseq['order'].iloc[i]='NA'

                    try:
                        predseq['family'].iloc[i]=tax_family
                    except:
                        predseq['family'].iloc[i]='NA'

                    try:
                        predseq['genus'].iloc[i]=tax_genus
                    except:
                        predseq['genus'].iloc[i]='NA'

                    try:
                        predseq['GBSeq_organism'].iloc[i]=tax_organism
                    except:
                        predseq['GBSeq_organism'].iloc[i]='NA'
                    
                    #Update the LUT
                    LUT = LUT.append({'Accession':genome_accession,'TaxId':taxid,'GBSeq_organism':tax_organism,'phylum':tax_phylum,'class':tax_class,'order':tax_order,'family':tax_family,'genus':tax_genus}, ignore_index = True)
                
                    if i % 10 == 0:
                        LUT.to_csv(LUTfile, index = False, header = True)
                
                except:
                    print('Missing TAXID information for '+genome_accession)
    
    #Save the LUT
    LUT.to_csv(LUTfile, index = False, header = True)
    return predseq

def add_dsgene(predseq):

    LUTfile = "LUTs/dsgeneLUT.csv"
    dsgeneLUT = pd.read_csv(LUTfile)

    #Initialize
    if 'downstream_protein' not in predseq:
        predseq['downstream_protein'] = None
    if 'downstream_protein_id' not in predseq:
        predseq['downstream_protein_id'] = None
    if 'downstream_protein_EC' not in predseq:
        predseq['downstream_protein_EC'] = None
    
    for i in range(0, len(predseq)):
        if pd.isna(predseq['downstream_protein'].iloc[i]):
            #Get required fields from input file
            genome_accession=(predseq['Name'].iloc[i])[0:((predseq['Name'].iloc[i]).find('.',0,-1))]
            locus_s=int(predseq['locus_start'].iloc[i])
            locus_e=int(predseq['locus_end'].iloc[i])
            
            try:
                row = dsgeneLUT.loc[(dsgeneLUT['accession_name'] == genome_accession) & (dsgeneLUT['locus_start'] == locus_s) & (dsgeneLUT['locus_end'] == locus_e)]
                predseq.loc[pd.Index([i]),['downstream_protein','downstream_protein_id','downstream_protein_EC']] = row.values[0,3:6]
                #print("found LUT " + genome_accession)
            
            except IndexError:#KeyError: #not found in LUT
                if locus_s<locus_e:
                    chr_s = str(locus_e) #Start at end of Tbox
                    chr_e = str(locus_e + 500) #500 bp after end of Tbox
                if locus_s>locus_e:
                    chr_s = str(locus_e - 500)
                    chr_e = str(locus_e)  #end 500bp before, maybe dont care about strand here
            
                handle = Entrez.efetch(db="nuccore", rettype = "ft", id =  genome_accession, seq_start = chr_s, seq_stop = chr_e)
                data = handle.readlines()
            
                #Default entries
                EC_number = 'NA' #Default null entry
                product = 'NA'
                protein_id = 'NA'
            
                for line in data:
                    features = line.strip().split('\t')
                    if features[0] == "product":
                        product = features[1]
                    if features[0] == "protein_id":
                        protein_id = features[1]
                    if features[0] == "EC_number":
                        EC_number = features[1]

          
                predseq['downstream_protein'].iloc[i]=product
                predseq['downstream_protein_id'].iloc[i]=protein_id
                predseq['downstream_protein_EC'].iloc[i]=EC_number
                
                dsgeneLUT = dsgeneLUT.append({'accession_name':genome_accession,'locus_start':locus_s,'locus_end':locus_e,'downstream_protein':product,'downstream_protein_id':protein_id,'downstream_protein_EC':EC_number}, ignore_index = True)
                if i % 10 == 0:
                    dsgeneLUT.to_csv(LUTfile, index = False, header = True)
                       
            print('Progress adding downstream genes: '+str(i)+' of '+str(len(predseq)))

    dsgeneLUT.to_csv(LUTfile,index = False)
    return predseq
    
def add_gene_desc(predseq):
    LUTfile = "LUTs/proteinLUT.csv"
    proteinLUT = pd.read_csv(LUTfile)

    url_conv = 'http://rest.kegg.jp/conv/genes/ncbi-proteinid:'
    url_get = 'http://rest.kegg.jp/get/'
    
    if 'protein_desc' not in predseq:
        predseq['protein_desc'] = None

    for i in range(0, len(predseq)):
        if pd.isna(predseq['protein_desc'].iloc[i]):
            proteinid_string = predseq['downstream_protein_id'].iloc[i]
            #Check if the string is already in the LUT
            try:
                row =  proteinLUT[proteinLUT['downstream_protein_id'] == proteinid_string]
                predseq['protein_desc'].iloc[i] = row.values[0,1]
            except IndexError: #not found in LUT, need to use API
                protein_desc = ""
                try: #use UniProt API first
                    proteinid = proteinid_string.split('|')[1].split('.')[0]
                    try:
                        uniprotURL = "https://www.uniprot.org/uniprot/?query="+proteinid+"&format=tab&columns=protein_names,go(biological_process)"
                        with urllib.request.urlopen(uniprotURL) as f:
                            protein_desc = f.readlines()[1].decode('utf-8').replace('\t',' | ')
                    
                    except Exception as e:
                        print(e)
                        print("Uniprot retrieval failed, trying other databases")
                        if proteinid_string.startswith('gb|'): #use KEGG
                            with urllib.request.urlopen(url_conv + proteinid) as f:
                                KEGG_id = f.read().decode('utf-8').split()[1]

                            with urllib.request.urlopen(url_get + KEGG_id) as f:
                                #parsingPathway = False
                                #pathway = ""
                                orthology = ""
                                module = ""
                                for line in f.readlines():
                                    response = line.decode('utf-8')
                                    if response.startswith('ORTHOLOGY'):
                                        orthology = ' '.join(response.split()[1:])
                                    if response.startswith('MODULE'):
                                        module = ' '.join(response.split()[1:])
                                        break
                                    if response.startswith('BRITE'): break
                                    #if parsingPathway:
                                    #    pathway += ' '+' '.join(response.split())
                                    #if response.startswith('PATHWAY'):
                                    #    pathway += ' '.join(response.split()[1:])
                                    #    parsingPathway = True
                            protein_desc = orthology + module

                        elif proteinid_string.startswith('emb|'): #use the ENA database
                            ENA_url = "https://www.ebi.ac.uk/ena/data/view/"+proteinid+"&display=text"
                            gene_abbr = ""
                            product = ""
                            EC_number = ""
                            with urllib.request.urlopen(ENA_url) as f:
                                for line in f.readlines():
                                    response = line.decode('utf-8')
                                    if "/gene=" in response:
                                        gene_abbr = response.split("\"")[1]
                                    if "/product=" in response:
                                        product = response.split("\"")[1]
                                    #if "/EC_number=" in response:
                                    #    EC_number = response.split("\"")[1]
                            protein_desc = gene_abbr + ' ' + product

                        elif proteinid_string.startswith('dbj|'):
                            proteinid = proteinid_string.split('|')[1]
                            DBJ_url = "http://getentry.ddbj.nig.ac.jp/getentry/dad/"+proteinid#+"/?format=flatfile&limit=1"
                            gene_title = ""
                            with urllib.request.urlopen(DBJ_url) as f:
                                parsingDEF = False
                                for line in f.readlines():
                                    response = line.decode('utf-8')
                                    #print(response)
                                    if response.startswith('ACCESSION'):
                                        break
                                    if parsingDEF:
                                        gene_title += ' '.join(response.split())
                                    if response.startswith('DEFINITION'):
                                        gene_title += ' '.join(response.split()[1:])
                                        parsingDEF = True
                            protein_desc = gene_title
                        
                except Exception as e: #not found in a database
                    print(e)
                    print(proteinid_string)
                    print("not found")
                    protein_desc = ''
                
                #Save after every time (even when not found)
                proteinLUT = proteinLUT.append({'downstream_protein_id':proteinid_string,'protein_desc':protein_desc}, ignore_index = True)
                if i % 10 == 0:
                    proteinLUT.to_csv(LUTfile, index = False, header = True)
                
                print('Progress adding protein descriptions: '+str(i)+' of '+str(len(predseq)))
                print(protein_desc)
                predseq['protein_desc'].iloc[i] = protein_desc.strip()
                

    proteinLUT.to_csv(LUTfile, index = False, header = True)
    return predseq
                  
infile = pd.read_csv(sys.argv[1])

infile = fix_name(infile)

infile = add_hash(infile)

infile = clean_sequences(infile)

infile = clean_values(infile)

infile = add_accession(infile)

infile = add_thermocalc(infile)

infile.to_csv('checkpoint.csv', index = False)

#Get organism and downstream gene

infile = add_organism_and_taxid(infile)

infile.to_csv('checkpoint.csv', index = False)

print("Adding downstream genes")

infile = add_dsgene(infile)

infile.to_csv('checkpoint.csv', index = False)

print("Adding gene descriptions")

infile = add_gene_desc(infile)

#Write output
infile.to_csv(sys.argv[2], index = False)