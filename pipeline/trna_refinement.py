#trna_refinement.py
#Code by Jorge Marchand
#Refines codon predictions according to tRNA and downstream gene matching

import pandas as pd
import sys

pd.options.mode.chained_assignment = None

def trna_refinement(infile,outfile):

    predseq=pd.read_csv(infile, low_memory=False)
    lut = pd.read_csv('LUTs/rccodonLUT.csv', low_memory=False)
    codelut = pd.read_csv('LUTs/LogicTreeLUT.csv')

    predseq['refine_codon'] = None
    predseq['refine_codon_top'] = None
    predseq['refine_codon_alt_1'] = None
    predseq['refine_codon_alt_2'] = None
    predseq['refine_codon_num'] = None
    predseq['refine_codon_io'] = None
    
    for i in range(0,len(predseq)):
        try:
            pr=predseq.iloc[i]
            
            #Get reading frames and convert to codons
            rfs = codelut.loc[codelut['Code']==pr['refine_codon_code'],['refine_rf']].values[0][0].split(',')
            
            
            #Prioirtize frames, resort based on order of best to worse
            refine_codon_list=[]
            priority_frames=[1,0,2]
            
            ###PRIORITY FRAMES CAN COME FROM LUT TO RESOLVE EXAMPLES
            
            j=0
            for ii in range(0,len(priority_frames)):
                rfslist=''.join(rfs)
                if  str(priority_frames[ii]) in rfslist:
                    rf=priority_frames[ii]
                    codon=pr['codon_region'][rf:rf+3].replace('T','U')
                    refine_codon_list.append(codon)
                    if j == 0:
                        predseq['refine_codon_top'].iloc[i] = codon
                        predseq['refine_codon_num'].iloc[i] = 1
                        predseq['refine_codon_io'].iloc[i] =predseq['codon'].iloc[i]==codon
                        j=1
                    elif j == 1:
                        predseq['refine_codon_alt_1'].iloc[i] = codon
                        predseq['refine_codon_num'].iloc[i] = 2

                        j=2
                    elif j == 2:
                        predseq['refine_codon_alt_2'].iloc[i] = codon
                        predseq['refine_codon_num'].iloc[i] = 3

                        j=3
            predseq['refine_codon'].iloc[i] = refine_codon_list

            if i%100==0:
                print('Progress: '+str(i))
       
        except:
            print('failed on '+str(i))
            predseq['refine_codon_top'].iloc[i] = predseq['codon'].iloc[i]
            predseq['refine_codon_io'].iloc[i] = True

    predseq.to_csv(outfile, index=False, header=True)

trna_refinement(sys.argv[1],sys.argv[2])