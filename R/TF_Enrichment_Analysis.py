"""
Purpose: Enrichment Analysis for Transcription Factor Targets
Workflow:
#1) define targets
#2) perform differential gene analysis (using only DESEq2 for now)
#3) perform Hypergeometric Test using TF binding targets
#4) output results
Example: JAK2 targets or not, SIK targets or not
Edited: 2019.05.31
Modified: 0.0.2
By: Nathan T. Johnson

Input:
Output:

"""

#[TODO] Make folder structure automatically
#[TODO] Incorporate master differential gene expression pipeline for multiple methods
#[TODO] Allow for dynamic usage

# <editor-fold desc="Library Loading">
print("Loading Libraries")
import scipy.stats as stats
import pandas as pd
import os
import glob
# </editor-fold>

# <editor-fold desc="Functions">
"""
def geneInterest(): #number of genes

def gene_sets_overlapping_geneInterest(): #number of TFs within gene sets

def #

def pathway_enrichment(total_gene: int,
                       total_func_gene: int,
                       total_nonfunc_gene: int,
                       obs: int) -> float:    
    
    #:param total_gene: 
    #:param total_func_gene: 
    #:param total_nonfunc_gene: 
    #:param obs: 
    #:return: 
   

    #total number of genes with A known function
    #total number of genes with FUNCTION of INTEREST
    #total number of genes WITHOUT function of INTEREST
    #number of observations/genes to select
    result=scipy.stats.hypergeom(total_gene,total_func_gene,total_nonfunc_gene,obs)
    return(result)


def hello_name(name: str) -> str:
        

        #:param name:
        #:return:
        
        return("Hello {name}")

hello_name('Nathan")
"""
# </editor-fold>

#assign which drugs are for different targets
drug_class = pd.read_csv('./data/TAS_Profile/drug_tas.csv', sep=",")
drug_class['name'] = drug_class['name'].str.lower()
#binding strength of drugs can be YES = 1,2,3 or NO = 10
#using all concentrations, alternate is to split
#string match is case sensitive, ensure both df's are the same case as there are a mixture of cases

#DGE1
#load in meta table
meta = pd.read_csv('./data/count_data/DGE1_postqc-meta.csv', sep=",")
meta['Drug'] = meta['Drug'].str.lower()

#filter for drugs that are JAK2 or SIK3 or not
target = drug_class['name'][(drug_class['symbol']=='JAK2') & (drug_class['tas'].isin([1,2,3]))]
target_not = drug_class['name'][(drug_class['symbol']=='JAK2') & (drug_class['tas'].isin([10]))]

#add additional columns to compare on
meta_output=meta
meta_output['Target_JAK2']='NA'
meta_output['Target_SIK3']='NA'

#assign JAK2 labels for target or not
filter=["dsRNAmi", "Lipo","DMSO","Drug_control", "naked dsRNA","dsRNA", "LPS"]# remove all lipo, dsRNA, DMSO combinations
meta_output['Target_JAK2'][(~meta_output['Drug'].str.contains('|'.join(filter))) &
            (meta_output['Drug'].isin(target))] ='JAK2'

meta_output['Target_JAK2'][(~meta_output['Drug'].str.contains('|'.join(filter))) &
            (meta_output['Drug'].isin(target_not))] ='NOT_JAK2'

#assign SIK3 labels for target or not
target = drug_class['name'][(drug_class['symbol']=='SIK3') & (drug_class['tas'].isin([1,2,3]))]
target_not = drug_class['name'][(drug_class['symbol']=='SIK3') & (drug_class['tas'].isin([10]))]
filter=["dsRNAmi", "Lipo","DMSO","Drug_control", "naked dsRNA","dsRNA", "LPS"]# remove all lipo, dsRNA, DMSO combinations

meta_output['Target_SIK3'][(~meta_output['Drug'].str.contains('|'.join(filter))) &
            (meta_output['Drug'].isin(target))] ='SIK3'

meta_output['Target_SIK3'][(~meta_output['Drug'].str.contains('|'.join(filter))) &
            (meta_output['Drug'].isin(target_not))] ='NOT_SIK3'

meta_output.to_csv('./data/count_data/DGE1_modified_postqc-meta.csv',sep=',')

#DGE2
meta = pd.read_csv('./data/count_data/DGE2_postqc-meta.csv', sep=",")
meta['Drug'] = meta['Drug'].str.lower()

#filter for drugs that are JAK2 or SIK3 or not
target = drug_class['name'][(drug_class['symbol']=='JAK2') & (drug_class['tas'].isin([1,2,3]))]
target_not = drug_class['name'][(drug_class['symbol']=='JAK2') & (drug_class['tas'].isin([10]))]

#add additional columns to compare on
meta_output=meta
meta_output['Target_JAK2']='NA'
meta_output['Target_SIK3']='NA'

#assign JAK2 labels for target or not
filter=["dsRNAmi", "Lipo","DMSO","Drug_control", "naked dsRNA","dsRNA", "LPS"]# remove all lipo, dsRNA, DMSO combinations
meta_output['Target_JAK2'][(~meta_output['Drug'].str.contains('|'.join(filter))) &
            (meta_output['Drug'].isin(target))] ='JAK2'

meta_output['Target_JAK2'][(~meta_output['Drug'].str.contains('|'.join(filter))) &
            (meta_output['Drug'].isin(target_not))] ='NOT_JAK2'

#assign SIK3 labels for target or not
target = drug_class['name'][(drug_class['symbol']=='SIK3') & (drug_class['tas'].isin([1,2,3]))]
target_not = drug_class['name'][(drug_class['symbol']=='SIK3') & (drug_class['tas'].isin([10]))]
filter=["dsRNAmi", "Lipo","DMSO","Drug_control", "naked dsRNA","dsRNA", "LPS"]# remove all lipo, dsRNA, DMSO combinations

meta_output['Target_SIK3'][(~meta_output['Drug'].str.contains('|'.join(filter))) &
            (meta_output['Drug'].isin(target))] ='SIK3'

meta_output['Target_SIK3'][(~meta_output['Drug'].str.contains('|'.join(filter))) &
            (meta_output['Drug'].isin(target_not))] ='NOT_SIK3'

meta_output.to_csv('./data/count_data/DGE2_modified_postqc-meta.csv',sep=',')

#DEG

#flag for how groups should be handled
# 1 individual classes
# 2 redesign meta table to fit the samples that do or do not bind
#

# DESeq2
command = "Rscript ./bin/DESeq2-analysis.R ./data/count_data/DGE1_postqc-counts.csv " \
          "./data/count_data/DGE1_modified_postqc-meta.csv Target_JAK2/JAK2 Target_JAK2/NOT_JAK2"
print(command)
os.system(command)

command = "Rscript ./bin/DESeq2-analysis.R ./data/count_data/DGE2_postqc-counts.csv " \
          "./data/count_data/DGE2_modified_postqc-meta.csv Target_JAK2/JAK2 Target_JAK2/NOT_JAK2"
print(command)
os.system(command)

command = "Rscript ./bin/DESeq2-analysis.R ./data/count_data/DGE1_postqc-counts.csv " \
          "./data/count_data/DGE1_modified_postqc-meta.csv Target_SIK3/SIK3 Target_SIK3/NOT_SIK3"
print(command)
os.system(command)

command = "Rscript ./bin/DESeq2-analysis.R ./data/count_data/DGE2_postqc-counts.csv " \
          "./data/count_data/DGE2_modified_postqc-meta.csv Target_SIK3/SIK3 Target_SIK3/NOT_SIK3"
print(command)
os.system(command)

##########
#Analysis#
##########

#Positive Control
files=glob.glob('./results/DGE1_fGSEA_DRUG_NOTDRUG/DESeq2/*JAK2*activation_&_repression*.tsv')
positive_control = pd.read_table('./data/isg.txt',header=None) #stat, IRF, NRKB-RELA should be good hits
#25 ISG == TFs

for i in files:
    df = pd.read_table(i,header='infer')
    #number of TF in positive control list
    #print(len(df[(df['padj'] <= 0.5)]))
    #print(len(df[(df['padj'] <= 0.5) & (df['pathway'].isin(positive_control.iloc[:, 0]))]))
    #\tmp=df.sort_values(['padj']).head(n=50)['pathway']
    #print(tmp[tmp.isin(positive_control.iloc[:, 0])])
    df['pathway'][df['leadingEdge'].isin(positive_control.iloc[:, 0])]

    print(df['pathway'])

    #which TF's are consider significant
    print(df['pathway'][(df['padj'] <= 0.5) & (df['pathway'].isin(positive_control.iloc[:, 0]))])

#checking number of samples labelled with/without the drug
df=pd.read_table('./data/count_data/DGE2_modified_postqc-meta.csv',sep=',')
df[df['Target_JAK2'] == 'JAK2'].__len__()
df[df['Target_JAK2'] == 'NOT_JAK2'].__len__()

#take into account fold change status: analysis ranked by p-val
#profile whether TF targets are up/down
#


#SIK3 Analysis
files=glob.glob('./results/DGE1_fGSEA_DRUG_NOTDRUG/DESeq2/*SIK3*activation_&_repression*.tsv')
DGE1= pd.read_table(files[1],header='infer')
files=glob.glob('./results/DGE2_fGSEA_DRUG_NOTDRUG/DESeq2/*SIK3*activation_&_repression*.tsv')
DGE2= pd.read_table(files[1],header='infer')

DGE2.sort_values(['padj']).head(n=50)['pathway'].isin(DGE1.sort_values(['padj']).head(n=50)['pathway'])


    DGE1.sort_values(['padj']).head(n=50)['pathway']





    #number of TF in positive control list
    print(len(DGE1[(DGE1['padj'] <= 0.5)]))
    print(len(DGE2[(DGE2['padj'] <= 0.5)]))

    #print(len(df[(df['padj'] <= 0.5) & (df['pathway'].isin(positive_control.iloc[:, 0]))]))
    #tmp=df.sort_values(['padj']).head(n=50)['pathway']
    #print(tmp[tmp.isin(positive_control.iloc[:, 0])])

    #which TF's are consider significant
    print(df['pathway'][(df['padj'] <= 0.5) & (df['pathway'].isin(positive_control.iloc[:, 0]))])





#save list of TFs
df = pd.read_table(files[1],header='infer')
output=pd.Series(df['pathway'][df['padj'] <= 0.5])
output.to_csv('SIK3_Targeted_TFs.txt',index=False)

#Does SIK3 Modulate ISGs?
for i in files:
    df = pd.read_table(i,header='infer')
    #number of TF in positive control list
    print(len(df[(df['padj'] <= 0.5) & (df['pathway'].isin(positive_control.iloc[:, 0]))]))
    #which TF's are consider significant
    print(df['pathway'][(df['padj'] <= 0.5) & (df['pathway'].isin(positive_control.iloc[:, 0]))])

#Do SIK3 Modulating TF's overlap with Indra's interactions?
df = pd.read_table('./data/IndraBot_SIK3_Results.tsv',sep='\t',header=None)
SIK3_indra_interactions=[]
for i in df.iloc[:,0]: #extract all nouns (genes for the most part?) from column
    if len(i.split('(')) == 3: #handle if none in INDRA output
        SIK3_indra_interactions.append(i.split('(')[1].split(',')[0])
        SIK3_indra_interactions.append(i.split('(')[1].split(',')[1].split(' ')[1])
    if len(i.split('(')) == 4: #handle if NOT none in INDRA output
        SIK3_indra_interactions.append(i.split('(')[1])
        SIK3_indra_interactions.append(i.split('(')[2].split(',')[1].split(' ')[1])

#many non-gene related terms, but keeping as genes shouldn't hit as it has to be an exact match
SIK3_indra_interactions = list(dict.fromkeys(SIK3_indra_interactions)) #creating dictionary will remove duplicates (fastest? don't care as its only ~300 examples)

#how many significant results are in INDRA results
for i in files:
    df = pd.read_table(i,header='infer')
    #number of TF in positive control list
    print(len(df[(df['padj'] <= 0.5) & (df['pathway'].isin(SIK3_indra_interactions))]))
    #which TF's are consider significant
    print(df['pathway'][(df['padj'] <= 0.5) & (df['pathway'].isin(SIK3_indra_interactions))])




#Artem doesn't like filtering, so pause implementation unless results aren't expected
"""
#Hypergeometric

#load in results and select for significant results
df = pd.read_csv('./data/DGE1_DGE_DRUG_DMSO/DESeq2/Drug_a443654,Concentration_0.3_Drug_Drug_control_DESeq2results.tsv',
                 sep="\t")

#targets=df[df['padj']<=0.05].index #selecting for significant DGE

pathway_list = pd.read_csv('data/TF_GSEA_GMT_FILES/activation_&_repression.gmt',sep='\t')

import sys
import scipy
import scipy.stats as stats

scipy.stats.cdf(13588,611,13588-611,59)


total_gene, total_func_gene, total_nonfunc_gene, obs


stats.hypergeom.cdf(

    # total number of genes with A known function
    # total number of genes with FUNCTION of INTEREST
    # total number of genes WITHOUT function of INTEREST
    # number of observations/genes to select
    result=scipy.stats.hypergeom(total_gene, total_func_gene, total_nonfunc_gene, obs)
return (result)



#filter results to significant results

#find

#fGSEA
"""





