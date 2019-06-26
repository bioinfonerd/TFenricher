#Title: Gene Set Enrichment Analysis in 3 different ways
#Purpose:
#Input:
#Created: Nathan T. Johnson
#Version: 0.0.3
#

#[TODO] Add functionality for dynamic changes
#[TODO] Add ability for different location input from command line
#[TODO] Include:  fGSEA (log fold change ranked, p-value ranked), Hypergeom (padj)
#[TODO] Add Documentation
#[TODO] Update diffent diff gene using new functions
#[TODO] Update combination of methods available within function call (EBSeq, EdgeR, DESeq2 updated, but only DESeq2 tested)
#[TODO] Change GSEA function to work for multiple methods

print("Starting GSEA Analysis: Hypergeometric & fGSEA")

#if running on local windows machine
#setwd('C://Users//Nathan//Documents//@Harvard//Artem//Transcription Factor Analysis//git//Transcription-Factor-Analysis')

# Libraries---------------------------------------------------
print("Loading Libraries")
list.of.packages <- c("fgsea","data.table","ggplot2","GSA","dplyr","stats")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){BiocManager::install(new.packages, update = TRUE)} else {lapply(list.of.packages, require, character.only = TRUE)}

# GMT Files for what Transcription Factor each Binding Partners are--------------
print("Finding GMT Files")
gmt_location<-"./data/TF_GSEA_GMT_FILES/"
gmt <- list.files(gmt_location,pattern = "\\.gmt$")

# Functions ---------------------------------------------------------------

diff_method<-function(method_input,type_input){
  "
  Since different Differential Gene Expression Methods organize output files differently, need condition how to split
  Organization, two lists that organize method and type with what split column should be used
  All possible combinations of both method and type are compiled with a master list then the split column is assigned as needed
  
  Input: Differential Expression Method, Type of Input
   method_input: str (Possible: DESeq2,EdgeR,EBSeq,Limma)
   type_input: str   (Possible: log2fold,p_adjusted,p_value)
  Output: Column Name on what to split
    sep: str 
  "
  method=c("DESeq2","EdgeR","EBSeq","Limma") #list of methods
  type=c("log2fold","p_adjusted","p_value") #list of data types to use
  master_list <- paste(rep(method, each = length(type)), type, sep = ".") #All possible combinations of both method and type
  
  #assign values
  master_list["DESeq2.log2fold"] = "log2FoldChange"
  master_list["DESeq2.p_adjusted"] = "padj" 
  master_list["DESeq2.p_value"] = "pvalue"
  master_list["EdgeR.log2fold"] = "logFC"
  master_list["EBSeq.log2fold"] = "RealFC"
  
  #check if available and assign how to separate 
  if (paste(method_input,type_input,sep=".") %in% names(master_list)){
    sep = master_list[paste(method_input,type_input,sep=".")]
  }  else {print('not available')}
  
  return(sep)
}

fGSEA<-function(files,gmt,output_loc,method_input,type_input){
  "
  Since different Differential Gene Expression Methods organize output files differently, need condition how to split
  Organization, two lists that organize method and type with what split column should be used
  All possible combinations of both method and type are compiled with a master list then the split column is assigned as needed
  
  Input: Differential Expression Method, Type of Input
   method_input: str (Possible: DESeq2,EdgeR,EBSeq,Limma)
   type_input: str   (Possible: log2fold,p_adjusted,p_value)
  Output: Column Name on what to split
    sep: str 
  "
  print(diff_method(method_input,type_input))
  
  for (i in 1:length(files)) {
    print(files[i])
    table<-read.csv(paste(c(DGA_location,files[i]),collapse = ''),sep='\t') #read in data
    #table %>% arrange(padj) #sort by ascending order
    rank_gene_pval <- as.vector(t(table[diff_method(method_input,type_input)])) #fill vector with numeric
    names(rank_gene_pval) <- row.names(table) #id each element  
    rank_gene_pval<-na.omit(rank_gene_pval) #remove NAs
    
    for (n in 1:length(gmt)){
      print(n)
      gsea<-gmtPathways(paste(c(gmt_location,gmt[n]),collapse = ''))
      
      output <- fgsea(pathways = gsea,
                      stats = rank_gene_pval,
                      #minSize = 10, [TODO]: Needs to be explored its effect on ranking
                      nperm=10000) #included all TF targets independent of size
      
      #write data
     fwrite(output, file=paste(c(output_loc,files[i],gmt[n],
                                 names(diff_method(method_input,type_input)),
                                 '.tsv'),collapse = ''), sep="\t", sep2=c("", " ", ""))
      #fwrite(output, file=paste(c(output_loc,files[i],gmt[n],'.tsv'),collapse = ''), sep="\t", sep2=c("", " ", ""))
      
    }
  }
  
} 

hypergeometric<-function(files,gmt,output_loc,method_input,type_input){
  
}

# DESeq2 Analysis ---------------------------------------------------------------
print("Starting DESeq2 Analysis")

#DGE1 data
DGA_location <- "./data/DGE1_DGE_DRUG_NOTDRUG/DESeq2/"
output_loc<-"./results/DGE1_fGSEA_DRUG_NOTDRUG/DESeq2/"
files <- list.files(DGA_location,pattern = "\\.tsv$")

#Run fGSEA on all types
#fGSEA(files,gmt,output_loc) #assume DESeq2 data for now, sort by p-value adjusted, +/- by log2 fold change

fGSEA(files,gmt,output_loc,"DESeq2","log2fold") #Log 2 Fold Analysis
fGSEA(files,gmt,output_loc,"DESeq2","p_adjusted") #P Adjusted Analysis
fGSEA(files,gmt,output_loc,"DESeq2","p_value") #Pval

#DGE2 data
DGA_location <- "./data/DGE2_DGE_DRUG_NOTDRUG/DESeq2/"
output_loc<-"./results/DGE2_fGSEA_DRUG_NOTDRUG/DESeq2/"
files <- list.files(DGA_location,pattern = "\\.tsv$")

#Run fGSEA on all types
#fGSEA(files,gmt,output_loc) #assume DESeq2 data for now, sort by p-value adjusted, +/- by log2 fold change

fGSEA(files,gmt,output_loc,"DESeq2","log2fold") #Log 2 Fold Analysis
fGSEA(files,gmt,output_loc,"DESeq2","p_adjusted") #P Adjusted Analysis
fGSEA(files,gmt,output_loc,"DESeq2","p_value") #Pval

#Hypergeometric 
DGA_location <- "./data/DGE2_DGE_DRUG_NOTDRUG/DESeq2/"
output_loc<-"./results/DGE2_fGSEA_DRUG_NOTDRUG/"
files <- list.files(DGA_location,pattern = "\\.tsv$")

table<-read.csv(paste(c(DGA_location,files[i]),collapse = ''),sep='\t',row.names = NULL) #named vector of ranked gene by p-value 

padj

test <- table %>% filter(padj <= 0.05)


rank_gene_pval <- as.vector(t(table$log2FoldChange)) #fill vector with numeric
names(rank_gene_pval) <- row.names(table) #id each element  
rank_gene_pval<-na.omit(rank_gene_pval) #remove NAs


n=2
gsea<-gmtPathways(paste(c(gmt_location,gmt[n]),collapse = ''))


length(unlist(gsea[3]))

length(unlist(gsea[3]))



#calculating based on probability of getting exactly X 

phyper(possible_func_genes, total_func_gene, total_non_func_gene, genes_within_list)

"
genes_within_list   = length of diff exp genes
possible_func_genes = number of genes associated with particular function
total_func_gene     = total number of genes annotated with function
total_non_func_gene = total number annotated genes within genome - total_func_gene  
"



# edgeR -------------------------------------------------------------------
# output_loc<-'/home/bionerd/Harvard/Artem/Transcription Factor Analysis/DGE2_GSEA_DRUG_DMSO/edgeR/'
# location <- "/home/bionerd/Harvard/Artem/Alzheimers/edgeR/DGE2_Drug_DMSO/"
# files <- list.files(location,pattern = "\\.tsv$")
# 
# for (i in 1:length(files)) {
#   print(files[i])
#   table<-read.csv(paste(c(location,files[i]),collapse = ''),sep='\t') #named vector of ranked gene by p-value 
#   rank_gene_pval <- as.vector(t(table$logFC)) #fill vector with numeric
#   names(rank_gene_pval) <- row.names(table) #id each element  
#   rank_gene_pval<-na.omit(rank_gene_pval) #remove NAs
#   
#   for (n in 1:length(gmt)){
#     print(n)
#     gsea<-gmtPathways(paste(c(gmt_location,gmt[n]),collapse = ''))
#     
#     #run analysis
#     # output <- fgsea(pathways = gsea,
#     #                 stats = rank_gene_pval,
#     #                 minSize=15,
#     #                 maxSize=500,
#     #                 nperm=10000)
#     
#     output <- fgsea(pathways = gsea,
#                     stats = rank_gene_pval,
#                     nperm=10000) #included all TF targets independent of size
#     
#     #write data
#     fwrite(output, file=paste(c(output_loc,files[i],gmt[n],'.tsv'),collapse = ''), sep="\t", sep2=c("", " ", ""))
#   }
# }



# EBSeq PostFC-------------------------------------------------------------------
# output_loc<-'/home/bionerd/Harvard/Artem/Transcription Factor Analysis/DGE2_GSEA_DRUG_DMSO/EBSeq_PostFC/'
# location <- "/home/bionerd/Harvard/Artem/Alzheimers/EBSeq/DGE2_Drug_DMSO/"
# files <- list.files(location,pattern = "\\.tsv$")
# 
# for (i in 1:length(files)) {
#   print(files[i])
#   table<-read.csv(paste(c(location,files[i]),collapse = ''),sep='\t') #named vector of ranked gene by p-value 
#   rank_gene_pval <- as.vector(t(table$PostFC)) #fill vector with numeric
#   names(rank_gene_pval) <- row.names(table) #id each element  
#   rank_gene_pval<-na.omit(rank_gene_pval) #remove NAs
#   
#   for (n in 1:length(gmt)){
#     print(n)
#     gsea<-gmtPathways(paste(c(gmt_location,gmt[n]),collapse = ''))
#     
#     #run analysis
#     # output <- fgsea(pathways = gsea,
#     #                 stats = rank_gene_pval,
#     #                 minSize=15,
#     #                 maxSize=500,
#     #                 nperm=10000)
#     
#     output <- fgsea(pathways = gsea,
#                     stats = rank_gene_pval,
#                     nperm=10000) #included all TF targets independent of size
#     
#     #write data
#     fwrite(output, file=paste(c(output_loc,files[i],gmt[n],'.tsv'),collapse = ''), sep="\t", sep2=c("", " ", ""))
#   }
# }


# EBSeq RealFC-------------------------------------------------------------------
# output_loc<-'/home/bionerd/Harvard/Artem/Transcription Factor Analysis/DGE2_GSEA_DRUG_DMSO/EBSeq_RealFC/'
# location <- "/home/bionerd/Harvard/Artem/Alzheimers/EBSeq/DGE2_Drug_DMSO/"
# files <- list.files(location,pattern = "\\.tsv$")
# 
# for (i in 1:length(files)) {
#   print(files[i])
#   table<-read.csv(paste(c(location,files[i]),collapse = ''),sep='\t') #named vector of ranked gene by p-value 
#   rank_gene_pval <- as.vector(t(table$RealFC)) #fill vector with numeric
#   names(rank_gene_pval) <- row.names(table) #id each element  
#   rank_gene_pval<-na.omit(rank_gene_pval) #remove NAs
#   
#   for (n in 1:length(gmt)){
#     print(n)
#     gsea<-gmtPathways(paste(c(gmt_location,gmt[n]),collapse = ''))
#     
#     #run analysis
#     # output <- fgsea(pathways = gsea,
#     #                 stats = rank_gene_pval,
#     #                 minSize=15,
#     #                 maxSize=500,
#     #                 nperm=10000)
#     
#     output <- fgsea(pathways = gsea,
#                     stats = rank_gene_pval,
#                     nperm=10000) #included all TF targets independent of size
#     
#     #write data
#     fwrite(output, file=paste(c(output_loc,files[i],gmt[n],'.tsv'),collapse = ''), sep="\t", sep2=c("", " ", ""))
#   }
# }