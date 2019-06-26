# Library & Loading-----------------------------------------------------------------
list.of.packages <- c("indRa","limma","tidyverse","EBSeq","edgeR","stringr","tibble","readr",'BiocManager','DESeq2',"dplyr",
                      "data.table","tximport","tximportData","pasilla","IHW","tidyr", "fgsea","ggplot2","GSA","stats",
                      "reticulate")
suppressMessages(lapply(list.of.packages, require, character.only = TRUE))

#load custom functions
R_loc <- './R/'
functions <- list.files(R_loc,pattern = "*.R")
suppressMessages(lapply(paste('./R/',functions,sep=""), source)) #load in custom functions


# Read in Data-------------------------------------------------------------
#possible meta & counts data
meta_loc <- './data/DGE_data/'
meta_files <- list.files(meta_loc,pattern = "*meta*.csv")
counts_loc <- './data/DGE_data/'
counts_files <- list.files(counts_loc,pattern = "*counts*.csv")
#GMT files
GMT_loc<-'./data/TF_GSEA_GMT_FILES/'
GMT_files <- list.files(GMT_loc,pattern = "*.gmt")

# Quality Control Counts & Meta -------------------------------------------
#[TODO]


# DGE Run -----------------------------------------------------------------
#DESeq2_run(counts_loc,counts_files[1],meta_loc,meta_files[1],column='Target_JAK2')

# TF GSEA Analysis --------------------------------------------------------
Diff_Gene_Expression_Results_loc<-"./Results/Differential_Gene_Analysis/DESeq2/"
Diff_Gene_Expression_Results<-list.files(Diff_Gene_Expression_Results_loc,pattern = "*.tsv",recursive = FALSE)

#Rank GSEA Analysis for DESeq2 & P Adjusted Analysis
#somtimes an error: 'package:stats' not available when loading, weird Rstudio bug that restarting Rstudio seems to fix
rank_GSEA(Diff_Gene_Expression_Results_loc,Diff_Gene_Expression_Results,GMT_loc,GMT_files,
          output_loc='./Results/Gene_Set_Enrichment_Analysis/',"DESeq2","p_adjusted")

# TF Expression Visualiation ----------------------------------------------
GSEA_results=



Gene_Set_Enrichment_loc<-"./Results/Gene_Set_Enrichment_Analysis/"
Gene_Set_Enrichment_Results<-list.files(Gene_Set_Enrichment_loc,pattern = "*.tsv",recursive = FALSE)

p <- TF_Violin_Plot(Gene_Set_Enrichment_loc,Gene_Set_Enrichment_Results[1],counts_loc,counts_files[1],
                       meta_loc,meta_files[1],1,'Target_JAK2')

# Pathway Analysis --------------------------------------------------------
## If using virtualenv
#reticulate::use_virtualenv(INDRA, required=TRUE) #while possible to make virtual envs are not recommended for Windows Environment
reticulate::use_condaenv('INDRA',required=TRUE)
Gene_Set_Enrichment_loc<-"./Results/Gene_Set_Enrichment_Analysis/"
Gene_Set_Enrichment_Results<-list.files(Gene_Set_Enrichment_loc,pattern = "*.tsv",recursive = FALSE)

## Access to INDRA is available through indra() function
indra() # Module(indra)

#[TODO]:  Needs updating
dijkstra()
PW <- dijkstra( "JAK2", trgts=c("NFKB1", "STAT1", "STAT2", "STAT3", "IRF1", "IRF3") )
P <- with(PW, setNames(Path, Gene))
plotPaths( PW$Path ) + ggthemes::scale_color_few()












