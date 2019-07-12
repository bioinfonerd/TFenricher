#DGE analysis across toxicity and concentration for JAK2 & SIK2
#
list.of.packages <- c("indRa","limma","tidyverse","EBSeq","edgeR","stringr","tibble","readr",'BiocManager','DESeq2',"dplyr",
		                            "data.table","tximport","tximportData","pasilla","IHW","tidyr", "fgsea","ggplot2","GSA","stats",
					                          "reticulate","synapser")
suppressMessages(lapply(list.of.packages, require, character.only = TRUE))

#until build reflects functions loaded
devtools::load_all()

#run analysis
DESeq2_run_all('JAK2')
DESeq2_run_all('SIK3')
