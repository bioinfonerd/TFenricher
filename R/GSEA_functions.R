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

rank_GSEA<-function(DGA_location,files,gmt_location,gmt,output_loc,method_input,type_input){
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
      print(gmt(n))
      gsea<-gmtPathways(paste(c(gmt_location,gmt[n]),collapse = ''))
      output <- fgsea(pathways = gsea,
                      stats = rank_gene_pval,
                      minSize = 5,
                      nperm=10000) #included all TF targets independent of size
      #write data
      fwrite(output, file=paste(c(output_loc,files[i],gmt[n],
                                  names(diff_method(method_input,type_input)),
                                  '.tsv'),collapse = ''), sep="\t", sep2=c("", " ", ""))
      #fwrite(output, file=paste(c(output_loc,files[i],gmt[n],'.tsv'),collapse = ''), sep="\t", sep2=c("", " ", ""))

    }
  }

}
