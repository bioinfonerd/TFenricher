TF_Violin_Plot<-function(GSEA_result_loc,GSEA_result,counts_loc,counts_files,
                         meta_loc,meta_files,TF_index,meta_column){
df <- read.csv(paste(GSEA_result_loc,GSEA_result,sep=""),sep='\t',colClasses = "character")
counts<-counts_load(counts_loc,counts_files)
meta<-meta_load(meta_loc,meta_files)
TF=df[TF_index,]$pathway #TF
TARGETS=unlist(strsplit(df[TF_index,]$leadingEdge,split=" "))#TF targets
counts<-counts[TARGETS,] #filter to targets
meta <- meta %>% filter(meta[[meta_column]] != 'NA')
counts <- counts[,meta$Well] #filter to targets
samples <- colnames(counts)
counts <- data.table(t(counts))
counts$Well <- samples
meta <- meta[c("Well",meta_column)]
final <- merge(counts,meta,by="Well")
final <- subset(final, select = -c(Well) )
name<-names(final)
final <- final %>% gather(Gene,Expression,name[1:length(name)-1])
p<-ggplot(final, aes_string(x="Gene", y="Expression", color=meta_column)) +
  geom_violin(trim=FALSE)
return(p)
}
