TF_enriched_targets_Violin_Plot<-function(DGA_location='../Results/Differential_Gene_Analysis/',
                                          GSEA_result_loc='../Results/Gene_Set_Enrichment_Analysis/EdgeR/',
                                          GSEA_result,
                                          count_file_name,
                                          meta_file_name,
                                          TF_name,
                                          DGA_method,
                                          DGA_type='p_adjusted',
                                          DGA_file){

  #load in data
  counts<-counts_load(count_file_name)
  meta<-meta_load(meta_file_name)
  TAS<-TAS_load()
  Drug_Target<-unlist(str_split(GSEA_result,pattern="_"))[1]
  meta<-TAS_drug_binder(TAS_tibble=TAS,drug_target=Drug_Target,meta_tibble=meta) #filter for samples that are binders or nonbinders for drug target
  df <- read_tsv(paste(GSEA_result_loc,GSEA_result,'.tsv',sep=""))
  TF= df %>% filter(pathway==TF_name) %>% select(pathway) #TF
  TARGETS=unlist(strsplit(unlist(c(df %>% filter(pathway==TF_name) %>% select(leadingEdge))),split = " "))#TF targets
  counts <- counts %>% filter(HUGO %in% TARGETS) #filter to targets
  genes<-counts$HUGO
  counts<-counts %>% select(meta$Well)#select samples
  samples <- colnames(counts)
  counts <- data.table(t(counts))
  colnames(counts)<-genes
  counts$Well <- samples
  meta <- meta[c("Well","binding")]
  final <- merge(counts,meta,by="Well")
  final <- subset(final, select = -c(Well) )
  name<-names(final)
  final <- final %>% gather(Gene,Expression,name[1:length(name)-1])
  p<-ggplot(final, aes_string(x="Gene", y="Expression", color="binding")) +
    geom_violin(trim=FALSE) +
    ggtitle(TF_name) +
    xlab('TF Targets Indicating Enrichment') +
    theme(axis.text.x = element_text(angle = 90))

  #DGE Results For Targets
  print(diff_method(DGA_method,DGA_type))
  tmp<-unlist(str_split(DGA_file,pattern='_'))
  tmp[length(tmp)]
  table<-read.csv(paste(DGA_location,tmp[length(tmp)],'/',DGA_file,'.tsv',sep=""),sep='\t') #read in data
  print('Target Genes Significance')
  if(DGA_method=='EdgeR'){
    print(table %>% filter(genes %in% TARGETS) %>% select (genes,FDR))
    return(p)
  }
  if(DGA_method=='DESeq2'){
    print(table %>% filter(genes %in% TARGETS) %>% select (genes,padj))
    return(p)
  }
  print('DGA Method Not Available')
}


TF_all_targets_heatmap_Plot<-function(GSEA_result_loc='../Results/Gene_Set_Enrichment_Analysis/EdgeR/',
                                                                       GSEA_result,
                                                                       count_file_name,
                                                                       meta_file_name,
                                                                       TF_name,
                                                                       meta_column){


  #heatmap of binder vs non binder

  #heatmap of binder vs non binder as all individual drugs



   #load in data
  counts<-counts_load(count_file_name)
  meta<-meta_load(meta_file_name)
  TAS<-TAS_load()
  meta<-TAS_drug_binder(TAS,meta_column,meta) #filter for samples that are binders or nonbinders for drug target
  df <- read_tsv(paste(GSEA_result_loc,GSEA_result,'.tsv',sep=""))
  TF= df %>% filter(pathway==TF_name) %>% select(pathway) #TF
  TARGETS=unlist(strsplit(unlist(c(df %>% filter(pathway==TF_name) %>% select(leadingEdge))),split = " "))#TF targets
  counts <- counts %>% filter(HUGO %in% TARGETS) #filter to targets
  genes<-counts$HUGO
  counts<-counts %>% select(meta$Well)#select samples
  samples <- colnames(counts)
  counts <- data.table(t(counts))
  colnames(counts)<-genes
  counts$Well <- samples
  meta <- meta[c("Well",meta_column)]
  final <- merge(counts,meta,by="Well")
  final <- subset(final, select = -c(Well) )
  name<-names(final)
  final <- final %>% gather(Gene,Expression,name[1:length(name)-1])
  p<-ggplot(final, aes_string(x="Gene", y="Expression", color=meta_column)) +
    geom_violin(trim=FALSE) +
    ggtitle(TF_name)


  heatmap.2(final)


  data(mtcars)
  x  <- as.matrix(mtcars)
  rc <- rainbow(nrow(x), start=0, end=.3)
  cc <- rainbow(ncol(x), start=0, end=.3)





  return(p)
}



#drugs that bind this, but what about strong affinity with other targets?  At least a report
#strong negative control for same target
#other ISG targets
#dirty drugs vs clean
