location_info<-function(){
  locations <- data.frame(meta_loc=NA,
                          meta_files=NA,
                          counts_loc=NA,
                          counts_files=NA,
                          GMT_loc=NA,
                          GMT_files=NA)
  #assign where meta, count, and GMT information is stored
  #meta
  locations$meta_loc <-'../data/DGE_data/'
  locations$meta_files <- list(list.files(locations$meta_loc,pattern = "*meta*.csv")) #assign a list of data
  #counts
  locations$counts_loc <- '../data/DGE_data/'
  locations$counts_files <- list(list.files(locations$counts_loc,pattern = "*counts*.csv"))
  #GMT files
  locations$GMT_loc<-'../data/TF_GSEA_GMT_FILES/'
  locations$GMT_files <- list(list.files(locations$GMT_loc,pattern = "*.gmt"))
  return(locations)
}

possible_drugs<-function(){
  TAS<-TAS_load()
  print(sort(unique(TAS$name)))
}

possible_targets<-function(){
  TAS<-TAS_load()
  print(sort(unique(TAS$symbol)))
}

merge_DGE_data<-function(){
  #load in data
  setwd('../data/DGE_data')
  dge1<-read_csv('DGE1_postqc-counts.csv')
  dge2<-read_csv('DGE2_postqc-counts.csv')
  meta1<-read_csv('DGE1_postqc-meta.csv')
  meta2<-read_csv('DGE2_postqc-meta.csv')
  #rename count them
  tmp<-colnames(dge1)
  for (i in 1:length(tmp)){
    tmp[i] <- paste('DGE1_',tmp[i],sep="")
  }
  tmp[1]<-"HUGO"
  colnames(dge1)<-tmp
  #rename them
  tmp<-colnames(dge2)
  for (i in 1:length(tmp)){
    tmp[i] <- paste('DGE2_',tmp[i],sep="")
  }
  tmp[1]<-"HUGO"
  colnames(dge2)<-tmp
  final<-inner_join(dge1,dge2,by='HUGO')
  write.csv(final,file='DGE1_DGE2_counts.csv',quote = FALSE,row.names=FALSE)
  #rename meta
  tmp<-meta1$Well
  for (i in 1:length(tmp)){
    tmp[i] <- paste('DGE1_',tmp[i],sep="")
  }
  meta1$Well<-tmp
  #rename them
  tmp<-meta2$Well
  for (i in 1:length(tmp)){
    tmp[i] <- paste('DGE2_',tmp[i],sep="")
  }
  meta2$Well<-tmp
  final<-bind_rows(meta1,meta2)
  write.csv(final,file='DGE1_DGE2_meta.csv',quote = FALSE,row.names=FALSE)
}

add_ISG_data<-function(){

  df<-read_csv('../data/TF_GSEA_GMT_FILES/all.gmt')

  ../data/isg.

}
