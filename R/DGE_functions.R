counts_load<-function(counts_file_path,file_name){
  counts=read.csv(paste(counts_file_path,file_name,sep=""),sep=',',row.names = 1)
  return(counts)
}

meta_load<-function(meta_file_path,file_name){
  meta=read.csv(paste(meta_file_path,file_name,sep=""),sep=',',colClasses = "character")
  return(meta)
}

deseq2_statistics<-function(dsTxi){ #input DESeqDataSet

  #statistic test:WaldTest
  dsTxi <- estimateSizeFactors(dsTxi)
  dsTxi <- estimateDispersions(dsTxi)
  dsTxi <- nbinomWaldTest(dsTxi)
  return(dsTxi)
}

counts_qc<-function(counts){
  #remove low count wells

  counts



}


DESeq2_run<-function(counts_file_path,count_file_name,meta_file_path,meta_file_name,column,condition=NA,QC=FALSE,
                     variable_filtering=FALSE,variable_filters=NA){
  print('Loading Data')
  #load data
  counts<-counts_load(counts_file_path,count_file_name)
  meta<-meta_load(meta_file_path,meta_file_name)
  print('Selecting for Column')
  #modify counts & meta data
  meta <- meta %>% filter(meta[[column]] != 'NA')
  counts <- counts[ , (names(counts) %in% meta$Well)] #remove unused samples

  if(QC==TRUE){

    meta <- meta[ , (names(meta) %in% meta$Well)] #remove unused samples
    counts <- counts[ , (names(counts) %in% meta$Well)] #remove unused samples
  }

  if(variable_filtering==TRUE){

    meta <- meta %>% filter(variable_filters)
    counts <- counts[ , (names(counts) %in% meta$Well)] #remove unused samples

  }


  #formula
  if(is.na(condition)){
    forumula=as.formula(paste(c('~',column,collapse = ' ')))

  } else {
    forumula=as.formula(paste(c('~',paste(c(column,paste(condition,sep=",")),collapse = " + ")),collapse = ' '))
  }

  #DEG
  print('Formatting DESeq2 Object')
  dds <- DESeqDataSetFromMatrix(countData=counts,colData= meta, design= forumula, tidy=FALSE)
  print('Performing Statistics')
  dsTxi_response<-deseq2_statistics(dds)
  output<-data.frame(results(dsTxi_response))
  # Create the differential expression directory, if needed [TODO]: make into function
  if( !dir.exists( "Results")) dir.create("Results")
  if( !dir.exists( "./Results/Differential_Gene_Analysis")) dir.create("./Results/Differential_Gene_Analysis")
  if( !dir.exists( "./Results/Differential_Gene_Analysis/DESeq2")) dir.create("./Results/Differential_Gene_Analysis/DESeq2")
  name <- unlist(strsplit(meta_file_name,'[.]'))[1] #add meta table name to file name output
  name <- unlist(strsplit(meta_files[2],'[.]'))[1] #add meta table name to file name output
  print('writing results')
  write.table(output,file=paste("./Results/Differential_Gene_Analysis/DESeq2",paste(column,name,".tsv",sep='_'),sep='/'),sep = "\t",quote=FALSE)
}

