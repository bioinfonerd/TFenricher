counts_load<-function(file_name){
  #' Load in a Counts .CSV File
  #'
  #' This function allows you to load in a CSV file and makes use of the global variable locations.
  #' Using this variable, checks if the file exists
  #' If it does exist, loads in csv with progress bar and returns a tibble
  #' @param file_name File name to load in Defaults to TRUE.
  #' @keywords load
  #' @export counts tibble count matrix where columns should be sample names and first column
  #' @examples
  #' counts<-counts_load(count_file_name)
  #make sure file exists
  if(file.exists(paste(locations$counts_loc,file_name,'.csv',sep=""))){
    print("Count Table Exists")
    counts=read_csv(paste(locations$counts_loc,file_name,'.csv',sep=""),progress=show_progress())
    return(counts)
  } else{
    print("Count File Does Not Exist")
    print("Is it a CSV file?")
    print("Put just the file name without '.csv'")
  }
}

meta_load<-function(file_name){
  #' Load in a Meta .CSV File
  #'
  #' This function allows you to load in a CSV file and makes use of the global variable locations.
  #' Using this variable, checks if the file exists
  #' If it does exist, loads in csv with progress bar and returns a tibble
  #' @param file_name File name to load in Defaults to TRUE.
  #' @keywords load
  #' @export meta tibble count matrix where columns should be meta data for sample where at least one column should be sample information
  #' @examples
  #' meta<-meta_load(meta_file_name)
  #'
  #make sure file exists
  if(file.exists(paste(locations$meta_loc,file_name,'.csv',sep=""))){
    print("Meta Table Exists")
    meta=read_csv(paste(locations$meta_loc,file_name,'.csv',sep=""),progress = show_progress())
    return(meta)
  } else{
    print("Meta File Does Not Exist")
    print("Is it a CSV file?")
    print("Put just the file name without '.csv'")
  }
}

TAS_load<-function(){
  #' Load in a TAS .CSV File
  #'
  #' This function allows you to load in a CSV file and makes use of the global variable locations.
  #' Using this variable, checks if the file exists
  #' If it does exist, loads in csv with progress bar and returns a tibble
  #' @param file_name File name to load in Defaults to TRUE.
  #' @keywords load
  #' @export TAS tibble count matrix where columns should be meta data for sample where at least one column should be sample information
  #' @examples
  #' meta<-meta_load(meta_file_name)
  #'

  #make sure file exists
  if(file.exists(paste('../data/TAS_Profile/','drug_tas.csv',sep=""))){
    print("TAS Table Exists")
    TAS=read_csv(paste('../data/TAS_Profile/','drug_tas.csv',sep=""),progress=show_progress())
    return(TAS)
  } else{
    print("TAS Table Does Not Exist")
  }
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

TAS_drug_binder<-function(TAS_tibble,drug_target,meta_tibble){

  #set Drugs to lowercase
  TAS_tibble<-TAS_tibble %>% mutate_at( "name", str_to_lower )
  meta_tibble<-meta_tibble %>%  mutate_at( "Drug", str_to_lower )

  #select TAS data for drug
  cat('Selecting for Drugs that target:',drug_target,'\n')
  binders<-TAS_tibble %>% filter(symbol==drug_target & tas %in% c(1,2,3)) %>% pull(name) %>% unique()
  cat('Number of Drugs that DO bind:',length(binders),'\n')
  nonbinders<-TAS_tibble %>% filter(symbol==drug_target & tas == 10) %>% pull(name) %>% unique()
  cat('Number of Drugs that do NOT bind:',length(nonbinders),'\n')

  #select meta and counts table that matches drug independent of number of counts or any other meta data
  print('Filtering meta table for drugs')
  well_binders<- meta_tibble %>% filter(Drug %in% binders) %>% select(Well) %>% mutate( Binder = "Yes" )
  well_nonbinders<-meta_tibble %>% filter(Drug %in% nonbinders) %>% select(Well) %>% mutate( Binder = "No" )
  cat('Found',nrow(well_binders),'wells that match a binding Drug \n')
  cat('Found',nrow(well_nonbinders),'wells that match a nonbinding Drug \n')
  output<-bind_rows(well_binders,well_nonbinders)
  return(output)
}

DESeq2_run<-function(count_file_name,meta_file_name,target,condition=NA,QC=FALSE,
                     variable_filtering=FALSE,variable_filters=NA){

  #what is the input data
  cat('Count File Name:',count_file_name,
      'Meta File Name:',meta_file_name,
      'Gene Target:',target,
      'Any Conditions?:',condition=NA,
      'Any QC?:',QC=FALSE,
      'Any Variable Filtering?:',variable_filtering=FALSE,
      'Any Variable Filters?:',variable_filters=NA)

  #load data
  print('Loading Data')
  counts<-counts_load(count_file_name)
  meta<-meta_load(meta_file_name)
  TAS<-TAS_load()

  meta<-TAS_drug_binder(TAS,target,meta)

  #modify counts & meta data [TODO]: add ability to modify
  HUGO<-counts$HUGO#removes HUGO name
  counts <- counts[ , (names(counts) %in% meta$Well)] #remove unused samples
  rownames(counts)<-HUGO #depracted [TODO: fix]

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
    formula=as.formula(paste(c('~',target,collapse = ' ')))

  } else {
    formula=as.formula(paste(c('~',paste(c(target,paste(condition,sep=",")),collapse = " + ")),collapse = ' '))
  }
  cat( "Using the following formula for setting up contrasts:", str_c(formula), "\n" )

  #Differential Gene Expression
  print('Formatting DESeq2 Object')
  dds <- DESeqDataSetFromMatrix(countData=counts,colData= meta, design= formula, tidy=FALSE)

  #add gene id information
  #featureData <- data.frame(gene=HUGO)
  #mcols(dds) <- DataFrame(mcols(dds), featureData)

  cat('\n')
  print('Performing Statistics')
  dsTxi_response<-deseq2_statistics(dds)
  output<-data.frame(results(dsTxi_response))
  #make gene names a separate column
  output$genes<-rownames(output)

  # Create the differential expression directory, if needed [TODO]: make into function
  if( !dir.exists( "../Results")) dir.create("Results")
  if( !dir.exists( "../Results/Differential_Gene_Analysis")) dir.create("../Results/Differential_Gene_Analysis")
  if( !dir.exists( "../Results/Differential_Gene_Analysis/DESeq2")) dir.create("../Results/Differential_Gene_Analysis/DESeq2")
  print('writing results')
  write.table(output,file=paste("../Results/Differential_Gene_Analysis/DESeq2",paste(target,meta_file_name,"DESeq2.tsv",sep='_'),sep='/'),sep = "\t",row.names=FALSE,quote=FALSE)
}

edgeR_run<-function(count_file_name,meta_file_name,target,condition=NA,QC=FALSE,
                    variable_filtering=FALSE,variable_filters=NA){

  #what is the input data
  cat('Count File Name:',count_file_name,
      'Meta File Name:',meta_file_name,
      'Gene Target:',target,
      'Any Conditions?:',condition=NA,
      'Any QC?:',QC=FALSE,
      'Any Variable Filtering?:',variable_filtering=FALSE,
      'Any Variable Filters?:',variable_filters=NA)

  #load data
  print('Loading Data')
  counts<-counts_load(count_file_name)
  meta<-meta_load(meta_file_name)
  TAS<-TAS_load()

  meta<-TAS_drug_binder(TAS_tibble=TAS,
                        drug_target=target,
                        meta_tibble=meta)

  #modify counts & meta data [TODO]: add ability for QC to filter automatically
  HUGO<-counts$HUGO#removes HUGO name
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
    formula=as.formula(paste(c('~','Binder',collapse = ' ')))

  } else {
    formula=as.formula(paste(c('~',paste(c('Binder',paste(condition,sep=",")),collapse = " + ")),collapse = ' '))
  }

  #Differential Gene Expression
  print('Formatting EdgeR Object')
  ## Initialize the DGEList object
  y <- DGEList( counts = counts, samples = meta, genes = HUGO,remove.zeros = FALSE) #count matrix by sample group (remove rows with zeros)

  y <- calcNormFactors(y) #normalizes for RNA composition by finding a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for most genes

  ## Compose the design matrix
  #cat( "Using the following formula for setting up contrasts:", str_c(formula), "\n" )
  #mmx <- model.matrix( as.formula( formula ), data = y$samples )

  mmx <- model.matrix(formula, data = y$samples )
  cat('\n')
  print('Performing Statistics')
  y <- estimateDisp( y, mmx )
  glf <- glmFit( y, mmx ) %>% glmLRT( coef=2) #consistency with Artem's data
  ## Perform quasi-likelihood F-test
  #qlf <- glmQLFit( y, mmx ) %>% glmQLFTest( coef=2:ncol(mmx) ) #possibly better for DGE data due to high 0 data

  # Create the differential expression directory, if needed [TODO]: make into function
  if( !dir.exists( "../Results")) dir.create("Results")
  if( !dir.exists( "../Results/Differential_Gene_Analysis")) dir.create("../Results/Differential_Gene_Analysis")
  if( !dir.exists( "../Results/Differential_Gene_Analysis/EdgeR")) dir.create("../Results/Differential_Gene_Analysis/EdgeR")
  print('writing results')
  write.table(as.data.frame(topTags(glf,n=nrow(glf$table),p.value = 1)),file=paste("../Results/Differential_Gene_Analysis/EdgeR",paste(target,meta_file_name,"EdgeR.tsv",sep='_'),sep='/'),sep = "\t",row.names=FALSE,quote=FALSE)
}

