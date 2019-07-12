#filter wells by plate index specific threshold (toxicity)
filterplateindexWells <- function(meta,counts){

  #filter based on plate index ('DGE1') specific threshold
  DGE1_thresh = 3e4
  wKeep_DGE1 <- counts %>% summarize_at( vars(-HUGO), sum ) %>% gather( Well, TotalCounts ) %>%
    filter( TotalCounts >= DGE1_thresh ) %>% filter(Well %in%  filter(meta,plate=='DGE1')$Well) %>% pull(Well)

  #filter based on plate index ('DGE2') specific threshold
  DGE2_thresh = 8.5e4
  wKeep2_DGE2 <- counts %>% summarize_at( vars(-HUGO), sum ) %>% gather( Well, TotalCounts ) %>%
    filter( TotalCounts >= DGE2_thresh ) %>% filter(Well %in%  filter(meta,plate=='DGE2')$Well) %>% pull(Well)
  return(c(wKeep_DGE1,wKeep2_DGE2))
}

#output information given to DGE function
input_info<-function(count_file_name,
                     meta_file_name,
                     target,condition,
                     variable_filtering,
                     variable_filters){
    #what is the input data
  cat('Count File Name:',count_file_name,
      'Meta File Name:',meta_file_name,
      'Gene Target:',target,
      'Any Conditions?:',condition=NA,
      'Any QC?:',QC=FALSE,
      'Any Variable Filtering?:',variable_filtering=FALSE,
      'Any Variable Filters?:',variable_filters=NA)
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

filter_meta_by_TAS<-function(TAS_tibble,drug_target,meta_tibble,condition){

  #set Drugs to lowercase
  TAS_tibble<-TAS_tibble %>% mutate_at( "name", str_to_lower )
  meta_tibble<-meta_tibble %>%  mutate_at( "Drug", str_to_lower )

  #select TAS data for drug
  cat('Selecting for Drugs that target:',drug_target,'\n')
  binders<-TAS_tibble %>% filter(symbol==drug_target & tas %in% c(1,2,3)) %>% pull(name) %>% unique()
  cat('Number of Drugs that DO bind:',length(binders),'\n')
  nonbinders<-TAS_tibble %>% filter(symbol==drug_target & tas == 10) %>% pull(name) %>% unique()
  cat('Number of Drugs that do NOT bind:',length(nonbinders),'\n')

  #if there is a condition added to formula keep that column, otherwise remove
  if(condition==FALSE){
    #select meta and counts table that matches drug independent of number of counts or any other meta data
    print('Filtering meta table for drugs')
    well_binders<- meta_tibble %>% filter(Drug %in% binders) %>% select(Well) %>% mutate( Binder = "Yes" )
    well_nonbinders<-meta_tibble %>% filter(Drug %in% nonbinders) %>% select(Well) %>% mutate( Binder = "No" )
    cat('Found',nrow(well_binders),'wells that match a binding Drug \n')
    cat('Found',nrow(well_nonbinders),'wells that match a nonbinding Drug \n')
    output<-bind_rows(well_binders,well_nonbinders)
    return(as.data.frame(output))
  } else {
    #select meta and counts table that matches drug independent of number of counts or any other meta data
    print('Filtering meta table for drugs')
    well_binders<- meta_tibble %>% filter(Drug %in% binders) %>% select(Well,condition) %>% mutate( Binder = "Yes" )
    well_nonbinders<-meta_tibble %>% filter(Drug %in% nonbinders) %>% select(Well,condition) %>% mutate( Binder = "No" )
    cat('Found',nrow(well_binders),'wells that match a binding Drug \n')
    cat('Found',nrow(well_nonbinders),'wells that match a nonbinding Drug \n')
    output<-bind_rows(well_binders,well_nonbinders)
    return(BiocGenerics::as.data.frame(output))
  }
}

#DESeq2 run with different conditions posisble for input
DESeq2_synrun<-function(count_file_syn_id='syn20494174',
                     meta_file_syn_id='syn20495336',
                     TAS_file_syn_id='syn18268627',
                     drug_target,
                     condition=FALSE,
                     concentration=FALSE,
                     toxicity=FALSE){

  #Logging into Synapse
  print('Logging Into Synapse')
  synapse_login()
  #load data
  print('Loading Data')
  counts<-syn_csv(count_file_syn_id)
  meta<-syn_csv(meta_file_syn_id)
  TAS<-syn_csv(TAS_file_syn_id)
  #filter by toxicity if label set
  if(toxicity==TRUE){meta <- filter(meta,Well %in% filterplateindexWells(meta,counts))}
  #filter by concentration if label set
  if(concentration==TRUE){condition="Concentration"}
  #filtering for drug target
  meta<-filter_meta_by_TAS(TAS,drug_target,meta,condition)

  #filter which counts data based on meta & and change to dataframe
  counts<-BiocGenerics::as.data.frame(counts)
  rownames(counts)<-counts$HUGO
  counts <- counts[ , (names(counts) %in% meta$Well)]
  #if to add condition to statistical test
  if(is.na(condition)){formula=as.formula(paste(c('~','Binder',collapse = ' ')))} else
    {formula=as.formula(paste(c('~',paste(c('Binder',paste(condition,sep=",")),collapse = " + ")),collapse = ' '))}
  cat( "Using the following formula for setting up contrasts:", str_c(formula), "\n" )

  #Differential Gene Expression
  print('Formatting DESeq2 Object')
  dds <- DESeqDataSetFromMatrix(countData=counts,colData= meta, design= formula, tidy=FALSE)
  cat('\n')
  print('Performing Statistics')
  dsTxi_response<-deseq2_statistics(dds)
  output<-data.frame(results(dsTxi_response))
  #make gene names a separate column
  output$Genes<-rownames(output)

  # Create the differential expression directory, if needed [TODO]: make into function
  if( !dir.exists( "../Results")) dir.create("Results")
  if( !dir.exists( "../Results/Differential_Gene_Analysis")) dir.create("../Results/Differential_Gene_Analysis")
  if( !dir.exists( "../Results/Differential_Gene_Analysis/DESeq2")) dir.create("../Results/Differential_Gene_Analysis/DESeq2")

  print('writing results')
  #output file name: '"JAK2_conT_toxF" == JAK2 drug target, concentration true, toxicity false
  if(concentration==FALSE){concentration='conF'}else{concentration='conT'}
  if(toxicity==FALSE){toxicity='toxF'}else{toxicity='toxT'}
  filename<-paste(drug_target,concentration,toxicity,"DESeq2.tsv",sep='_')
  write.table(output,file=paste("../Results/Differential_Gene_Analysis/DESeq2",filename,sep='/'),sep = "\t",row.names=FALSE,quote=FALSE)
}

#[TODO]: change to synapse based systems
edgeR_run<-function(count_file_name,meta_file_name,target,condition=NA,QC=FALSE,
                    variable_filtering=FALSE,variable_filters=NA){

  #what is the input data ([TODO]: everything above DGE analysis could be a function)
  cat('Count File Name:',count_file_name,
      'Meta File Name:',meta_file_name,
      'Gene Target:',target,
      'Any Conditions?:',condition=NA,
      'Any QC?:',QC=FALSE,
      'Any Variable Filtering?:',variable_filtering=FALSE,
      'Any Variable Filters?:',variable_filters=NA)

  #load data
  print('Loading Data')
  counts_N<-counts_load(count_file_name)
  meta<-meta_load(meta_file_name)
  TAS<-TAS_load()
  print('Filter for Drug')
  meta<-TAS_drug_binder(TAS_tibble=TAS,drug_target=target,meta_tibble=meta)
  meta<-as.data.frame(meta)
  #modify counts & meta data [TODO]: add ability for QC to filter automatically
  HUGO_N<-counts_N$HUGO#removes HUGO name
  counts_N <- counts_N[ , (names(counts_N) %in% meta$Well)] #remove unused samples
  #rownames(counts_N)<-HUGO_N #setting rownames on a tibble is depracated
  counts_N<-as.data.frame(counts_N) #to match Artem
  rownames(counts_N)<-HUGO_N

  #modify counts & meta data [TODO]: add ability for QC and condition filtering automatically
  #if(QC==TRUE){

  #  meta <- meta[ , (names(meta) %in% meta$Well)] #remove unused samples
  #  counts <- counts[ , (names(counts) %in% meta$Well)] #remove unused samples
  #}

  #if(variable_filtering==TRUE){

  #  meta <- meta %>% filter(variable_filters)
  #  counts <- counts[ , (names(counts) %in% meta$Well)] #remove unused samples

  #}

  #will assign conditions to consider dispersion analysis based on condition input (must be column names)
  #[TODO]: add column name verifiction from meta data
  #if(is.na(condition)){
  #  formula=as.formula(paste(c('~','Binder',collapse = ' ')))
#
 # } else {
  #  formula=as.formula(paste(c('~',paste(c('Binder',paste(condition,sep=",")),collapse = " + ")),collapse = ' '))
  #}

  #formula
  formula=as.formula(paste(c('~','Binder',collapse = ' ')))
  ## Initialize the DGEList object
  y <- DGEList( counts = counts_N, samples = meta, remove.zeros = FALSE) #count matrix by sample group (remove rows with zeros and changed when gene names are added) (Artem's)

  #y <- DGEList( counts = counts_N, samples = meta, genes = HUGO_N,remove.zeros = TRUE) #count matrix by sample group (remove rows with zeros) (original)
  y <- calcNormFactors(y) #normalizes for RNA composition by finding a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for most genes
  mmx_N <- model.matrix(formula, data = y$samples )
  y <- estimateDisp( y, mmx_N )
  #glf <- glmFit( y, mmx_N ) %>% glmLRT( coef=2) #consistency with Artem's data
  ## Perform quasi-likelihood F-test
  qlf <- glmQLFit( y, mmx ) %>% glmQLFTest( coef=2:ncol(mmx) ) #possibly better for DGE data due to high 0 data

  # Create the differential expression directory, if needed [TODO]: make into function
  #output_N<-as.data.frame(topTags(glf,n=nrow(glf$table),p.value = 1))
  output_N<-edgeR::topTags( gf, nrow(counts_N) ) %>% as.data.frame %>% rownames_to_column( "Gene" )

  # Create the differential expression directory, if needed [TODO]: make directory structure into function
  if( !dir.exists( "../Results")) dir.create("Results")
  if( !dir.exists( "../Results/Differential_Gene_Analysis")) dir.create("../Results/Differential_Gene_Analysis")
  if( !dir.exists( "../Results/Differential_Gene_Analysis/EdgeR")) dir.create("../Results/Differential_Gene_Analysis/EdgeR")
  print('writing results')
  write.table(output_N,file=paste("../Results/Differential_Gene_Analysis/EdgeR",paste(target,meta_file_name,"EdgeR.tsv",sep='_'),sep='/'),sep = "\t",row.names=FALSE,quote=FALSE)
  return(output_N)

}

#[TODO]: change to synapse based systems
# DESeq2_run<-function(count_file_name,meta_file_name,target,condition=NA,QC=FALSE,
#                      variable_filtering=FALSE,variable_filters=NA){
#
#   #what is the input data
#   cat('Count File Name:',count_file_name,
#       'Meta File Name:',meta_file_name,
#       'Gene Target:',target,
#       'Any Conditions?:',condition=NA,
#       'Any QC?:',QC=FALSE,
#       'Any Variable Filtering?:',variable_filtering=FALSE,
#       'Any Variable Filters?:',variable_filters=NA)
#
#   #load data
#   print('Loading Data')
#   counts<-counts_load(count_file_name)
#   meta<-meta_load(meta_file_name)
#   TAS<-TAS_load()
#
#   meta<-TAS_drug_binder(TAS,target,meta)
#
#   #modify counts & meta data [TODO]: add ability to modify
#   HUGO<-counts$HUGO#removes HUGO name
#   counts <- counts[ , (names(counts) %in% meta$Well)] #remove unused samples
#   rownames(counts)<-HUGO #depracted [TODO: fix]
#
#   if(QC==TRUE){
#
#     meta <- meta[ , (names(meta) %in% meta$Well)] #remove unused samples
#     counts <- counts[ , (names(counts) %in% meta$Well)] #remove unused samples
#   }
#
#   if(variable_filtering==TRUE){
#
#     meta <- meta %>% filter(variable_filters)
#     counts <- counts[ , (names(counts) %in% meta$Well)] #remove unused samples
#
#   }
#
#   #formula
#   if(is.na(condition)){
#     formula=as.formula(paste(c('~',target,collapse = ' ')))
#
#   } else {
#     formula=as.formula(paste(c('~',paste(c(target,paste(condition,sep=",")),collapse = " + ")),collapse = ' '))
#   }
#   cat( "Using the following formula for setting up contrasts:", str_c(formula), "\n" )
#
#   #Differential Gene Expression
#   print('Formatting DESeq2 Object')
#   dds <- DESeqDataSetFromMatrix(countData=counts,colData= meta, design= formula, tidy=FALSE)
#
#   #add gene id information
#   #featureData <- data.frame(gene=HUGO)
#   #mcols(dds) <- DataFrame(mcols(dds), featureData)
#
#   cat('\n')
#   print('Performing Statistics')
#   dsTxi_response<-deseq2_statistics(dds)
#   output<-data.frame(results(dsTxi_response))
#   #make gene names a separate column
#   output$genes<-rownames(output)
#
#   # Create the differential expression directory, if needed [TODO]: make into function
#   if( !dir.exists( "../Results")) dir.create("Results")
#   if( !dir.exists( "../Results/Differential_Gene_Analysis")) dir.create("../Results/Differential_Gene_Analysis")
#   if( !dir.exists( "../Results/Differential_Gene_Analysis/DESeq2")) dir.create("../Results/Differential_Gene_Analysis/DESeq2")
#   print('writing results')
#   write.table(output,file=paste("../Results/Differential_Gene_Analysis/DESeq2",paste(target,meta_file_name,"DESeq2.tsv",sep='_'),sep='/'),sep = "\t",row.names=FALSE,quote=FALSE)
# }
