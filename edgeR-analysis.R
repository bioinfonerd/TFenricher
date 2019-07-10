#edgeR: DEG on each subset ------------------------------------------------------------------

#counts=read.csv('../data/postqc-counts.csv',sep=',',row.names = 1)
#meta=read.csv('../data/postqc-meta.csv',sep=',')

#filter meta and counts to selection
#meta <- meta %>% filter(Drug == 'DMSO' | Drug == 'dsRNAmi') #select meta table
counts <- counts[ , (names(counts) %in% meta$Well)] #remove unused columns

#remove unused factor levels
meta<-droplevels(meta) 
counts<-droplevels(counts)


## Initialize the DGEList object
y <- DGEList( counts = counts, samples = meta, remove.zeros = TRUE) #count matrix by sample group (remove rows with zeros)
y <- calcNormFactors(y) #normalizes for RNA composition by finding a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for most genes

## Compose the design matrix
f <- str_c( "~", 'Drug' )
cat( "Using the following formula for setting up contrasts:", f, "\n" )
mmx <- model.matrix( as.formula( f ), data = y$samples )
y <- estimateDisp( y, mmx )

## Perform quasi-likelihood F-test
qlf <- glmQLFit( y, mmx ) %>% glmQLFTest( coef=2:ncol(mmx) )

## Retrieve the results matrix
RR <- topTags(qlf) %>% as.data.frame()

## Create the differential expression directory, if needed
if( !dir.exists( "edgeR" ) ) dir.create( "edgeR" )

write.table(RR,file=paste('./edgeR',gsub("/","_",paste(group1,group2,"edgeR.results.tsv",sep='_')),sep='/'),sep = "\t",quote=FALSE)
write.table(as.data.frame(qlf$table),file=paste('./edgeR',gsub("/","_",paste(group1,group2,"_full_edgeR.results.tsv",sep='_')),sep='/'),sep = "\t",quote=FALSE)


edgeR_run<-function(count_file_name,meta_file_name,target,condition=NA,QC=FALSE,
                     variable_filtering=FALSE,variable_filters=NA){

  #load data
  print('Loading Data')
  counts<-counts_load(count_file_name)
  meta<-meta_load(meta_file_name)
  TAS<-TAS_load()

  #select TAS data for drug
  cat('Selecting for Drugs that target:',target,'\n')
  binders<-TAS %>% filter(symbol==target & tas %in% c(1,2,3))
  cat('Number of Drugs that DO bind:',nrow(binders),'\n')
  nonbinders<-TAS %>% filter(symbol==target & tas == 10) %>% select(name)
  cat('Number of Drugs that do NOT bind:',nrow(nonbinders),'\n')

  #select meta and counts table that matches drug independent of number of counts or any other meta data
  cat('Found',nrow(meta %>% filter(Drug %in% binders$name)),'wells that match a binding Drug \n')
  cat('Found',nrow(meta %>% filter(Drug %in% nonbinders$name)),'wells that match a nonbinding Drug \n')
  print('Filtering meta table for drugs') #[TODO]: should be doable using dyplr...
  meta_mod<- meta %>% filter(Drug %in% binders$name) %>% select(Well)
  meta_mod[target]<-'binder'
  tmp<-meta %>% filter(Drug %in% nonbinders$name) %>% select(Well)
  tmp[target]<-'nonbinder'
  meta_mod<-rbind(meta_mod,tmp)
  meta<-merge(meta,meta_mod,by='Well')

  #modify counts & meta data [TODO]: add ability for QC to filter automatically
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
    formula=as.formula(paste(c('~',target,collapse = ' ')))

  } else {
    formula=as.formula(paste(c('~',paste(c(target,paste(condition,sep=",")),collapse = " + ")),collapse = ' '))
  }

  #Differential Gene Expression
  print('Formatting DESeq2 Object')
  dds <- DESeqDataSetFromMatrix(countData=counts,colData= meta, design= formula, tidy=FALSE)
  print('\n')
  print('Performing Statistics')
  dsTxi_response<-deseq2_statistics(dds)
  output<-data.frame(results(dsTxi_response))
  # Create the differential expression directory, if needed [TODO]: make into function
  if( !dir.exists( "Results")) dir.create("Results")
  if( !dir.exists( "./Results/Differential_Gene_Analysis")) dir.create("./Results/Differential_Gene_Analysis")
  if( !dir.exists( "./Results/Differential_Gene_Analysis/EdgeR")) dir.create("./Results/Differential_Gene_Analysis/DESeq2")
  name <- unlist(strsplit(meta_file_name,'[.]'))[1] #add meta table name to file name output
  name <- unlist(strsplit(meta_files[2],'[.]'))[1] #add meta table name to file name output
  print('writing results')
  write.table(output,file=paste("./Results/Differential_Gene_Analysis/DESeq2",paste(column,name,".tsv",sep='_'),sep='/'),sep = "\t",quote=FALSE)
}

