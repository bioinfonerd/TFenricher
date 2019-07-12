#returns dataframe of locations relative to vignette folder
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

#prints possible list of drugs to select
possible_drugs<-function(){
  TAS<-TAS_load()
  print(sort(unique(TAS$name)))
}

#prints out possible list of drug targets
possible_targets<-function(){
  TAS<-TAS_load()
  print(sort(unique(TAS$symbol)))
}

#merge DGE1 & DGE2 counts & meta data
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
  meta1$plate<-rep('DGE1',times = length(tmp))
  #rename them
  tmp<-meta2$Well
  for (i in 1:length(tmp)){
    tmp[i] <- paste('DGE2_',tmp[i],sep="")
  }
  meta2$Well<-tmp
  meta2$plate<-rep('DGE2',times = length(tmp))
  final<-bind_rows(meta1,meta2)
  write.csv(final,file='DGE1_DGE2_meta.csv',quote = FALSE,row.names=FALSE)
}

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

##Synapse Login
synapse_login<-function(){
  synapser::synLogin(email='nathan_johnson@hms.harvard.edu',
                     apiKey='07rRy7r0BjxoRPF6YNLFpJBj10LAekFHeCaPwpg2FYNOn4YbINW7RAVkjy7Bf8JqsxED/S5VaJYg14Wi/tDgmQ==')
}

## Define Synapse downloader(s)
syn_csv <- function( id ) {
  syn <- synExtra::synDownloader("./data", ifcollision="overwrite.local")
  syn(id) %>% read_csv(col_types=cols())}

#test all DGE conditions for DESeq2
DESeq2_run_all<-function(drug_target){
  DESeq2_synrun(drug_target=drug_target,concentration=TRUE,toxicity=TRUE)
  DESeq2_synrun(drug_target=drug_target,concentration=FALSE,toxicity=FALSE)
  DESeq2_synrun(drug_target=drug_target,concentration=TRUE,toxicity=FALSE)
  DESeq2_synrun(drug_target=drug_target,concentration=FALSE,toxicity=TRUE)
}
