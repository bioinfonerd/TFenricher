txt_to_gmt <- function(df,description){
  g <- vector(mode="character", length=length(unique(df[,1]))) #initialize vector with plus 1 for header space
  d <- rep(description,length(unique(df[,1]))) #add description column
  y <- vector(mode="character", length=length(unique(df[,1]))) #initialize vector space for

  #column headers
  #x[1]<- 'Transcription Factor'
  #y[1]<- 'Binding Targets'

  n=1 #starting point
  for (i in unique(df[,1])){
    tmp<-subset(df,TF==i)
    g[n] <- i
    d[n] <- 'Transcription Factor'
    y[n] <- c(paste(tmp$Target, collapse = "\t"))
    n=n+1
  }

  return(data.frame(g,d, y, stringsAsFactors=FALSE))
}

gmx_to_gmt <- function(){
  tmp <- read.csv('../data/isg.gmx')
  tmp<-as.data.frame(t(tmp))
  write.table(tmp,file='../data/isg.gmt',sep="\t",col.names=FALSE,quote=FALSE)
}

append_gmt<-function(){
  df <- read.csv('C:/Users/Nathan/Documents/@Harvard/Artem/Transcription_Factor_Analysis/git/Transcription-Factor-Databases/Ttrust_v2/trrust_rawdata.human.tsv',
                 head = FALSE,sep='\t')
  colnames(df) <- c('TF','Target','Interaction','Pubmed_ID')

  df <- df %>% distinct(TF,Target) #drop duplicates

    g <- vector(mode="character", length=length(unique(df[,1]))+1) #initialize vector with plus 1 for header space
    d <- rep('Transcription Factor',length(unique(df[,1]))+1) #add description column
    y <- vector(mode="character", length=length(unique(df[,1]))+1) #initialize vector space for

    n=1 #starting point
    for (i in unique(df[,1])){
      tmp<-subset(df,TF==i)
      g[n] <- i
      d[n] <- 'Transcription Factor'
      y[n] <- c(paste(tmp$Target, collapse = "\t"))
      n=n+1
    }

    #add ISG


    #add ISG to gmt
    add <- fgsea::gmtPathways('./data/isg.gmt')
    g[n] <- 'ISG'
    d[n] <- 'Interferon Stimulated Genes'
    y[n] <- c(paste(unlist(add[1]), collapse = "\t"))

    write.table(data.frame(g,d, y, stringsAsFactors=FALSE),file='all_plusISG.gmt',sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

}
