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
    y[n] <- c(paste(tmp$Target, collapse = "\t"))
    n=n+1
  }

  return(data.frame(g,d, y, stringsAsFactors=FALSE))
}
