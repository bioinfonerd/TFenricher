
hello <- function() {
  print("Hello, world!")
}

pathway_analysis <- function(){
  reticulate::use_virtualenv( "/home/bionerd/Harvard/Artem/INDRA/env/", required=TRUE )
}

counts_qc<-function(counts){
  #remove low count wells

  counts

  #remove

}
