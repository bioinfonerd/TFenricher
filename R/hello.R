# Hello, world!
#`
#`  This is an example function named 'hello'
#` which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


#' Add together two numbers.
#'
#' @param x A number.
#' @param y A number.
#' @return The sum of \code{x} and \code{y}.
#' @examples
#' add(1, 1)
#' add(10, 1)
add <- function(x, y) {
  x + y
}



hello <- function() {
  print("Hello, world!")
}


pathway_analysis <- function(){
  reticulate::use_virtualenv( "/home/bionerd/Harvard/Artem/INDRA/env/", required=TRUE )
}






library( indRa )

## If using virtualenv
reticulate::use_virtualenv( "/path/to/virtual/env", required=TRUE )

## Access to INDRA is available through indra() function
indra()
# Module(indra)


dijkstra()

PW <- dijkstra( "JAK2", trgts=c("NFKB1", "STAT1", "STAT2", "STAT3", "IRF1", "IRF3") )
P <- with(PW, setNames(Path, Gene))

plotPaths( PW$Path ) + ggthemes::scale_color_few()

