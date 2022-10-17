#' Load motif counts from regulatory regions of interest
#'
#' @param motifcountsfile your motif counts file in csv form
#'
#' @return the output is simply your csv counts file loaded into R
#' @export

load_motif_counts <- function(motifcountsfile){
  motif_counts <- read.csv(motifcountsfile)
  row.names(motif_counts) <- motif_counts[,1]
  motif_counts <- motif_counts[,-1]
  return(motif_counts)
}
