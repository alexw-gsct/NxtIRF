#' @useDynLib rIRFinder, .registration = TRUE 
#' @import zlibbioc
#' @importFrom Rcpp evalCpp
#' @import data.table
#' @importFrom dplyr %>%
#' @import BiocParallel
#' @importFrom rtracklayer import
#' @importFrom IRanges IRanges
#' @importFrom Biostrings getSeq DNAStringSet startIndex endIndex vmatchPattern translate
NULL