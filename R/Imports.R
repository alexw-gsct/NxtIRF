#' @useDynLib NxtIRF, .registration = TRUE 
#' @import zlibbioc
#' @import data.table
#' @import shiny
#' @import shinyFiles
#' @import shinyWidgets
#' @import rhandsontable
#' @import ggplot2
#' @import plotly
#' @importFrom Rcpp evalCpp
#' @importFrom dplyr %>%
#' @importFrom AnnotationHub AnnotationHub cache
#' @importFrom assertthat assert_that
#' @importFrom BiocGenerics start end width
#' @importFrom Biostrings getSeq readDNAStringSet DNAStringSet translate replaceAmbiguities
#' @importFrom BiocParallel bpparam bplapply SnowParam MulticoreParam SerialParam
#' @importFrom DT dataTableOutput renderDataTable datatable selectRows 
#' @importFrom fst read.fst write.fst
#' @importFrom genefilter rowttests
#' @importFrom grDevices colorRampPalette
#' @importFrom GenomeInfoDb sortSeqlevels seqinfo seqlengths
#' @importFrom GenomicRanges GRanges reduce findOverlaps makeGRangesFromDataFrame mcols split strand flank setdiff seqnames psetdiff disjoin 
#' @importFrom heatmaply heatmaply
#' @importFrom matrixStats rowSds
#' @importFrom parallel detectCores
#' @importFrom rtracklayer import export TwoBitFile
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom stringr str_locate
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData assay
#' @importFrom S4Vectors coolcat metadata Rle new_SimpleList_from_list
#' @importFrom IRanges IRanges Views
NULL