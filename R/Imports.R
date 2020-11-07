#' @useDynLib NxtIRF, .registration = TRUE 
#' @import data.table
#' @import shiny
#' @import shinyFiles
#' @import shinyWidgets
#' @import rhandsontable
#' @import ggplot2
#' @importFrom boot logit inv.logit
#' @importFrom methods as is
#' @importFrom stats as.formula model.matrix qt runif
#' @importFrom utils download.file packageVersion setTxtProgressBar txtProgressBar
#' @importFrom Rcpp evalCpp
#' @importFrom dplyr %>%
#' @importFrom AnnotationHub AnnotationHub cache
#' @importFrom assertthat assert_that
#' @importFrom BiocGenerics start end width
#' @importFrom Biostrings getSeq readDNAStringSet DNAStringSet translate replaceAmbiguities
#' @importFrom BiocParallel bpparam bplapply SnowParam MulticoreParam SerialParam
#' @importFrom DT datatable selectRows 
#' @importFrom fst read.fst write.fst
#' @importFrom genefilter rowttests
#' @importFrom grDevices colorRampPalette
#' @importFrom GenomeInfoDb sortSeqlevels seqinfo seqlengths
#' @importFrom GenomicRanges GRanges reduce findOverlaps makeGRangesFromDataFrame makeGRangesListFromDataFrame mcols split strand flank setdiff seqnames psetdiff disjoin mcols<- strand<-
#' @importFrom heatmaply heatmaply
#' @importFrom matrixStats rowSds
#' @importFrom parallel detectCores
#' @importFrom plotly config layout plotlyOutput event_data ggplotly plotlyProxy plotlyProxyInvoke renderPlotly subplot
#' @importFrom rtracklayer import export TwoBitFile
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom stringr str_locate
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData assay
#' @importFrom S4Vectors coolcat metadata Rle metadata<- SimpleList
#' @importFrom IRanges IRanges Views
NULL