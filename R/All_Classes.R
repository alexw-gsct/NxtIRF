
#' The NxtAnnotation Class
#'
#' The class that stores NxtIRF reference
#'
#' @name NxtAnnotation
#' @field chrOrder Chromosome names
#' @field Genes Gene annotations (GRanges)
#' @field Trascripts Transcript annotations (GRanges)
#' @field Exons Exon annotations (GRanges)
#' @field Introns Introns (data.table). Contains duplicate introns for each transcripts that share the same genomic interval. For unique introns, use Introns_dedup
#' @field Introns_dedup Introns, unique for each genomic interval (data.table).
#' @field Splice Splice annotations (data.table)
#' @importFrom GenomicRanges GRanges
#' @importFrom data.table data.table
setClass("NxtAnnotation",
    slots = c(RefPath = "character",
        RefFasta = "character",
        RefGTF = "character",
		chrOrder = "character",
		Genes = "GRanges",
		Transcripts ="GRanges",
		Exons = "GRanges",
		Introns = "data.table",
		Introns_dedup = "data.table",
		Splice = "data.table")
)

setMethod("show", "NxtAnnotation",
    function(object)
{
    cat("class:", class(object), "\n\n")

    cat("Reference Path:", object@RefPath, "\n")
	S4Vectors::coolcat("Chromosomes(%d): %s\n", object@chrOrder)
	S4Vectors::coolcat("Genes(%d): %s\n", object@Genes$gene_name)
	S4Vectors::coolcat("Transcripts(%d): %s\n", object@Transcripts$transcript_name)
	S4Vectors::coolcat("Annotated Introns: %d\n", object@Introns_dedup$EventName)

	S4Vectors::coolcat("Mutually Exclusive Exons: %d\n", object@Splice$EventName[object@Splice$EventType == "MXE"])
	S4Vectors::coolcat("Skipped Exons: %d\n", object@Splice$EventName[object@Splice$EventType == "SE"])
	S4Vectors::coolcat("Alternative 5'-Splice Site: %d\n", object@Splice$EventName[object@Splice$EventType == "A5SS"])
	S4Vectors::coolcat("Alternative 3'-Splice Site: %d\n", object@Splice$EventName[object@Splice$EventType == "A3SS"])
	S4Vectors::coolcat("Alternative First Exons: %d\n", object@Splice$EventName[object@Splice$EventType == "AFE"])
	S4Vectors::coolcat("Alternative Last Exons: %d\n", object@Splice$EventName[object@Splice$EventType == "ALE"])	
})


#' The NxtProject Class
#'
#' The main class used by NxtIRF to store project information.
#'
#' @name NxtProject
#' @field metadata list, contains Project, Data, and Bam directory paths, as well as given name of data Header
#' @field metrics list, contains metrics that describe the nature of the data (stranded vs unstranded) as well as various parameters from which the raw IRFinder data is compiled. See NxtIRF.ProcessData()
#' @field Expr list, contains list of sample annotations. The entire data set is retrieved as Expr$global or Expr[["global"]]. Subset analysis is also possible by sub-setting the Expr$global data frame.
#' @field Data list, contains compiled raw data loaded into memory.
#' @field nds list, contains refined countData with accompanying colData, for parsing onto other package methods for differential analysis such as DESeq2, edgeR or limma
#' @field benchmarks list, contains system times for benchmark purposes
#' @importFrom S4Vectors metadata metadata<- SimpleList
setClass("NxtProject",
    slots = c(metadata = "SimpleList",
		metrics = "SimpleList",
		Annotation="NxtAnnotation",
		Data = "SimpleList",
		nds = "SimpleList",
		benchmarks = "SimpleList")
)

setMethod("show", "NxtProject",
    function(object)
{
    cat("class:", class(object), "\n\n")

    cat("Project name:", object@metadata$ProjectName, "\n")
    cat("Source IRFinder output:", object@metadata$SamplePath, "\n")

	cat("\n")
	
	cat("Metrics:\n")
	
    cat("Reference Path:", object@metrics$ReferencePath, "\n")
    cat("Stranded data:", object@metrics$Stranded, "\n")
    cat("Unannotated Junction Threshold (minimum splice reads):", object@metrics$UnannotatedJnThreshold, "\n")
	
	cat("\n")
	S4Vectors::coolcat("Samples(%d): %s\n", object@Data$Expr$SampleName)
})

