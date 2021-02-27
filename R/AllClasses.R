# NxtSE Class functions

#' @export
#' @rdname NxtSE
#' @importFrom utils packageVersion
#' @importFrom stats setNames
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom S4Vectors DataFrame
setClass("NxtSE",
    slots=c(int_elementMetadata = "DataFrame",
        int_colData = "DataFrame",
        int_metadata = "list"),
    contains = "SummarizedExperiment"
)

setMethod("show", "NxtSE",
    function(object)
{
    callNextMethod(object)
})

#' @exportMethod coerce
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setAs("SummarizedExperiment", "NxtSE", function(from) {
    .se_to_nxtse(from)
})

#' The NxtSE class
#'
#' The NxtSE class inherits from the \linkS4class{SummarizedExperiment} 
#' class and is defined to 
#' represent that it is constructed from `MakeSE()`.
#' @param ... Additional parameters to pass into `SummarizedExperiment()`
#' @return A validated NxtSE object
#' @docType class
#' @aliases
#' coerce,SummarizedExperiment,NxtSE-method
#' @md
#' @export
NxtSE <- function(...) {
    se <- SummarizedExperiment(...)
    
    .se_to_nxtse(se)
}

#' @importFrom S4Vectors DataFrame SimpleList
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom methods new
#' @importFrom BiocGenerics nrow ncol
.se_to_nxtse <- function(se) {
    # Simple test to make sure se is a SummarizedExperiment made by NxtIRF::MakeSE()
    if(!all(c("Included", "Excluded", "Depth", "Coverage", "minDepth") %in% 
        names(assays(se)))) {
        warning(paste("Object was not created by MakeSE(),",
            "returning SummarizedExperiment object instead"))
        return(se)
    }
    if(!all(c("Up_Inc", "Down_Inc", "Up_Exc", "Down_Exc") %in%
            names(S4Vectors::metadata(se)))) {
        warning(paste("Object was not created by MakeSE(),",
            "returning SummarizedExperiment object instead"))
        return(se)
    }
    out <- new("NxtSE", se, 
        int_elementMetadata=new("DFrame", nrows=nrow(se)),
        int_colData=new("DFrame", nrows=ncol(se)))

    return(out)
}