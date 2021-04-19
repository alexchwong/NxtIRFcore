# NxtSE Class functions

#' The NxtSE class methods
#'
#' The NxtSE class inherits from the \linkS4class{SummarizedExperiment} 
#' class and is defined to represent that it is constructed from `MakeSE()`.
#' @name NxtSE-class
#' @export
setClass("NxtSE",
    slots=c(int_elementMetadata = "DataFrame",
        int_colData = "DataFrame",
        int_metadata = "list"),
    contains = "SummarizedExperiment"
)


