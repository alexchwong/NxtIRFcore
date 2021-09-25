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


#' The NxtFilterVars class
#'
#' This is a simple class to contain a number of filter settings
#' for the filtering of alternative splicing and intron retention
#' @name NxtFilter-class
#' @export
setClass("NxtFilter",
    slots = c(
        filterClass = "character",
        filterType = "character",
    
        pcTRUE = "numeric",
        minimum = "numeric",
        maximum = "numeric",
        minDepth = "numeric",
        condition = "character",
        minCond = "numeric",
        EventTypes = "character"
    )
)