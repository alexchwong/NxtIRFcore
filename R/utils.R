
#' Retrieves a list of default recommended filters
#'
#' This function returns the recommended filters. These are:\cr\cr
#' (1) Depth filter of 20,\cr\cr
#' (2) Coverage filter requiring 90% coverage in IR events\cr\cr
#' (3) Coverage filter requiring 60% coverage in AS events
#'   (i.e. Included + Excluded isoforms must cover at least 60% of all junction
#'   events across the given region)\cr\cr
#' (4) Consistency filter requring log difference of 2 (for skipped exon and
#'  mutually exclusive exon events, each junction must comprise at least 1/(2^2)
#'  = 1/4 of all reads associated with each isoform)\cr\cr
#' In all filters, we require at least 80% samples `pcTRUE = 80` from the entire
#'   dataset `minCond = "All"`.
#' Events with read depth (reads supporting either included or excluded 
#'   isoforms) lower than 20 `minDepth = 20` are excluded from filters 2,3,4.
#' @return A list of filters to be used in `apply_filters()`. Alternatively,
#'   individual filters can be run using `runFilter()`
#' @examples
#' filters = get_default_filters()
#' @md
#' @export
get_default_filters <- function() {
    filterUnit <- list()
    filterUnit$filterVars = list(
        1, 20, "All", "(none)", 80
    )
    names(filterUnit$filterVars) = 
        c("maximum", "minDepth", "minCond", "condition", "pcTRUE")
    filterUnit$filterClass = "(none)"
    filterUnit$filterType = "(none)"
    
    filters = list()
    for(i in seq_len(8)) {
        filters[[i]] = filterUnit
    }

    filters[[1]]$filterClass = "Data"    
    filters[[1]]$filterType = "Depth"    
    filters[[1]]$filterVars$minimum =  20

    filters[[2]]$filterClass = "Data"    
    filters[[2]]$filterType = "Coverage"
    filters[[2]]$filterVars$minimum =  90
    filters[[2]]$filterVars$minDepth =  5
    filters[[2]]$filterVars$EventTypes =  "IR"

    filters[[3]]$filterClass = "Data"    
    filters[[3]]$filterType = "Coverage"
    filters[[3]]$filterVars$minimum =  60
    filters[[3]]$filterVars$minDepth =  20
    filters[[3]]$filterVars$EventTypes =
        c("MXE", "SE", "AFE", "ALE", "A5SS", "A3SS")

    filters[[4]]$filterClass = "Data"    
    filters[[4]]$filterType = "Consistency"
    filters[[4]]$filterVars$maximum =  2
    filters[[4]]$filterVars$minDepth =  20
    filters[[4]]$filterVars$EventTypes = c("MXE", "SE")

    return(filters)
}


#' Constructs a matrix containing PSI values of the given ASE events
#'
#' This function takes an input SummarizedExperiment `se`, a list of 
#'   alternative splicing events `event_list`, and (optionally) a list of
#'   sample names `sample_list`. 
#' It returns a matrix containing PSI values with columns as samples and rows
#'   as ASE events.
#' 
#' @param se A NxtIRF SummarizedExperiment
#' @param event_list A character vector containing the row names of ASE events
#'   (as given by the `EventName` column of differential ASE results table 
#'   using `limma_ASE()` or `DESeq_ASE()`)
#' @param sample_list (default = `colnames(se)`) A list of sample names
#'   referring to the subset of samples in the given experiment to be included
#'   in the returned matrix
#' @param method The values to be returned (default = "PSI"). It can
#'   alternately be "logit" which returns logit-transformed PSI values, or 
#'   "Z-score" which returns Z-score-transformed PSI values
#' @param depth_threshold (default = 10) If any PSI is derived from raw values
#'   with both isoforms represented by reads below this value, it is given as
#'   `NA` as the uncertainty of PSI would be deemed too highlight
#' @param logit_max The max or min logit values to be capped at, because 
#'   `logit(0) == -Inf` and `logit(1) = Inf`,
#'   this function caps logit values using logit_max. 
#'   (Only used if `method = "logit"`)
#' @param na.percent.max (default = 0.1) The maximum number of NA values in 
#'   an event for the PSI values for each event to be returned. Most heatmap
#'   functions will spring an error if there are too many NA values in any
#'   given row. This option caps the number of NA values to avoid returning
#'   this error.
#' @return A matrix of PSI (or alternate) values, with columns as samples and
#'   rows as ASE events.
#' @md
#' @examples
#' se = NxtIRF_example_NxtSE()
#' 
#' colData(se)$treatment = rep(c("A", "B"), each = 3)
#' 
#' event_list = rowData(se)$EventName
#'
#' mat = make_matrix(se, event_list[1:10])
#' @export
make_matrix <- function(se, event_list, sample_list = colnames(se), 
    method = c("PSI", "logit", "Z-score"), 
    depth_threshold = 10, logit_max = 5, na.percent.max = 0.1) {

    method = match.arg(method)
    inc = as.matrix(assay(se, "Included")[event_list, sample_list, drop = FALSE])
    exc = as.matrix(assay(se, "Excluded")[event_list, sample_list, drop = FALSE])
    mat = inc/(inc + exc)
    mat[inc + exc < depth_threshold] = NA
    mat = mat[rowSums(is.na(mat)) < na.percent.max * ncol(mat),, drop = FALSE]
    if(method == "PSI") {
        # essentially M/Cov
        return(mat)
    } else if(method == "logit") {
        mat = qlogis(mat)
        mat[mat > logit_max] = logit_max
        mat[mat < -logit_max] = -logit_max
        return(mat)
    } else if(method == "Z-score") {
        mat = mat - rowMeans(mat)
        mat = mat / rowSds(mat)
        return(mat)
    }
    
}

#' Constructs a data frame containing average PSI values of the two contrasted 
#'   conditions
#'
#' This function takes an input SummarizedExperiment `se`, a list of 
#'   alternative splicing events `event_list`, the `condition` column 
#'   containing the experimental annotation, and the nominator `nom_DE` and
#'   denominator `denom_DE`.\cr\cr
#' Note that this function takes the geometric mean of PSI, by first converting
#'   all values to logit(PSI), taking the average logit(PSI) values of each
#'   condition, and then converting back to PSI using inv.logit.
#' 
#' @param se A NxtIRF SummarizedExperiment
#' @param event_list A character vector containing the row names of ASE events
#'   (as given by the `EventName` column of differential ASE results table
#'   using `limma_ASE()` or `DESeq_ASE()`)
#' @param condition The name of the column containing the condition values in 
#'   `colData(se)`
#' @param nom_DE The condition to be contrasted, e.g. `nom_DE = "treatment"`
#' @param denom_DE The condition to be contrasted against, e.g. 
#'   `denom_DE = "control"`
#' @param depth_threshold (default = 10) Samples with the number of reads 
#'   supporting either included or excluded isoforms below this values are 
#'   excluded
#' @param logit_max (default = 5). This function works by converting all PSI 
#'   values as logit(PSI), taking the average, then converting back to PSI 
#'   using inverse logit. Because `logit(0) == -Inf` and `logit(1) = Inf`,
#'   this function caps logit values using logit_max
#' @return A a 3 column data frame, with the first column containing 
#'   `event_list` list of ASE events, and the last 2 columns containing the 
#'   average PSI values of the nominator and denominator conditions.
#' @examples
#' se = NxtIRF_example_NxtSE()
#'
#' event_list = rowData(se)$EventName
#'
#' colData(se)$treatment = rep(c("A", "B"), each = 3)
#'
#' diag_values <- make_diagonal(se, event_list,
#'   condition = "treatment", nom_DE = "A", denom_DE = "B"
#' )
#' @md
#' @export
make_diagonal <- function(se, event_list = rownames(se), 
        condition, nom_DE, 
        denom_DE, depth_threshold = 10, logit_max = 5) {

    inc = assay(se, "Included")[event_list, ]
    exc = assay(se, "Excluded")[event_list, ]
    mat = inc/(inc + exc)
    mat[inc + exc < depth_threshold] = NA

    # use logit method to calculate geometric mean

    mat.nom = qlogis(mat[, colData(se)[,condition] == nom_DE])
    mat.denom = qlogis(mat[, colData(se)[,condition] == denom_DE])
    
    mat.nom[mat.nom > logit_max] = logit_max
    mat.denom[mat.denom > logit_max] = logit_max
    mat.nom[mat.nom < -logit_max] = -logit_max
    mat.denom[mat.denom < -logit_max] = -logit_max

    df = data.frame(EventName = event_list, 
        nom = plogis(rowMeans(mat.nom, na.rm = TRUE)),
        denom = plogis(rowMeans(mat.denom, na.rm = TRUE)))
    
    return(df)
}
