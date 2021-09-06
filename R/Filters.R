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
#'  = 1/4 of all reads associated with each isoform).
#'  For retained introns, the exon-intron overhangs must not differ by 1/4\cr\cr
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
        c("MXE", "SE", "AFE", "ALE", "A5SS", "A3SS", "RI")

    filters[[4]]$filterClass = "Data"    
    filters[[4]]$filterType = "Consistency"
    filters[[4]]$filterVars$maximum =  2
    filters[[4]]$filterVars$minDepth =  20
    filters[[4]]$filterVars$EventTypes = c("MXE", "SE", "RI")

    return(filters)
}


#' Filtering for IR and Alternative Splicing Events
#'
#' This function implements filtering of IR or AS events based on customisable
#' criteria
#' 
#' @details
#'   \strong{Annotation Filters}\cr\cr 
#'     \strong{Protein_Coding}: Filters for alternative splicing or IR events 
#'       within protein reading frames. No additional parameters required.\cr\cr
#'     \strong{NMD_Switching}: Filters for events in which one isoform is a
#'       predicted NMD substrate.\cr\cr 
#'     \strong{Transcript_Support_Level}: filters for events in which both
#'       isoforms have a TSL level below or equal to filterVars$minimum\cr\cr 
#'   \strong{Data Filters}\cr\cr 
#'     \strong{Depth}: Filters IR or alternative splicing events of transcripts
#'       that are "expressed" with adequate \code{Depth} as calculated by the
#'       sum of all splicing and IR reads spanning the event. Events with 
#'       \code{Depth} below filterVars$minimum are excluded\cr\cr
#'     \strong{Coverage}: Coverage means different things to IR and alternative
#'       splicing.\cr\cr 
#'     For \emph{IR}, Coverage refers to the percentage of the measured intron
#'       covered with reads. Introns of samples with an IntronDepth above 
#'       \code{filterVars$minDepth} are assessed, with introns with coverage 
#'       below \code{filterVars$minimum} are excluded.\cr\cr 
#'     For \emph{Alternative Splicing}, Coverage refers to the percentage of all
#'       splicing events observed across the genomic region that is compatible
#'       with either the included or excluded event. This prevents NxtIRF from 
#'       doing differential analysis between two minor isoforms. Instead of
#'       IntronDepth, in AS events NxtIRF considers events where the spliced
#'       reads from both exonic regions exceed \code{filterVars$minDepth}.
#'       Then, events with a splicing coverage below \code{filterVars$minimum}
#'       are excluded. We recommend testing IR events for > 90% coverage and AS
#'       events for > 60% coverage as given in the default filters which can be
#'       accessed using \code{\link{get_default_filters}}\cr\cr  
#'     \strong{Consistency}: Skipped exons (SE) and mutually exclusive exons
#'       (MXE) comprise reads of two contiguous splice junctions (for the
#'       included casette exon). Summating counts from both junctions is
#'       misleading as there may be overlapping events (e.g. alternate first 
#'       / last exons) that only rely on one splice event. To ensure the SE /
#'       MXE is the dominant event, we require both splice junctions to have
#'       comparable counts.\cr\cr
#'     Events are excluded if either of the upstream or downstream
#'       event is lower than total splicing events by a log-2 magnitude 
#'       above filterVars$maximum. For example, if 
#'       \code{filterVars$maximum = 2}, we require both upstream and downstream
#'       events to represent at least 1/(2^2) = 1/4 of the sum of upstream
#'       and downstream event.
#'       This is considered for each isoform of each event, as long as the
#'       total counts belonging to the considered isoform is above 
#'       \code{filterVars$minDepth}.
#'
#'   We highly recommend using the default filters, which can be acquired 
#'     using \code{\link{get_default_filters}}
#' 
#' @param filterClass One of either \code{Annotation} or \code{Data}
#' @param filterType For filterClass \code{Annotation}, either one of 
#'   \code{Protein_Coding}, \code{NMD_Switching}, 
#'   \code{Transcript_Support_Level}.
#'   For filterClass \code{Data}, either one of
#'   \code{Depth}, \code{Coverage}, \code{Consistency}.
#' @param filterVars A list of parameters, as explained below
#' @param filterObject the SummarizedExperiment to filter
#'
#' @return A vector of type \code{logical} designating which events to retain 
#'   \code{TRUE} and which to remove \code{FALSE}.
#' @examples
#' # see ?MakeSE on example code of generating this NxtSE object
#' se = NxtIRF_example_NxtSE()
#' 
#' # Get NxtIRF recommended filters
#' filters = get_default_filters()
#' 
#' # Filter the SummarizedExperiment using the first default filter ("Depth")
#' se.depthfilter = se[runFilter(
#'         filterClass = filters[[1]]$filterClass,
#'         filterType = filters[[1]]$filterType,
#'         filterVars = filters[[1]]$filterVars,
#'         filterObject = se
#'     ), ]
#' @seealso \code{\link{get_default_filters}}, 
#' \code{\link{apply_filters}}
#' @export
runFilter <- function(filterClass, filterType, filterVars, filterObject) {
    # filterClass: can be one of 'Annotation', 'Data', 'Runtime'
    # filterType:
    # - Annotation:
    # - Data:
    # -     Depth: 1-minimum, 2-minCond, 3-pcTRUE
    # -     Coverage: 1-minimum, 1a-minDepth, 2-minCond, 3-pcTRUE
    # -     Consistency: 1-maximum, 1a-minDepth, 2-minCond, 3-pcTRUE
    # - for Consistency, maximum is the max(abs(log2_delta)) between comparison
    #       and calculated value
    if(filterClass == "Data") {
        if(filterType == "Depth") {
            message("Running Depth filter")
            return(.runFilter_data_depth(filterObject, filterVars))
        } else if(filterType == "Coverage") {
            message("Running Coverage filter")
            return(.runFilter_data_coverage(filterObject, filterVars))
        } else if(filterType == "Consistency") {    # requires: 
            message("Running Consistency filter")
            return(.runFilter_data_consistency(filterObject, filterVars))
        }
    } else if(filterClass == "Annotation") {
        if(filterType == "Protein_Coding") {
            # returns if any of included or excluded is protein_coding
            rowSelected = as.data.table(rowData)
            rowSelected = rowSelected[
                get("Inc_Is_Protein_Coding") == TRUE | 
                get("Exc_Is_Protein_Coding") == TRUE]
            rowSelected = rowSelected[get("EventType") != "IR" | 
                get("Inc_Is_Protein_Coding") == TRUE] # filter for CDS introns
            res = rowData$EventName %in% rowSelected$EventName
        } else if(filterType == "NMD_Switching") {
            rowSelected = as.data.table(rowData)
            rowSelected = rowSelected[!is.na(get("Inc_Is_NMD")) & 
                !is.na(get("Exc_Is_NMD"))]
            rowSelected = rowSelected[get("Inc_Is_NMD") != get("Exc_Is_NMD")]
            res = rowData$EventName %in% rowSelected$EventName
        } else if(filterType == "Transcript_Support_Level") {
            if(!("minimum" %in% names(filterVars))) {
                minimum = 1
            } else {
                minimum = filterVars$minimum
            }
            rowSelected = as.data.table(rowData)
            rowSelected = rowSelected[get("Inc_TSL") != "NA" & 
                get("Exc_TSL") != "NA"]
            rowSelected[, c("Inc_TSL") := as.numeric(get("Inc_TSL"))]
            rowSelected[, c("Exc_TSL") := as.numeric(get("Exc_TSL"))]
            rowSelected = rowSelected[get("Inc_TSL") <= minimum & 
                get("Exc_TSL") <= minimum]
            res = rowData$EventName %in% rowSelected$EventName
        }
        if("EventTypes" %in% names(filterVars)) {
            res[!(rowData$EventType %in% filterVars$EventTypes)] = TRUE
        }
        return(res)
    } else {
        return(rep(TRUE, nrow(filterObject)))
    }
}

################################################################################
# Individual functions:

.runFilter_data_depth <- function(se, filterVars) {
    usePC = ifelse("pcTRUE" %in% names(filterVars),
        filterVars$pcTRUE, 100)
    use_cond = ifelse(
            !is.null(names(filterVars)) && 
            "condition" %in% names(filterVars) && 
            filterVars$condition %in% colnames(colData), 
        TRUE, FALSE)

    colData = as.data.frame(colData(se))
    rowData = as.data.frame(rowData(se))

    minimum = ifelse("minimum" %in% names(filterVars),
        filterVars$minimum, 20)
        
    if(use_cond == TRUE) {
        cond_vec = unlist(
            colData[, which(colnames(colData) == filterVars$condition)])
        cond_vars = unique(cond_vec)
    }
    depth = as.matrix(assay(se, "Depth"))
    sum_res = rep(0, nrow(se))
    if(use_cond == TRUE) {
        for(cond in cond_vars) {
            depth.subset = depth[, which(cond_vec == cond)]
            sum = rowSums(depth.subset >= minimum)
            sum_res = sum_res + 
                ifelse(sum * 100 / ncol(depth.subset) >= usePC, 1, 0)
        }
        n_TRUE = ifelse(
            !is.na(suppressWarnings(as.numeric(filterVars$minCond))),
            as.numeric(filterVars$minCond), -1
        )
        if(n_TRUE == -1) n_TRUE = length(cond_vars)
        res = (sum_res >= n_TRUE)
    } else {
        sum = rowSums(depth >= minimum)
        res = ifelse(sum * 100 / ncol(depth) >= usePC, TRUE, FALSE)
    }
    if("EventTypes" %in% names(filterVars)) {
        res[!(rowData$EventType %in% filterVars$EventTypes)] = TRUE
    }
    return(res)
}

.runFilter_data_coverage <- function(se, filterVars) {
    usePC = ifelse("pcTRUE" %in% names(filterVars),
        filterVars$pcTRUE, 100)
    use_cond = ifelse(
            !is.null(names(filterVars)) && 
            "condition" %in% names(filterVars) && 
            filterVars$condition %in% colnames(colData), 
        TRUE, FALSE)
    minDepth = ifelse("minDepth" %in% names(filterVars),
        filterVars$minDepth, 0)

    colData = as.data.frame(colData(se))
    rowData = as.data.frame(rowData(se))

    minimum = ifelse("minimum" %in% names(filterVars),
        filterVars$minimum, 20)
    if(use_cond == TRUE) {
        cond_vec = unlist(
            colData[, which(colnames(colData) == filterVars$condition)])
        cond_vars = unique(cond_vec)
    }
    cov = as.matrix(assay(se, "Coverage"))
    depth = as.matrix(assay(se, "minDepth"))
    
    # do not test if depth below threshold
    cov[depth < minDepth] = 1    

    sum_res = rep(0, nrow(se))
    if(use_cond == TRUE) {
        for(cond in cond_vars) {
            cov.subset = cov[, which(cond_vec == cond)]
            sum = rowSums(cov.subset >= minimum / 100)
            sum_res = sum_res + 
                ifelse(sum * 100 / ncol(cov.subset) >= usePC, 1, 0)
        }
        n_TRUE = ifelse(
            !is.na(suppressWarnings(as.numeric(filterVars$minCond))), 
            as.numeric(filterVars$minCond), -1)
        if(n_TRUE == -1) n_TRUE = length(cond_vars)
        res = (sum_res >= n_TRUE)
    } else {
        sum = rowSums(cov >= minimum / 100)
        res = ifelse(sum * 100 / ncol(cov) >= usePC, TRUE, FALSE)
    }
    if("EventTypes" %in% names(filterVars)) {
        res[!(rowData$EventType %in% filterVars$EventTypes)] = TRUE
    }
    return(res)
}

.runFilter_data_consistency <- function(se, filterVars) {
    usePC = ifelse("pcTRUE" %in% names(filterVars),
        filterVars$pcTRUE, 100)
    use_cond = ifelse(
            !is.null(names(filterVars)) && 
            "condition" %in% names(filterVars) && 
            filterVars$condition %in% colnames(colData), 
        TRUE, FALSE)
    minDepth = ifelse("minDepth" %in% names(filterVars),
        filterVars$minDepth, 0)

    colData = as.data.frame(colData(se))
    rowData = as.data.frame(rowData(se))

    maximum = ifelse("maximum" %in% names(filterVars),
        filterVars$maximum, 1)

    if(use_cond == TRUE) {
        cond_vec = unlist(colData[, 
            which(colnames(colData) == filterVars$condition)])
        cond_vars = unique(cond_vec)
    }
    Up_Inc = as.matrix(S4Vectors::metadata(se)$Up_Inc)[
        rowData(se)$EventName[
            rowData(se)$EventType %in% c("IR", "MXE", "SE", "RI")
        ],
    ]
    Down_Inc = as.matrix(S4Vectors::metadata(se)$Down_Inc)[
        rowData(se)$EventName[
            rowData(se)$EventType %in% c("IR", "MXE", "SE", "RI")
        ],
    ]
    IntronDepth = as.matrix(assay(se, "Included"))
    IntronDepth = IntronDepth[
        rowData$EventType %in% c("IR", "MXE", "SE", "RI"),]
    minDepth.Inc = Up_Inc + Down_Inc
    # do not test if depth below threshold
    Up_Inc[minDepth.Inc < minDepth] = IntronDepth[minDepth.Inc < minDepth]
    Down_Inc[minDepth.Inc < minDepth] = IntronDepth[minDepth.Inc < minDepth]

    Excluded = as.matrix(assay(se, "Excluded"))
    Excluded = Excluded[rowData$EventType %in% c("MXE"),]
    Up_Exc = as.matrix(S4Vectors::metadata(se)$Up_Exc)
    Down_Exc = as.matrix(S4Vectors::metadata(se)$Down_Exc)
    minDepth.Exc = Up_Exc + Down_Exc
    # do not test if depth below threshold
    Up_Exc[minDepth.Exc < minDepth] = Excluded[minDepth.Exc < minDepth]    
    Down_Exc[minDepth.Exc < minDepth] = Excluded[minDepth.Exc < minDepth]

    sum_res = rep(0, nrow(se))
    if(use_cond == TRUE) {
        for(cond in cond_vars) {
            Up_Inc.subset = Up_Inc[, which(cond_vec == cond)]
            Down_Inc.subset = Down_Inc[, which(cond_vec == cond)]
            IntronDepth.subset = IntronDepth[, which(cond_vec == cond)]
            Up_Exc.subset = Up_Exc[, which(cond_vec == cond)]
            Down_Exc.subset = Down_Exc[, which(cond_vec == cond)]
            Excluded.subset = Excluded[, which(cond_vec == cond)]

            sum_inc = rowSums(
                abs(log2(Up_Inc.subset + 1) - log2(IntronDepth.subset + 1)) 
                    < maximum &
                abs(log2(Down_Inc.subset + 1) - log2(IntronDepth.subset + 1)) 
                    < maximum
            )
            sum_exc = rowSums(
                abs(log2(Up_Exc.subset + 1) - log2(Excluded.subset + 1)) 
                    < maximum &
                abs(log2(Down_Exc.subset + 1) - log2(Excluded.subset + 1)) 
                    < maximum
            )
            sum_inc = c(sum_inc, rep(ncol(Up_Inc.subset), 
                sum(!(rowData$EventType %in% c("IR", "MXE", "SE", "RI")))))
            sum_exc = c(rep(ncol(Up_Inc.subset), 
                sum(rowData$EventType == "IR")),
            sum_exc, rep(ncol(Up_Inc.subset), 
                sum(!(rowData$EventType %in% c("IR", "MXE")))))
            sum = 0.5 * (sum_inc + sum_exc)
            sum_res = sum_res + 
                ifelse(sum * 100 / ncol(Up_Inc.subset) >= usePC, 1, 0)
        }
        n_TRUE = ifelse(
            !is.na(suppressWarnings(as.numeric(filterVars$minCond))), 
            as.numeric(filterVars$minCond), -1)
        if(n_TRUE == -1) n_TRUE = length(cond_vars)
        res = (sum_res >= n_TRUE)
    } else {
        sum_inc = rowSums(
            abs(log2(Up_Inc + 1) - log2(IntronDepth +1)) < filterVars$maximum &
            abs(log2(Down_Inc + 1) - log2(IntronDepth +1)) < filterVars$maximum
        )
        sum_exc = rowSums(
            abs(log2(Up_Exc + 1) - log2(Excluded +1)) < filterVars$maximum &
            abs(log2(Down_Exc + 1) - log2(Excluded +1)) < filterVars$maximum
        )
        sum_inc = c(sum_inc, rep(ncol(Up_Inc), 
            sum(!(rowData$EventType %in% c("IR", "MXE", "SE", "RI")))))
        sum_exc = c(
            rep(ncol(Up_Inc), sum(rowData$EventType == "IR")),
            sum_exc, 
            rep(
                ncol(Up_Inc), 
                sum(!(rowData$EventType %in% c("IR", "MXE")))
            )
        )
        sum = 0.5 * (sum_inc + sum_exc)
        res = ifelse(sum * 100 / ncol(Up_Inc) >= usePC, TRUE, FALSE)
    }
    if("EventTypes" %in% names(filterVars)) {
        res[!(rowData(se)$EventType %in% filterVars$EventTypes)] = TRUE
    }
    return(res)
}

################################################################################



#' Convenience function to apply a list of filters to a SummarizedExperiment
#'   object
#'
#' See [runFilter()] for details regarding filters
#' 
#' @param se A SummarizedExperiment object created by `MakeSE()`
#' @param filters A list of filters to apply. Each filter must contain the
#'   elements `filterClass`, `filterType` and `filterVars`. 
#'   See `?runFilter` for details
#' @return A vector of logicals, with `TRUE` indicating events to be retained,
#'   and `FALSE` for events to be filtered out
#' @examples
#' # see ?MakeSE on example code of generating this NxtSE object
#' se = NxtIRF_example_NxtSE()
#' 
#' # Get NxtIRF recommended filters
#' filters = get_default_filters()
#' 
#' # Apply all recommended filters:
#' se.filtered = apply_filters(
#'     se = se, 
#'     filters = filters
#' )
#' @seealso [get_default_filters()], [runFilter()]
#' @md
#' @export
apply_filters <- function(se, filters = get_default_filters()) {
    # filters are a list of filters to apply on se
    # returns a vector of TRUE / FALSE
    # a filtered se can be made using:
    #       se.filtered = se[apply_filters(se, filters),]
    if(!is(filters, "list")) {
        .log(paste("In apply_filters(),",
            "filters must be a list"))
    }
    for(i in seq_len(length(filters))) {
        if(!("filterVars" %in% names(filters[[i]]))) {
            .log(paste("In apply_filters(),",
                "filterVars is missing from filters @ index #", i))
        }
        if(!("filterClass" %in% names(filters[[i]]))) {
            .log(paste("In apply_filters(),",
                "filterClass is missing from filters @ index #", i))
        }
        if(!("filterType" %in% names(filters[[i]]))) {
            .log(paste("In apply_filters(),",
                "filterType is missing from filters @ index #", i))
        }
    }
    if(!is(se, "NxtSE")) {
        .log(paste("In apply_filters(),",
            "se must be a NxtSE object"))
    }
    filterSummary = rep(TRUE, nrow(se))
    for(i in seq_len(length(filters))) {
        filterSummary = filterSummary & runFilter(
            filters[[i]]$filterClass,
            filters[[i]]$filterType,
            filters[[i]]$filterVars,
            se
        )
    }
    
    return(filterSummary)
}