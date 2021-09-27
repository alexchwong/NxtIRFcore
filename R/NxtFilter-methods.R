setMethod("initialize", "NxtFilter", function(.Object,
        filterClass = c("Data", "Annotation"),
        filterType = c(
            "Depth", "Coverage", "Consistency",
            "Protein_Coding", "NMD", "TSL"
        ),
        pcTRUE = 100, minimum = 20, maximum = 1, minDepth = 5,
        condition = "", minCond = -1,
        EventTypes = c("IR", "MXE", "SE", "A3SS", "A5SS", "AFE", "ALE", "RI")
) {
    .Object <- callNextMethod()

    filterClass = match.arg(filterClass)
    filterType = match.arg(filterType)

    if(filterClass == "") 
        .log("filterClass must be one of `Data` or `Annotation`")
    if(filterClass == "Data") {
        if(!(filterType %in% c("Depth", "Coverage", "Consistency")))
        .log("filterClass must be one of `Depth`, `Coverage` or `Consistency`")
    } else {
        if(!(filterType %in% c("Protein_Coding", "NMD", "TSL")))
        .log("filterClass must be one of `Protein_Coding`, `NMD` or `TSL`")    
    }
        

    pcTRUE = as.numeric(pcTRUE)
    pcTRUE = min(100, max(0, pcTRUE))
    minimum = as.numeric(minimum)
    minimum = max(0, minimum)
    maximum = as.numeric(maximum)
    maximum = max(0, maximum)
    minDepth = as.numeric(minDepth)
    minDepth = max(0, minDepth)
    minCond = as.numeric(minCond)
    
    # Only 1 condition allowed
    if(length(condition) != 1) .log("`condition` should be a character scalar")
    
    et_args = c("IR", "MXE", "SE", "A3SS", "A5SS", "AFE", "ALE", "RI")
    EventTypes = EventTypes[EventTypes %in% et_args]
    if(length(EventTypes) == 0) EventTypes = et_args

    .Object@filterClass = filterClass
    .Object@filterType = filterType
    .Object@pcTRUE = pcTRUE
    .Object@minimum = minimum
    .Object@maximum = maximum
    .Object@minDepth = minDepth
    .Object@condition = condition
    .Object@minCond = minCond
    .Object@EventTypes = EventTypes
    
    .Object
})

setMethod("show", "NxtFilter", function(object) {
    .cat_NxtFilter_common(object)
    .cat_NxtFilter_filterSpecific(object)
})

# Describe class, type, and conditions filter will be applied across
.cat_NxtFilter_common <- function(object) {
    .nxtcat(paste0("NxtFilter Class: ", .colourise("%s\t", "red")), 
        object@filterClass)
    .nxtcat(paste0("Type: ", .colourise("%s\n", "purple")), 
        object@filterType)
    .nxtcat(paste0("EventTypes: ", .colourise("%s\n", "green")), 
        paste(object@EventTypes, collapse = " "))
    if(object@filterClass == "Data") {
        cat("Filter must be passed in at least ")
        if(object@condition == "") {
            .nxtcat(paste0(.colourise("%.1f", "purple"), 
                " percent of all samples\n"), object@pcTRUE)
        } else {
            if(object@minCond == -1) {
            .nxtcat(paste0(.colourise("%.1f", "purple"), 
                " percent of all categories of ",
                .colourise("%s\n", "purple")), object@pcTRUE, object@condition)
            } else {
            .nxtcat(paste0(.colourise("%.1f", "purple"), 
                " percent of ", .colourise("%i","green"), 
                " categories of ", .colourise("%s\n", "purple")), 
                object@pcTRUE, as.integer(object@minCond), object@condition)
            }
        }
    }
}

# Describe filter-specific functions
.cat_NxtFilter_filterSpecific <- function(object) {
    if(object@filterType == "Depth") {
        .nxtcat(paste0("Minimum Event Depth: ", 
            .colourise("%i\n", "red")), as.integer(object@minimum))
        .cat_filter_info("minDepth")
    } else if(object@filterType == "Coverage") {
        .nxtcat(paste0("Minimum Coverage Percentage: ", 
            .colourise("%.1f\n", "red")), object@minimum)
        .cat_filter_info("Coverage minimum", object@EventTypes)
        .nxtcat(paste0("Event Depth below ", .colourise("%i", "purple"), 
            " are ignored\n"), as.integer(object@minDepth))
        .cat_filter_info("minDepth")
    } else if(object@filterType == "Consistency") {
        .nxtcat(paste0("Sub-event minimum proportion: ",
            .colourise("%.1f", "red"), " of total isoform splice count\n"), 
            2^(-object@maximum))
        .cat_filter_info("Consistency maximum", object@EventTypes)
        .nxtcat(paste0("Event Depth below ", .colourise("%i", "purple"), 
            " are ignored\n"), as.integer(object@minDepth))
        .cat_filter_info("minDepth")
    } else if(object@filterType == "Protein_Coding") {
        cat("Events of which neither isoform encodes protein are removed")
    } else if(object@filterType == "NMD") {
        cat("Events of which neither isoform are NMD substrates are removed")
    } else if(object@filterType == "TSL") {
        .nxtcat(paste0("Transcript Support Level: ", 
            .colourise("%i", "purple"), " or higher\n"), 
            as.integer(object@minimum))
        cat("Events with both isoforms belonging to ")
        cat("lower-ranking TSLs are removed")
    }
}

# Describe more info re specific functions
.cat_filter_info <- function(mode, EventTypes = "") {
    if(mode == "minDepth") {
        cat("Event Depth refers to number of aligned splice reads plus ")
        cat("effective depth of their introns\n")    
    } else if(mode == "Coverage minimum") {
        if(any(EventTypes %in% c("IR", "RI"))) {
            cat("In retained introns, coverage refers to the proportion of the")
            cat(" measured intron that is covered by at least 1 alignment\n")    
        }
        if(any(EventTypes %in% c("MXE", "SE", "ALE", "AFE", "A3SS", "A5SS"))) {
            cat("In splice events, coverage refers to the proportion of ")
            cat("junction reads that belong to either the included or excluded")
            cat(" isoforms of the given event\n")    
        }
    } else if(mode == "Consistency maximum") {
        if(any(EventTypes %in% c("IR", "RI"))) {
            cat("In retained introns, sub-events refer to 5'- and 3'-")
            cat("exon-intron spanning alignments\n")    
        }
        if(any(EventTypes %in% c("MXE", "SE"))) {
            cat("In mutually exclusive exons and skipped exons, ")
            cat("sub-events refer to the upstream and downstream splice ")
            cat("alignments\n")
        }
    }
}

#' NxtIRF filters for quality alternative splicing and intron retention events
#'
#' @param filterClass Must be either `"Data"` or `"Annotation"`. See details
#' @param filterType If `filterClass = "Data"`, then must be one of 
#'   `c("Depth", "Coverage", "Consistency")`. If `filterClass = "Annotation"`,
#'   must be one of `c("Protein_Coding", "NMD", "TSL")`. See details
#' @param pcTRUE If conditions are set, what percentage of all samples in each
#'   of the condition must the filter be satisfied for the event to pass the
#'   filter check. Must be between 0 and 100 (default 100)
#' @param minimum Filter-dependent argument. See details
#' @param maximum Filter-dependent argument. See details
#' @param minDepth Filter-dependent argument. See details
#' @param condition (default "") If set, must match the name of an experimental
#'   condition in the NxtSE object to be filtered, 
#'   i.e. a column name in `colData(se)`. Leave blank to disable filtering
#'   by condition
#' @param minCond (default -1) If condition is set, how many minimum number of
#'   conditions must pass criteria must the filters be passed. For example,
#'   if condition = "Batch", and batches are "A", "B", or "C", setting
#'   `minCond = 2` with `pcTRUE = 100` means that all samples within two of the
#'   three types of `Batch` must pass the filter criteria for the ASE event
#'   to pass the filter. Setting `-1` means all elements of `condition` must
#'   pass criteria. Set to `-1` when the number of elements in the experimental
#'   condition is unknown.
#' @param EventTypes What types of events are considered for filtering. Must be 
#'   one of `c("IR", "MXE", "SE", "A3SS", "A5SS", "AFE", "ALE", "RI")`. Events 
#'   not specified in `EventTypes` are not filtered (i.e. they will pass the
#'   filter without checks)
#' @details
#'   **Annotation Filters**\cr\cr 
#'     - **Protein_Coding**: Filters for alternative splicing or IR events 
#'       within protein reading frames. No additional parameters required.\cr\cr
#'     - **NMD**: Filters for events in which one isoform is a
#'       predicted NMD substrate.\cr\cr 
#'     - **TSL**: filters for events in which both
#'       isoforms have a TSL level below or equal to `minimum`\cr\cr 
#'   **Data Filters**\cr\cr 
#'     - **Depth**: Filters IR or alternative splicing events of transcripts
#'       that are "expressed" with adequate `Depth` as calculated by the
#'       sum of all splicing and IR reads spanning the event. Events with 
#'       `Depth` below `minimum` are filtered out\cr\cr
#'     - **Coverage**: Coverage means different things to IR and alternative
#'       splicing.\cr\cr 
#'     For *IR*, Coverage refers to the percentage of the measured intron
#'       covered with reads. Introns of samples with an IntronDepth above 
#'       `minDepth` are assessed, with introns with coverage 
#'       below `minimum` are filtered out.\cr\cr 
#'     For *Alternative Splicing*, Coverage refers to the percentage of all
#'       splicing events observed across the genomic region that is compatible
#'       with either the included or excluded event. This prevents NxtIRF from 
#'       doing differential analysis between two minor isoforms. Instead of
#'       IntronDepth, in AS events NxtIRF considers events where the spliced
#'       reads from both exonic regions exceed `minDepth`.
#'       Then, events with a splicing coverage below `minimum`
#'       are excluded. We recommend testing IR events for > 90% coverage and AS
#'       events for > 60% coverage as given in the default filters which can be
#'       accessed using [get_default_filters]\cr\cr  
#'     - **Consistency**: Skipped exons (SE) and mutually exclusive exons
#'       (MXE) comprise reads of two contiguous splice junctions (for the
#'       included casette exon). Summating counts from both junctions is
#'       misleading as there may be overlapping events (e.g. alternate first 
#'       / last exons) that only rely on one splice event. To ensure the SE /
#'       MXE is the dominant event, we require both splice junctions to have
#'       comparable counts.\cr\cr
#'     Events are excluded if either of the upstream or downstream
#'       event is lower than total splicing events by a log-2 magnitude 
#'       above `maximum`. For example, if 
#'       `maximum = 2`, we require both upstream and downstream
#'       events to represent at least 1/(2^2) = 1/4 of the sum of upstream
#'       and downstream event.
#'       This is considered for each isoform of each event, as long as the
#'       total counts belonging to the considered isoform is above 
#'       `minDepth`.
#'
#'   We highly recommend using the default filters, which can be acquired 
#'     using [get_default_filters]
#' @return A NxtFilter object with the specified parameters
#' @examples
#' # Create a NxtFilter that filters for protein-coding ASE
#' f1 = NxtFilter(filterClass = "Annotation", filterType = "Protein_Coding")
#'
#' # Create a NxtFilter that filters for Depth >= 20 in IR events
#' f2 = NxtFilter(
#'     filterClass = "Data", filterType = "Depth",
#'     minimum = 20, EventTypes = c("IR", "RI")
#' )
#'
#' # Create a NxtFilter that filters for Coverage > 60% in splice events
#' # that must be satisfied in at least 2 categories of condition "Genotype"
#' f3 = NxtFilter(
#'     filterClass = "Data", filterType = "Coverage",
#'     minimum = 60, EventTypes = c("MXE", "SE", "AFE", "ALE", "A3SS", "A5SS"),
#'     condition = "Genotype", minCond = 2
#' )
#'
#' # Create a NxtFilter that filters for Depth > 10 in all events
#' # that must be satisfied in at least 50% of each gender
#' f4 = NxtFilter(
#'     filterClass = "Data", filterType = "Depth",
#'     minimum = 10, condition = "gender", pcTRUE = 50
#' )
#'
#' # Get a description of what these filters do:
#' f1
#' f2
#' f3
#'
#' @seealso [Run_NxtIRF_Filters]
#' @md
#' @export
NxtFilter <- function(
        filterClass = c("Data", "Annotation"),
        filterType = c(
            "Depth", "Coverage", "Consistency",
            "Protein_Coding", "NMD", "TSL"
        ),
        pcTRUE = 100, minimum = 20, maximum = 1, minDepth = 5,
        condition = "", minCond = -1,
        EventTypes = c("IR", "MXE", "SE", "A3SS", "A5SS", "AFE", "ALE", "RI")
) {
    nfv = new("NxtFilter",
        filterClass = filterClass, filterType = filterType, pcTRUE = pcTRUE, 
        minimum = minimum, maximum = maximum, minDepth = minDepth,
        condition = condition, minCond = minCond,
        EventTypes = EventTypes
    )
    nfv
}