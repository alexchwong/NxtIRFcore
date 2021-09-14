################################################################################

#' Constructs a NxtSE object from the collated data
#'
#' MakeSE() creates a NxtSE object from the data collated
#' from IRFinder output using [CollateData()].
#'
#' @details NxtSE is a specialised class which inherits from the
#'   \linkS4class{SummarizedExperiment} class, NxtSE is used in NxtIRF
#'   for downstream analysis. It contains processed data collated from IRFinder
#'   output from multiple samples. Counts of included or excluded
#'   IR and alternative splicing events, as well as relevant QC parameters
#'   are included. See [NxtSE-methods] for further details.
#'
#' @param collate_path The output path given to CollateData() pointing to the
#'   collated data
#' @param colData A data frame containing the sample annotation information.
#'   Note that the first column must contain the sample names. If the names of 
#'   only a subset of samples are given, then `MakeSE()` will construct the SE 
#'   object based only on the samples given. Omit `colData` to generate an SE 
#'   object based on the whole dataset. The colData can be set later using 
#'   `colData()` for `SummarizedExperiment` objects
#' @param RemoveOverlapping (default = TRUE) Whether to filter out overlapping 
#'   introns of IR events belonging to minor isoforms. MakeSE will try to 
#'   identify which junctions belong to major isoforms, then select the 
#'   junctions from non-overlapping minor isoforms in an iterative approach, 
#'   until no non-overlapping introns remain. This is important
#'   to make sure IR events are not 'double-counted'
#'
#' @return A NxtIRF SummarizedExperiment (`NxtSE`) object 
#'
#' @examples
#' se = MakeSE(collate_path = file.path(tempdir(), "NxtIRF_output"))
#' @md
#' @export
MakeSE = function(collate_path, colData, RemoveOverlapping = TRUE) {
    # Includes iterative filtering for IR events with highest mean PSI
        # To annotate IR events of major isoforms

    colData <- .makeSE_validate_args(collate_path, colData)
    colData <- .makeSE_colData_clean(colData)

    N <- 8
    dash_progress("Loading NxtSE object from file...", N)
    .log("Loading NxtSE object from file...", "message", appendLF = FALSE)

    se = .MakeSE_load_NxtSE(file.path(collate_path, "NxtSE.rds"))
    
    # Encapsulate as NxtSE object
    se = as(se, "NxtSE")

    # Subset
    se = se[, colData$sample]
    if(ncol(colData) > 1) {
        colData_use <- colData[, -1, drop = FALSE]
        rownames(colData_use) <- colData$sample
        colData(se) <- as(colData_use, "DataFrame")    
    }
    
    message("done\n")
    
    if(RemoveOverlapping == TRUE) {
        dash_progress("Removing overlapping introns...", N)
        se <- .makeSE_iterate_IR(se, collate_path)
    }
    return(se)
}



################################################################################

# Checks:
# - whether the given path contains a valid CollateData() output
# - whether 
.makeSE_validate_args <- function(collate_path, colData) {
    item.todo = c("rowEvent", "Included", "Excluded", "Depth", "Coverage", 
        "minDepth", "Up_Inc", "Down_Inc", "Up_Exc", "Down_Exc")

    if(!file.exists(file.path(collate_path, "colData.Rds"))) {
        .log(paste("In MakeSE():",
            file.path(collate_path, "colData.Rds"),
            "was not found"))
    }
    colData.Rds = readRDS(file.path(collate_path, "colData.Rds"))
    if(!("df.anno" %in% names(colData.Rds))) {
        .log(paste("In MakeSE():",
            file.path(collate_path, "colData.Rds"),
            "must contain df.anno containing annotations"))
    }
    if(missing(colData)) {    
        colData = colData.Rds$df.anno
    } else {
        colData = as.data.frame(colData)
        if(!("sample" %in% colnames(colData))) {
            colnames(colData)[1] = "sample"
        }
        if(!all(colData$sample %in% colData.Rds$df.anno$sample)) {
            .log(paste("In MakeSE():",
                "some samples in colData were not found in given path"),
                "message")
            colData = colData[colData$sample %in% colData.Rds$df.anno$sample,]
        }
    }
    return(colData)
}

# Converts charactor vectors to factors, removes columns with all NA's
.makeSE_colData_clean <- function(colData) {
    remove_na = NULL
    if(ncol(colData) > 1) {
        for(i in seq(2, ncol(colData))) {
            if(is(colData[,i], "character")) {
                colData[,i] = factor(unlist(colData[,i]))      
            } else if(is(colData[,i], "logical")) {
                colData[,i] <- factor(unlist(
                    ifelse(colData[,i], "TRUE","FALSE")))                
            } else if(all(is.na(unlist(colData[,i])))) {
                remove_na = append(remove_na, i)
            }
        }
    }
    if(!is.null(remove_na)) {
        colData = colData[,-remove_na]
    }
    return(colData)
}

# Loads a NxtSE RDS
.MakeSE_load_NxtSE <- function(filepath) {
    se <- readRDS(filepath)
    path <- dirname(filepath)
    se@assays <- .collateData_expand_assay_paths(se@assays, path)
    se@metadata[["Up_Inc"]] <- .collateData_expand_assay_path(
        se@metadata[["Up_Inc"]], path)
    se@metadata[["Down_Inc"]] <- .collateData_expand_assay_path(
        se@metadata[["Down_Inc"]], path)
    se@metadata[["Up_Exc"]] <- .collateData_expand_assay_path(
        se@metadata[["Up_Exc"]], path)
    se@metadata[["Down_Exc"]] <- .collateData_expand_assay_path(
        se@metadata[["Down_Exc"]], path)
    return(se)
}

# Iterates through IRFinder introns; removes overlapping minor introns
.makeSE_iterate_IR <- function(se, collate_path) {

    junc_PSI <- HDF5Array(file.path(normalizePath(collate_path), 
        "data.h5"), "junc_PSI")[, colnames(se), drop = FALSE]

    se.IR = se[rowData(se)$EventType == "IR",,drop = FALSE]
    se.coords = rowData(se.IR)$EventRegion[
        rowData(se.IR)$EventRegion %in% rownames(junc_PSI)]
    
    if(length(se.coords) > 0) {
        .log(paste("Iterating through IR events to determine introns",
            "of main isoforms"), type = "message")
        include <- .makeSE_iterate_IR_select_events(se.coords, junc_PSI)
        se.coords.final = se.coords[include]
        se.coords.excluded = se.coords[!include]

        # Iteration to find events not overlapping with se.IR.final
        include <- .makeSE_iterate_IR_retrieve_excluded_introns(
            se.coords.final, se.coords.excluded)
        iteration = 0
        while(length(include) > 0 & length(se.coords.final) > 0) {
            iteration = iteration + 1
            .log(paste("Iteration", iteration), type = "message")
            dash_progress(paste("Iteration", iteration), 8)
            se.coords.excluded = se.coords.excluded[include]

            include <- .makeSE_iterate_IR_select_events(
                    se.coords.excluded, junc_PSI)

            if(length(include) > 0) {
                se.coords.final = c(se.coords.final, 
                    se.coords.excluded[include])
                se.coords.excluded = 
                    se.coords.excluded[!include]
                include <- .makeSE_iterate_IR_retrieve_excluded_introns(
                    se.coords.final, se.coords.excluded)
            } else {
                include = c()
            }
        }

        se = se[c(
            which(rowData(se.IR)$EventRegion %in% se.coords.final), 
            which(rowData(se)$EventType != "IR")
        ),]
    }
    return(se)
}

# Selects introns of major isoforms
.makeSE_iterate_IR_select_events <- function(se.coords, junc_PSI) {
    gr = NxtIRF.CoordToGR(se.coords)
    gr.reduced = reduce(gr)

    OL = findOverlaps(gr, gr.reduced)
    junc_PSI.group = as.data.table(junc_PSI[se.coords,, drop = FALSE])
    junc_PSI.group$means = rowMeans(junc_PSI.group)
    junc_PSI.group$group = to(OL)
    junc_PSI.group[, c("max_means") := max(get("means")), 
        by = "group"]
    return(junc_PSI.group$means == junc_PSI.group$max_means)
}

# Find excluded introns that does not overlap with given selection of introns
.makeSE_iterate_IR_retrieve_excluded_introns <- function(
        se.coords.final, se.coords.excluded) {
    if(length(se.coords.excluded) > 0) {
        final.gr = NxtIRF.CoordToGR(se.coords.final)
        excluded.gr = NxtIRF.CoordToGR(se.coords.excluded)

        OL = findOverlaps(excluded.gr, final.gr)
        include = which(!(
            seq_len(length(excluded.gr))) %in% sort(unique(from(OL))))
    } else {
        include = c()
    }
    return(include)
}
