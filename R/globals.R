globalVariables(c(":=","."))

is.nan.data.frame <- function(x) do.call(cbind, lapply(x, is.nan))

is_valid <- function(x) {
    !missing(x) &&
    !is.null(x) && length(x) > 0 && 
        (isS4(x) || !is.na(x)) && 
        (!is.character(x) || (x != "" && x != "(none)"))
}

#' Converts an IGV-style coordinate to a GenomicRanges object
#'
#' IGV-style coordinates typically have the syntax 
#'   `seqnames:start-end/strand`\cr\cr
#' For example: "chr3:10550-10730/+" or "X:51231-51330/-"
#' @param coordinates A vector of strings containing the coordinates
#'   to be converted
#' @return A GRanges object that corresponds to the given coordinates
#' @examples
#' se = NxtIRF_example_NxtSE()
#'
#' coordinates = rowData(se)$EventRegion
#'
#' coords.gr = NxtIRF.CoordToGR(coordinates)
#' @md
#' @export
NxtIRF.CoordToGR = function(coordinates) {
    temp = tstrsplit(coordinates,split="/")
    strand = as.character(temp[[2]])
    temp2 = tstrsplit(temp[[1]],split=":")
    seqnames = temp2[[1]]
    temp3 = tstrsplit(temp2[[2]],split="-")
    start = temp3[[1]]
    end = temp3[[2]]
    return(GRanges(seqnames = seqnames, ranges = IRanges(
        start = as.numeric(start), end = as.numeric(end)),
        strand = strand))
}

NxtIRF.CheckPackageInstalled <- function(
        package = "DESeq2", version = "1.0.0", 
        returntype = c("error", "warning", "silent")) {
    res = tryCatch(
        ifelse(packageVersion(package)>=version, TRUE, FALSE),
        error = function(e) FALSE)
    if(!res) {
        returntype = match.arg(returntype)
        if(returntype == "error") {
            stop(paste(
                package, "version", version, "is not installed;",
                "and is required for this function"
            ), call. = FALSE)
        } else if(returntype == "warning") {
            warning(paste(
                package, "version", version, "is not installed;",
                "and is required for this function"
            ))
            return(FALSE)
        } else {
            return(FALSE)
        }
    } else {
        return(TRUE)
    }
}

.validate_threads <- function(n_threads) {
    n_threads_to_use = as.numeric(n_threads)
    if(is.na(n_threads_to_use)) {
        stop(paste(
            "n_threads must be a numeric value"
        ), call. = FALSE)
    }
    if(is.na(n_threads_to_use)) {
        stop(paste(
            "n_threads must be a numeric value"
        ), call. = FALSE)
    }
    if(n_threads_to_use > (parallel::detectCores() - 2) ) {
        n_threads_to_use = max(1, parallel::detectCores() - 2)
    }
    if(Sys.info()["sysname"] == "Windows") {
        BPPARAM_mod = BiocParallel::SnowParam(n_threads_to_use)
        message(paste("Using SnowParam", BPPARAM_mod$workers, "threads"))
    } else {
        BPPARAM_mod = BiocParallel::MulticoreParam(n_threads_to_use)
        message(paste("Using MulticoreParam", BPPARAM_mod$workers, "threads"))
    }
    return(BPPARAM_mod)
}

NxtIRF.SplitVector <- function(vector = "", n_workers = 1) {
    if(!is.numeric(n_workers) || n_workers < 1) {
        stop(paste(
            "n_workers must be at least 1"
        ), call. = FALSE)
    }
    n_workers_use = as.integer(n_workers)
    if(length(vector) < 1) {
        stop(paste(
            "vector to split must be of length at least 1"
        ), call. = FALSE)
    }
    
    if(n_workers_use > length(vector)) n_workers_use = length(vector)
    vector_starts = round(seq(1, length(vector) + 1, 
        length.out = n_workers_use + 1))
    vector_starts = unique(vector_starts)
    
    return_val = list()
    for(i in seq_len(length(vector_starts) - 1)) {
        return_val[[i]] = vector[seq(vector_starts[i], vector_starts[i+1] - 1)]
    }
    return(return_val)
}

semi_join.DT = function(A, B, by, nomatch = 0) {
    A[A[B, on = by, which = TRUE, nomatch = nomatch]]
}

.grDT <- function(DT, ...) {
    # Converts data table to GRanges object, preserving info
    makeGRangesFromDataFrame(as.data.frame(DT), ...)
}

.grlGaps <- function(grl) {
    psetdiff(unlist(range(grl), use.names = TRUE), grl)
}



make.path.relative = function(base, target) {
    if(Sys.info()["sysname"] == "Windows") {
        base = normalizePath(base, winslash = "/")
    }
    common = sub('^([^|]*)[^|]*(?:\\|\\1[^|]*)$', '^\\1/?', 
        paste0(base, '|', target))
    
    paste0(gsub('[^/]+/?', '../', sub(common, '', base)),
        sub(common, '', target))
}

# SQLITE operations

.open_sqlite <- function(sqlite_file) {
    return(DBI::dbConnect(SQLite(), sqlite_file))
}

.close_sqlite <- function(conn) {
    DBI::dbDisconnect(conn)
}

.to_sqlite <- function(sqldb, df, DB_name,...) {
    tryCatch({
        RSQLite::dbWriteTable(sqldb, DB_name, df, ...)
    }, error = function(e) {
        .close_sqlite(sqldb)
        stop(paste(
            "Writing table:", DB_name, "to sqlite failed"
        ), call. = FALSE)
    })
}

.from_sqlite <- function(sqldb, DB_name, DT = FALSE, ...) {
    df = NULL
    tryCatch({
        df = RSQLite::dbReadTable(sqldb, DB_name, ...)
    }, error = function(e) {
        .close_sqlite(sqldb)
        stop(paste(
            "Reading table:", DB_name, "from sqlite failed"
        ), call. = FALSE)
    })
    if(DT == TRUE) {
        return(as.data.table(df))
    } else {
        return(df)
    }
}

#' GGPLOT themes
#'
#' A ggplot theme object for white background figures without a legend
#' @examples
#' library(ggplot2)
#' df <- data.frame(
#'   gp = factor(rep(letters[1:3], each = 10)),
#'   y = rnorm(30))
#' ggplot(df, aes(gp, y)) +
#'   geom_point() +
#'   theme_white
#' @export
theme_white = theme(axis.line.x = element_line(colour = "black"),
            panel.grid.major = element_line(size = rel(0.5), colour="grey"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            legend.position = "none",
            axis.title.x.top = element_blank(),
            # axis.text.x.bottom = element_blank(),
            # axis.text.y = element_blank(),
            axis.line.x.bottom = element_blank(),
            axis.text=element_text(size=rel(1.0)),
            plot.title = element_text(hjust = 0.5),
            # axis.title.x=element_blank(),
            # axis.title.y=element_blank()
            )

#' @describeIn theme_white White theme but with a figure legend (if applicable)
#' @export
theme_white_legend = theme(axis.line.x = element_line(colour = "black"),
            panel.grid.major = element_line(size = rel(0.5), colour="grey"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            # legend.position = "none",
            axis.title.x.top = element_blank(),
            # axis.text.x.bottom = element_blank(),
            # axis.text.y = element_blank(),
            axis.line.x.bottom = element_blank(),
            axis.text=element_text(size=rel(1.0)),
            plot.title = element_text(hjust = 0.5),
            # axis.title.x=element_blank(),
            # axis.title.y=element_blank()
            )

dash_progress <- function(message = "", total_items = 1, add_msg = FALSE) {
    if(total_items != round(total_items) | total_items < 1) {
        stop(paste(
            "dash_progress needs at least 1 item"
        ), call. = FALSE)
    }
    if(add_msg) {
        message(message)
    }
    if(!is.null(shiny::getDefaultReactiveDomain())) {
        shiny::incProgress(1/total_items, message = message)
    }
}

#' Retrieves a data frame containing names and paths of example BAM files
#'
#' Intended to be run using the NxtIRF mock reference genome and gene 
#' annotations
#' @return A data frame, with the columns being the sample names and paths
#'   of the example BAM dataset aligned to the mock chromosome
#' @examples
#' bams = NxtIRF_example_bams()
#' @export
NxtIRF_example_bams <- function() {
    Find_Bams(dirname(example_bams())[1])
}

#' Retrieves a NxtSE object created from NxtIRF's example bam and mock reference
#'
#' @return A NxtSE object returned after creating and using the mock reference
#'   to run IRFinder, CollateData and MakeSE on the example bam files
#' @examples
#' se = NxtIRF_example_NxtSE()
#' @seealso [NxtIRF_example_bams()], [mock_genome()], [mock_gtf()]
#' @seealso [BuildReference()], [IRFinder()], [CollateData()], [MakeSE()]
#' @export
NxtIRF_example_NxtSE <- function() {
    se = readRDS(system.file("extdata", 
        "example_NxtSE.Rds", package = "NxtIRF"))
    covs = FindSamples(system.file("extdata", package = "NxtIRF"), ".cov")
    covfile(se) <- covs$path
    se
}
