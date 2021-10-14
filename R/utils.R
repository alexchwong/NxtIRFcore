# Exported miscellaneous functions

#' Converts genomic coordinates into a GRanges object
#'
#' This function takes a string vector of genomic coordinates and converts it
#' into a GRanges object.
#'
#' @details
#' Genomic coordinates can take one of the following syntax:
#' * `seqnames:start`
#' * `seqnames:start-end`
#' * `seqnames:start-end/strand`
#'
#' The following examples are considered valid genomic coordinates:
#' * "chr1:21535"
#' * "chr3:10550-10730"
#' * "X:51231-51330/-"
#' * "chrM:2134-5232/+"
#' @param coordinates A string vector of one or more genomic coordinates
#'   to be converted
#' @return A GRanges object that corresponds to the given coordinates
#' @examples
#' se <- NxtIRF_example_NxtSE()
#'
#' coordinates <- rowData(se)$EventRegion
#'
#' gr <- CoordToGR(coordinates)
#' @md
#' @export
CoordToGR <- function(coordinates) {
    stopmsg <- paste(
        "Coordinates must take the form chrN:X, chrN:X-Y,",
        "chrN:X-Y/+ or chrN:X-Y/-"
    )
    temp <- tstrsplit(coordinates, split = "/")
    if (length(temp) == 2) {
        strand <- as.character(temp[[2]])
    } else if (length(temp) == 1) {
        strand <- "*"
    } else {
        .log(stopmsg)
    }
    temp2 <- tstrsplit(temp[[1]], split = ":")
    if (length(temp2) == 2) {
        seqnames <- temp2[[1]]
        temp3 <- tstrsplit(temp2[[2]], split = "-")
        start <- as.numeric(temp3[[1]])
        if (length(temp3) == 2) {
            end <- as.numeric(temp3[[2]])
        } else if (length(temp3) == 1) {
            end <- as.numeric(temp3[[1]])
        } else {
            .log(stopmsg)
        }
    } else {
        .log(stopmsg)
    }
    if(length(start) != length(end))
        .log("In CoordToGR, start and end lengths don't match")
    if(any(start > end))
        .log("In CoordToGR, some starts > end")
    return(GRanges(
        seqnames = seqnames, 
        ranges = IRanges(start = start, end = end),
        strand = strand
    ))
}

#' Validates the given file as a valid COV file
#'
#' This function takes the path of a possible COV file and checks whether its
#' format complies with that of the COV format defined by this package.
#'
#' @details
#' COV files are BGZF-compressed files. The first 4 bytes of the file must
#' always be 'COV\1', distinguishing it from BAM or other files in BGZF format.
#' This function checks whether the given file complies with this.
#'
#' @param coverage_files A vector containing the file names of files to be
#'   checked
#' @return `TRUE` if all files are valid COV files. `FALSE` otherwise
#' @examples
#' se <- NxtIRF_example_NxtSE()
#'
#' cov_files <- covfile(se)
#'
#' IsCOV(cov_files) # returns true if these are true COV files
#' @seealso [IRFinder] [CollateData]
#' @md
#' @export
IsCOV <- function(coverage_files) {
    for (i in coverage_files) {
        if (file.exists(i) && IRF_Check_Cov(normalizePath(i))) {
            # do nothing
        } else {
            return(FALSE)
        }
    }
    return(TRUE)
}
