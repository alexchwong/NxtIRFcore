# Miscellaneous internal wrappers to NxtIRFcore.dll

# Simple unzip function
# To check gunzip produces same output using NxtIRF vs other utilities
run_Gunzip = function(infile = "", outfile) {
    file_to_read = normalizePath(infile)
    if(!file.exists(file_to_read)) {
        .log(paste("In run_Gunzip(),",
            file_to_read, "does not exist"))
    }
    if(!dir.exists(dirname(outfile))) {
        .log(paste("In run_Gunzip(),",
            dirname(outfile), "does not exist"))
    }
    IRF_gunzip(file_to_read, outfile)
}

# Gets a specific data frame in a gzipped multi-tabular text file
# If getting a small data frame situated at the beginning of a large file
#   e.g. in IRFinder output, this is typically faster than data.table::fread
get_multi_DT_from_gz = function(infile = "", 
        block_headers = c("Header1", "Header2")) {
    file_to_read = normalizePath(infile)
    if(!file.exists(file_to_read)) {
        .log(paste("In get_multi_DT_from_gz(),",
            file_to_read, "does not exist"))
    }
    df.list = IRF_gunzip_DF(file_to_read, block_headers)
    for(i in seq_len(length(df.list))) {
        for(j in seq_len(length(df.list[[i]]))) {
            suppressWarnings({
                if(all(df.list[[i]][[j]] == "NA" | 
                        !is.na(as.numeric(df.list[[i]][[j]])))) {
                    df.list[[i]][[j]] = as.numeric(df.list[[i]][[j]])
                }
            })
        }
        df.list[[i]] = as.data.table(df.list[[i]])
    }
    return(df.list)
}


#' Validates the given file as a valid COV file
#' @param coverage_files A vector containing the file names of files to be
#'   checked
#' @return `TRUE` if all files are valid COV files. `FALSE` otherwise
#' @examples
#' se = NxtIRF_example_NxtSE()
#'
#' cov_files = covfile(se)
#'
#' IsCOV(cov_files) # returns true if these are true COV files
#' @md
#' @export
IsCOV = function(coverage_files) {
    for(i in coverage_files) {
        if(file.exists(i) && IRF_Check_Cov(normalizePath(i))) {
            # do nothing
        } else {
            return(FALSE)
        }
    }    
    return(TRUE)
}