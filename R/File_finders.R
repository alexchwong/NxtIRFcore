#' Convenience Function to (recursively) find all files in a folder.
#' 
#' Often, output files (whether it be raw sequencing files, aligned sequences)
#' in BAM files, or IRFinder output files, are stored in a single folder.
#' Sometimes, these are named according to the samples they represent. Other
#' times, they have generic names but are partitioned in sub-folders named
#' by sample names. This function (recursively) finds all files and extracts
#' sample names assuming the files are named by sample names 
#' (`level = 0`). Alternately, the names can be derived from the
#' subfolders one level higher (`level = 1`)
#' 
#' @param sample_path The path in which to recursively search for files
#'   that match the given `suffix`
#' @param suffix A vector of or or more strings that specifies the file suffix 
#'   (e.g. 'bam' denotes BAM files, whereas "txt.gz" denotes gzipped txt 
#'   files).
#' @param level Whether sample names can be found in the file names themselves
#'   (level = 0), or their parent directory (level = 1). Potentially parent
#'   of parent directory (level = 2). Support max level <= 3 (for sanity).
#' @param paired Whether to expect single FASTQ files (of the format
#'   "sample.fastq"), or
#'   paired files (of the format "sample_1.fastq", "sample_2.fastq")
#' @param fastq_suffix The name of the FASTQ suffix. Options are:
#'   "fastq", "fastq.gz", "fq", or "fq.gz"
#' @return A multi-column data frame with the first column containing
#'   the sample name, and subsequent columns being the file paths with suffix
#'   as determined by `suffix`.
#' @examples
#' # Retrieve all BAM files in a given folder, named by sample names
#' bam_path = tempdir()
#' example_bams(path = bam_path)
#' df.bams = Find_Samples(sample_path = bam_path, 
#'   suffix = ".bam", level = 0)
#' # equivalent to:
#' df.bams = Find_Bams(bam_path, level = 0)
#'
#' # Retrieve all IRFinder output files in a given folder, 
#' # named by sample names
#'
#' expr = Find_IRFinder_Output(file.path(tempdir(), "IRFinder_output"))
#' @name Find_Samples
#' @md
#' @export
Find_Samples <- function(sample_path, suffix = ".txt.gz", level = 0) {
    if(length(suffix) == 0)
        .log(paste("In Find_Samples(),",
            "suffix must be of length greater than zero"))
    if(!dir.exists(sample_path))
        .log(paste("In Find_Samples(),", 
            sample_path, "- given path does not exist"))
    if(length(level) > 1) .log("level must be a numeral of length 1")
    if(!is.numeric(level)) .log("In Find_Samples(), level must be numeric")
    if(level < 0) .log("In Find_Samples(), level must be non-negative")

    level = floor(level)
    suffix_name = "path"
    if(length(suffix) > 1) suffix_name = suffix

    DT.list = list()
    for(i in seq_len(length(suffix))) {
        pattern = paste0("\\", suffix[i], "$")
        files_found = list.files(pattern = pattern,
            path = normalizePath(sample_path), 
            full.names = TRUE, recursive = TRUE)
        if(length(files_found) > 0) {
            DT = data.table(sample = gsub(pattern,"",
                    files_found), 
                path = files_found)
            lvl = level
            while(lvl > 0) {
                DT$sample = dirname(DT$sample)
                if(any(DT$sample %in% c(".", "/", "~")))
                    .log(paste("Sample name points to", DT$sample, 
                        "- please check level of sample names"))
                lvl = lvl - 1
            }
            DT$sample = basename(DT$sample)
            colnames(DT)[2] = suffix_name[i]
            setkeyv(DT, "sample")
            DT.list[[i]] = DT
        } else {
            DT.list[[i]] = NULL
        }
    }
    if(length(DT.list) <= 1) return(as.data.frame(DT.list))
    final = DT.list[[1]]
    if(length(suffix) > 1) {
        # Check identity of sorted sample names
        samples = DT.list[[1]]$sample
        for(i in seq(2, length(suffix))) {
            if(!is.null(DT.list[[i]]) && 
                    identical(samples, DT.list[[i]]$sample)) {
                final = cbind(final, DT.list[[i]][, 2, with = FALSE])
            }
        }
    }
    cols = c("sample", suffix_name[suffix_name %in% colnames(final)])
    return(as.data.frame(final[, cols, with = FALSE]))
}

#' @describeIn Find_Samples Returns all FASTQ files in a given folder
#' @export
Find_FASTQ <- function(sample_path, paired = TRUE, 
        fastq_suffix = c("fastq", "fq", "fastq.gz", "fq.gz"), level = 0) {
    fastq_suffix = match.arg(fastq_suffix)
    suffix_use = ifelse(paired, paste0(c("_1.", "_2."), fastq_suffix), 
        paste0(".", fastq_suffix))
    DT = Find_Samples(sample_path, suffix_use, level = level)
    colnames(DT)[2] = "forward"
    if(paired) colnames(DT)[3] = "reverse"
    return(DT)
}

#' @describeIn Find_Samples Returns all BAM files in a given folder
#' @export
Find_Bams <- function(sample_path, level = 0) {
    return(Find_Samples(sample_path, ".bam", level = level))
}

#' @describeIn Find_Samples Returns all IRFinder output files in a given folder
#' @export
Find_IRFinder_Output <- function(sample_path, level = 0) {
    DT = Find_Samples(sample_path, c(".txt.gz", ".cov"), level = level)
    colnames(DT) = c("sample", "irf_file", "cov_file")
    return(DT)
}
