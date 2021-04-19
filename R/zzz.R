#' NxtIRF Examples
#'
#' This contains files that provides a workable example for the 
#' NxtIRF package.\cr\cr
#' A mock reference, with genome sequence (FASTA) and gene annotation (GTF)
#' files are provided, based on the genes SRSF1, SRSF2, SRSF3, TRA2A, TRA2B, 
#' TP53 and NSUN5, of which sequences are used to construct an artificial 
#' chromosome Z. This was generated based on release-94 of Ensembl GRCh38 (hg38)
#' reference.\cr\cr
#' Also, there are 6 example bam files based on samples from the 
#' Leucegene dataset (GSE67039). Bam files are constructed
#' based on the complete bam files of 6 samples from Leucegene,
#' subsetted by regions containing the 7 above genes. Then, the reads of these 
#' subsetted BAMs were realigned to the mock reference using STAR.\cr\cr
#' Additionally, there are files for Mappability exclusion regions generated
#' using NxtIRF, suitable for use in generating references based on hg38,
#' hg19, mm10 and mm9 genomes.
#' @param genome_type Either one of `hg38`, `hg19`, `mm10` or `mm9`
#' @param destination_path The directory to place the downloaded example files
#' @return See Examples section below.
#' @examples
#' # returns the location of the genome.fa file of the mock reference
#' mock_genome() 
#'
#' # returns the location of the transcripts.gtf file of the mock reference
#' mock_gtf() 
#'
#' # returns the locations of the 6 example bam files
#' example_bams() 
#'
#' # returns the location of the Mappability exclusion BED for hg38
#' get_mappability_exclusion("hg38") 
#' @references
#' Generation of the mappability files was performed using NxtIRF using
#' a method analogous to that described in:
#' 
#' Middleton R, Gao D, Thomas A, Singh B, Au A, Wong JJ, Bomane A, Cosson B, 
#' Eyras E, Rasko JE, Ritchie W.
#' IRFinder: assessing the impact of intron retention on mammalian gene 
#' expression.
#' Genome Biol. 2017 Mar 15;18(1):51.
#' \url{https://doi.org/10.1186/s13059-017-1184-4}
#' @name NxtIRF-example-data
#' @aliases 
#' mock_genome
#' mock_gtf
#' example_bams
#' get_mappability_exclusion
#' @keywords package
#' @md
NULL

res_url <- "https://raw.github.com/alexchwong/NxtIRFdata/main/inst/NxtIRF"

#' @export
mock_genome <- function(destination_path = tempdir())
{
    .cache_and_create_file(paste(res_url, "genome.fa", sep="/"),
        destination_path)
}

#' @export
mock_gtf <- function(destination_path = tempdir())
{
    .cache_and_create_file(paste(res_url, "transcripts.gtf", sep="/"),
        destination_path)
}

#' @export
example_bams <- function(destination_path = tempdir())
{
    bams = c("02H003.bam", "02H025.bam", "02H026.bam",
        "02H033.bam", "02H043.bam", "02H046.bam")
    .cache_and_create_file(
        paste(res_url, bams, sep="/"),
        destination_path)
}

#' @export
get_mappability_exclusion <- function(
        genome_type = c("hg38", "hg19", "mm10", "mm9"),
        destination_path = tempdir()) {
    genome_type = match.arg(genome_type)
    if(genome_type == "hg38") {
        .cache_and_create_file(paste(res_url, 
            "Mappability_Regions_hg38_v94.txt.gz", sep="/"),
            destination_path)
    } else if(genome_type == "hg19") {
        .cache_and_create_file(paste(res_url, 
            "Mappability_Regions_hg19_v75.txt.gz", sep="/"),
            destination_path)    
    } else if(genome_type == "mm10") {
        .cache_and_create_file(paste(res_url, 
            "Mappability_Regions_mm10_v94.txt.gz", sep="/"),     
            destination_path)    
    } else if(genome_type == "mm9") {
        .cache_and_create_file(paste(res_url, 
            "Mappability_Regions_mm9_v67.txt.gz", sep="/"),    
            destination_path)    
    } else {
        stop(paste("In get_mappability_exclusion():",
            "genome_type = ", genome_type, "is not recogised"
        ), call. = FALSE)
    }
}

.capture_error <- function(code, otherwise = NULL, quiet = TRUE) {
    tryCatch(
        list(result = code, error = NULL),
        error = function(e) {
            if (!quiet) message("Error: ", e$message)
            list(result = otherwise, error = e)
        },
        interrupt = function(e) {
            stop("Terminated by user", call. = FALSE)
        }
    )
}
.safely <- function(.f, otherwise = NULL, quiet = TRUE) {
    function(...) .capture_error(.f(...), otherwise, quiet)
}

.sHEAD <- .safely(httr::HEAD)
.sGET <- .safely(httr::GET)

.check_if_url_exists <- function(urls, ...) {
    # from https://stackoverflow.com/questions/52911812/check-if-url-exists-in-r
    non_2xx_return_value = FALSE
    quiet = FALSE
    
    # Try HEAD first since it's lightweight
    ret = c()
    for(url in urls) {
        res <- .sHEAD(url, ...)
        if (is.null(res$result) || 
                ((httr::status_code(res$result) %/% 200) != 1)) {
            res <- .sGET(url, ...)
            if (is.null(res$result)) {
                if (!quiet) warning("Invalid web link")
                ret = c(ret, FALSE)
            }
            # or whatever you want to return on "hard" errors
            if (((httr::status_code(res$result) %/% 200) != 1)) {
                if (!quiet) warning(
                    sprintf(paste(
                        "Requests for [%s] responded but without an",
                        "HTTP status code in the 200-299 range"
                    ), url))
                ret = c(ret, non_2xx_return_value)
            }
            ret = c(ret, TRUE)
        } else {
            ret = c(ret, TRUE)
        }
    }
    return(ret)
}

.fetch_files_from_urls <- function(urls, filenames, destination_path) {
    # No error checking
    cache <- rappdirs::user_cache_dir(appname = "NxtIRF")
    bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
    files = c()
    for(i in seq_len(length(urls))) {
        url = urls[i]
        file = filenames[i]
        test <- BiocFileCache::bfcquery(bfc, url, field = "rname")
        if(length(test$rpath) > 1) {
            # Corrupt cache, delete all and re-download
            BiocFileCache::bfcremove(bfc, test$rid)
            test <- BiocFileCache::bfcquery(bfc, url, field = "rname")
            if(length(test$rpath) > 0) {
                stop(paste("Corrupt BiocFileCache for NxtIRF:",
                    url, "multiple records exist and undeletable"
                ), call. = FALSE)
            }
        }
        if(length(test$rpath) < 1) {
            if(.check_if_url_exists(url)) {
                path <- tryCatch(BiocFileCache::bfcrpath(bfc, url),
                error = function(e) {
                    stop(paste("Download from url failed:", url, 
                        "- please check connection"
                    ), call. = FALSE)
                })
            } else {
                stop(paste("url not found:", url, 
                    "- please check connection"
                ), call. = FALSE)
            }
        } else {
            path <- test$rpath
        }
        
        if(destination_path == "") {
            files = c(files, path)
        } else {
            dest = file.path(destination_path, file)
            if(file.exists(dest) && 
                    identical(openssl::md5(file(path)), 
                        openssl::md5(file(dest)))) {
                # if md5 identical, do nothing
            } else {
                file.copy(path,
                    file.path(destination_path, file),
                    overwrite = TRUE
                )
            }
            files = c(files, file.path(destination_path, file))
        }
    }
    return(normalizePath(files))
}

.cache_and_create_file <- function(urls, destination_path = "",
        filenames = basename(urls)) {
    if(destination_path != "" && !dir.exists(dirname(destination_path))) {
        warning(paste(destination_path, "- parent directory does not exist"))
        return("")
    } else if(destination_path != "" && any(filenames == "")) {
        warning("Invalid return filename(s) for web resource download")
        return("")
    } else if(length(urls) == 0 && length(urls) != length(filenames)) {
        warning("length of 'urls' and 'filenames' must be the same")
        return("")
    }
    if(!any(startsWith(urls, "http") | startsWith(urls, "ftp"))) {
        warning(paste(urls, "some urls does not appear to be valid web links"))
    }
    if(destination_path != "" && !dir.exists(destination_path)) {
        tryCatch(dir.create(destination_path),
            error = function(e) {
                warning(paste("Unable to create", destination_path))
                return("")
            }
        )
    }
    return(.fetch_files_from_urls(urls, filenames, destination_path))
}

