#' NxtIRF Examples
#'
#' Contains files that provides a workable example for the 
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
#' @param destination_path The directory to place the downloaded example files.
#'   If left blank, gives the direct BiocFileCache path of the resource.
#' @return See Functions section below.
#' @examples
#' # returns the location of the genome.fa file of the mock reference
#'
#' mock_genome() 
#'
#' # returns the location of the transcripts.gtf file of the mock reference
#'
#' mock_gtf() 
#' 
#' # returns the location of the Mappability exclusion BED for hg38
#'
#' get_mappability_exclusion("hg38") 
#'
#' # returns the locations of the 6 example bam files
#'
#' example_bams() 
#'
#' # returns a data frame with the first column as sample names, and the 
#' # second column as BAM paths
#'
#' NxtIRF_example_bams() 
#'
#' # Returns a NxtSE object created by the example bams aligned to the
#'
#' # mock NxtSE reference
#'
#' NxtIRF_example_NxtSE() 
#'

#'
#' @references
#' Generation of the mappability files was performed using NxtIRF using
#' a method analogous to that described in:
#' 
#' Middleton R, Gao D, Thomas A, Singh B, Au A, Wong JJ, Bomane A, Cosson B, 
#' Eyras E, Rasko JE, Ritchie W.
#' IRFinder: assessing the impact of intron retention on mammalian gene 
#' expression.
#' Genome Biol. 2017 Mar 15;18(1):51.
#' \doi{10.1186/s13059-017-1184-4}
#' @name example-NxtIRF-data
#' @aliases 
#' mock_genome mock_gtf example_bams NxtIRF_example_bams NxtIRF_example_NxtSE
#' get_mappability_exclusion
#' @keywords package
#' @seealso [BuildReference()], [IRFinder()], [CollateData()], [MakeSE()]
#' @md
NULL


res_url <- "https://raw.github.com/alexchwong/NxtIRFdata/main/inst/NxtIRF"

#' @describeIn example-NxtIRF-data Makes a copy of the NxtIRF mock genome FASTA
#'   file, in the given destination path
#' @export
mock_genome <- function(destination_path = tempdir())
{
    .cache_and_create_file(paste(res_url, "genome.fa", sep="/"),
        destination_path, hub = "AnnotationHub")
}

#' @describeIn example-NxtIRF-data Makes a copy of the NxtIRF mock gene 
#'   annotation GTF file, in the given destination path
#' @export
mock_gtf <- function(destination_path = tempdir())
{
    .cache_and_create_file(paste(res_url, "transcripts.gtf", sep="/"),
        destination_path, hub = "AnnotationHub")
}

#' @describeIn example-NxtIRF-data Returns a copy of Mappability Exclusion 
#'   BED file for the specified genome, in the given destination path
#' @export
get_mappability_exclusion <- function(
        genome_type = c("hg38", "hg19", "mm10", "mm9"),
        destination_path = tempdir()) {
    genome_type = match.arg(genome_type)
    if(genome_type == "hg38") {
        .cache_and_create_file(paste(res_url, 
            "Mappability_Regions_hg38_v94.txt.gz", sep="/"),
            destination_path, hub = "AnnotationHub")
    } else if(genome_type == "hg19") {
        .cache_and_create_file(paste(res_url, 
            "Mappability_Regions_hg19_v75.txt.gz", sep="/"),
            destination_path, hub = "AnnotationHub")    
    } else if(genome_type == "mm10") {
        .cache_and_create_file(paste(res_url, 
            "Mappability_Regions_mm10_v94.txt.gz", sep="/"),     
            destination_path, hub = "AnnotationHub")    
    } else if(genome_type == "mm9") {
        .cache_and_create_file(paste(res_url, 
            "Mappability_Regions_mm9_v67.txt.gz", sep="/"),    
            destination_path, hub = "AnnotationHub")    
    } else {
        .log(paste("In get_mappability_exclusion():",
            "genome_type = ", genome_type, "is not recogised"))
    }
}

#' @describeIn example-NxtIRF-data Makes a copy of the Leucegene example 
#'   BAM files, aligned to the NxtIRF mock genome, in the given destination path
#' @export
example_bams <- function(destination_path = tempdir())
{
    bams = c("02H003.bam", "02H025.bam", "02H026.bam",
        "02H033.bam", "02H043.bam", "02H046.bam")
    .cache_and_create_file(
        paste(res_url, bams, sep="/"),
        destination_path, hub = "ExperimentHub")
}

#' @describeIn example-NxtIRF-data Returns a 2-column data frame, containing 
#'   sample names and sample paths (in tempdir()) of example BAM files
#' @export
NxtIRF_example_bams <- function() {
    Find_Bams(dirname(example_bams())[1])
}

#' @describeIn example-NxtIRF-data Returns a (in-memory / realized) NxtSE object 
#'   generated using the NxtIRF mock reference and example BAM files
#' @export
NxtIRF_example_NxtSE <- function() {
    se = readRDS(system.file("extdata", 
        "example_NxtSE.Rds", package = "NxtIRF"))
    covs = Find_Samples(system.file("extdata", package = "NxtIRF"), ".cov")
    covfile(se) <- covs$path
    se
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
            } else if (((httr::status_code(res$result) %/% 200) != 1)) {
                if (!quiet) warning(
                    sprintf(paste(
                        "Requests for [%s] responded but without an",
                        "HTTP status code in the 200-299 range"
                    ), url))
                ret = c(ret, non_2xx_return_value)
            } else {
                ret = c(ret, TRUE)
            }
        } else {
            ret = c(ret, TRUE)
        }
    }
    return(ret)
}

.fetch_files_from_urls <- function(urls, filenames, destination_path,
        bfc_only = FALSE) {
    # No error checking
    cache <- tools::R_user_dir(package = "NxtIRF", which="cache")
    bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
    files = c()
    paths = c()
    for(i in seq_len(length(urls))) {
        url = urls[i]
        file = filenames[i]
        test <- BiocFileCache::bfcquery(bfc, url, field = "rname")
        if(length(test$rpath) > 1) {
            # Corrupt cache, delete all and re-download
            BiocFileCache::bfcremove(bfc, test$rid)
            test <- BiocFileCache::bfcquery(bfc, url, field = "rname")
            if(length(test$rpath) > 0) {
                .log(paste("Corrupt BiocFileCache for NxtIRF:",
                    url, "multiple records exist and undeletable"))
            }
        }
        if(!bfc_only && length(test$rpath) < 1) {
            if(.check_if_url_exists(url)) {
                path <- tryCatch(BiocFileCache::bfcrpath(bfc, url),
                error = function(e) {
                    .log(paste("Download from url failed:", url, 
                        "- please check connection"))
                })
            } else {
                .log(paste("url not found:", url, 
                    "- please check connection"))
            }
        } else if(length(test$rpath) == 1) {
            path <- test$rpath
        } else {
            return(NULL)
        }
        paths = c(paths, path)
    }
    
    if(length(paths) == length(filenames)) {
        for(i in seq_len(length(paths))) {
            path = paths[i]
            file = filenames[i]
            if(destination_path == "") {
                files = c(files, path)
            } else {
                dest = file.path(destination_path, file)
                if(file.exists(dest) && 
                        identical(openssl::md5(file(path)), 
                            openssl::md5(file(dest)))) {
                    # if md5 identical, do nothing
                } else {
                    file.copy(path, dest, overwrite = TRUE)
                }
                files = c(files, dest)
            }
        }
    }    
    return(normalizePath(files))
}

.cache_and_create_file <- function(urls, destination_path = "",
        filenames = basename(urls), hub = c("AnnotationHub", "ExperimentHub")) {
    if(destination_path != "" && !dir.exists(dirname(destination_path))) {
        warning(paste(destination_path, "- parent directory does not exist"))
        return("")
    } else if(destination_path != "" && any(filenames %in% "")) {
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
    # Check NxtIRF hub first:
    files_from_own_hub <- .fetch_files_from_urls(
        urls, filenames, destination_path, bfc_only = TRUE)
    if(!is.null(files_from_own_hub)) return(files_from_own_hub)
    # Check if in AnnotationHub / ExperimentHub:   
    hub = match.arg(hub)
    if(hub == "AnnotationHub") {
        hub_obj <- AnnotationHub::AnnotationHub()
    } else if(hub == "ExperimentHub") {
        hub_obj <- ExperimentHub::ExperimentHub()    
    }
    hub_obj = hub_obj[match(urls, hub_obj$sourceurl)]
    if(length(hub_obj) == length(urls)) {
        return(.cache_and_create_file_hub(urls, hub_obj, 
            destination_path, filenames))
    } else {
        return(.fetch_files_from_urls(urls, filenames, destination_path))
    }
}

.cache_and_create_file_hub <- function(urls, 
        hub = ExperimentHub::ExperimentHub(),
        destination_path = "", filenames = basename(urls)) {
    hub = hub[match(urls, hub$sourceurl)]
    files = c()
    for(i in seq_len(length(hub))) {
        cache = unname(cache(hub[i])[1])     # first file is bam
        file = filenames[i]
        if(destination_path == "") {
            files = c(files, cache)
        } else {
            dest = file.path(destination_path, file)
            if(file.exists(dest) && 
                    identical(openssl::md5(file(dest)), 
                        openssl::md5(file(cache)))) {
                # if md5 identical, do nothing
            } else {
                file.copy(cache, dest, overwrite = TRUE)
            }
            files = c(files, dest)
        }
    }
    # Make a copy to local hub for quicker access next time
    # (as AnnotationHub and ExperimentHub are slow!)
    cache <- tools::R_user_dir(package = "NxtIRF", which="cache")
    bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
    for(i in seq_len(length(urls))) {
        url = urls[i]
        file = files[i]
        test <- BiocFileCache::bfcquery(bfc, url, field = "rname")
        if(length(test$rpath) == 0) {
            BiocFileCache::bfcadd(bfc, url, fpath = file)
        }
    }
    return(normalizePath(files))
}