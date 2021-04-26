#' @export
STAR_buildRef <- function(reference_path, 
        STAR_ref_path = file.path(reference_path, "STAR"),
        sjdbOverlap = 150,
        n_threads = 4) {
    .validate_reference(reference_path)
    .validate_STAR_version()
    .validate_reference_path(STAR_ref_path)
    # Unzip reference files
    genome.fa <- .STAR_get_FASTA(reference_path)
    transcripts.gtf <- .STAR_get_GTF(reference_path)
    # Build STAR using defaults
    res = system2(command = "STAR", args = c(
        "--runMode", "genomeGenerate",
        "--genomeDir", STAR_ref_path,
        "--genomeFastaFiles", genome.fa,
        "--sjdbGTFfile", transcripts.gtf,
        "--sjdbOverhang", sjdbOverlap,
        "--runThreadN", .validate_threads(n_threads, as_BPPARAM = FALSE),
        "&>>", file.path(reference_path, "STAR", "log-star-build-ref.log")
    ), stdout = TRUE)
    if(res > 0) {
        warning("STAR Build Reference appears to have failed")
    }
}

.validate_STAR_version <- function() {
    if(Sys.info()["sysname"] != "Linux") {
        stop(paste("STAR only works on Linux"), call. = FALSE)
    }
    star_version = NULL
    tryCatch({
        star_version = system2("STAR", "--version", stdout = TRUE)
    }, error = function(e) {
        star_version = NULL
    })
    if(is.null(star_version)) {
        stop(paste("STAR is not installed"), call. = FALSE)
    }
    if(star_version < "2.5.0") {
        stop(paste("STAR version < 2.5.0 is not supported;",
            "current version:", star_version
        ), call. = FALSE)
    }
}

.STAR_get_FASTA <- function(reference_path) {
    genome.fa = file.path(reference_path, "resource", "genome.fa")
    if(!file.exists(paste0(genome.fa, ".gz"))) {
        stop(paste(paste0(genome.fa, ".gz"), "not found"), call. = FALSE)
    }
    gunzip(paste0(genome.fa, ".gz"), remove = FALSE)
    return(genome.fa)
}

.STAR_get_GTF <- function(reference_path) {
    transcripts.gtf = file.path(reference_path, "resource", "transcripts.gtf")
    if(!file.exists(paste0(transcripts.gtf, ".gz"))) {
        stop(paste(paste0(transcripts.gtf, ".gz"), "not found"), call. = FALSE)
    }
    gunzip(paste0(transcripts.gtf, ".gz"), remove = FALSE)
    return(transcripts.gtf)
}