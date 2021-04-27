#' @export
STAR_buildRef <- function(reference_path, 
        STAR_ref_path = file.path(reference_path, "STAR"),
        sjdbOverlap = 150,
        n_threads = 4) {
    .validate_reference(reference_path)
    .validate_STAR_version()
    .validate_path(STAR_ref_path)
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

#' @export
STAR_align_fastq <- function(STAR_ref_path, BAM_output_path,
        fastq_1 = c("./sample_1.fq"), fastq_2 = NULL,
        trim_adaptor = "AGATCGGAAG",
        memory_mode = "NoSharedMemory",
        
        n_threads = 4) {
    .validate_STAR_version()
    STAR_ref_path = .validate_STAR_reference(STAR_ref_path)
    .validate_STAR_fastq_samples(fastq_1, fastq_2)
    
    paired = (length(fastq_1) == length(fastq_2))
    gzipped = all(grepl(paste0("\\", ".gz", "$"), fastq_1)) &&
        (!paired || all(grepl(paste0("\\", ".gz", "$"), fastq_2)))
    if(is_valid(trim_adaptor)) .validate_STAR_trim_sequence(trim_adaptor)

    BAM_output_path = .validate_path(BAM_output_path)
    # Run STAR 
    args = c(
        "--genomeLoad", memory_mode,
        "--runThreadN", .validate_threads(n_threads, as_BPPARAM = FALSE),
        "--genomeDir", STAR_ref_path,

        "--outFileNamePrefix", paste0(BAM_output_path, "/"),
        "--outStd", "Log",      # Not Bam_Unsorted
        
        "--outSAMtype", shQuote("BAM Unsorted"), 
        "--outSAMstrandField", "intronMotif",
        "--outSAMunmapped", "None",

        "--outFilterMultimapNmax", "1"
    )
    args = c(args,
        "--readFilesIn",
        paste(fastq_1, collapse = ",")
    )
    if(paired) {
        args = c(args, paste(fastq_2, collapse = ","))
    }
    if(gzipped) {
        args = c(args, "--readFilesCommand", shQuote("gzip -dc"))
    }
    if(is_valid(trim_adaptor)) {
        args = c(args, "--clip3pAdapterSeq", trim_adaptor)
    }
    
    res = system2(command = "STAR", args = args, stdout = TRUE)
    
    print(res)
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

.validate_STAR_reference <- function(STAR_ref_path) {
    return(.validate_path(STAR_ref_path))
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

.validate_STAR_fastq_samples <- function(fastq_1, fastq_2) {
    if(!is_valid(fastq_2)) {
        # assume single
        if(!all(file.exists(fastq_1))) {
            stop(paste("Some fastq files were not found"), call. = FALSE)
        }
    } else {
        if(length(fastq_2) != length(fastq_1)) {
            stop(paste("There must be equal numbers of",
                "forward and reverse fastq samples"
            ), call. = FALSE)
        }
        if(!all(file.exists(fastq_1)) || !all(file.exists(fastq_2))) {
            stop(paste("Some fastq files were not found"), call. = FALSE)
        }
    }
    paired = (length(fastq_1) == length(fastq_2))
    gzipped = all(grepl(paste0("\\", ".gz", "$"), fastq_1)) &&
        (!paired || all(grepl(paste0("\\", ".gz", "$"), fastq_2)))
    if(!gzipped && 
        (
            any(grepl(paste0("\\", ".gz", "$"), fastq_1)) ||
            (paired && any(grepl(paste0("\\", ".gz", "$"), fastq_2)))
        )
    ) {
        stop(paste("A mixture of gzipped and uncompressed",
            "fastq files found.", "You must supply either all",
            "gzipped or all uncompressed fastq files"), call. = FALSE)
    }
}

.validate_STAR_trim_sequence <- function(sequence) {
    if(length(sequence) != 1) {
        stop(paste("Multiple adaptor sequences are not supported"),
        call. = FALSE)
    }
    tryCatch({
        ACGT_sum = sum(Biostrings::letterFrequency(
            Biostrings::DNAString(sequence), 
            letters = "AGCT", OR = 0))
    }, error = function(e) ACGT_sum = 0)
    if(nchar(sequence) != ACGT_sum) {
        stop(paste("Adaptor sequence can only contain A, C, G or T"),
        call. = FALSE)
    }
}