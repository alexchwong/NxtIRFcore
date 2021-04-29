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
STAR_align_experiment <- function(Experiment, STAR_ref_path, BAM_output_path,
        trim_adaptor = "AGATCGGAAG", two_pass = FALSE, n_threads = 4) {
    .validate_STAR_version()
    STAR_ref_path = .validate_STAR_reference(STAR_ref_path)
    BAM_output_path = .validate_path(BAM_output_path)
    
    # Dissect Experiment:
    if(ncol(Experiment) < 2 || ncol(Experiment) > 3) {
        stop(paste(
            "Experiment must be a 2- or 3- column data frame,",
            "with the columns denoting sample name, fastq file (forward),",
            "and (optionally) fastq file (reverse)"
        ), call. = FALSE)
    } else if(ncol(Experiment) == 2) {
        colnames(Experiment) = c("sample", "forward")
        fastq_1 = Experiment[, "forward"]
        fastq_2 = NULL
        .validate_STAR_fastq_samples(fastq_1)
        paired = FALSE
    } else if(ncol(Experiment) == 3) {
        colnames(Experiment) = c("sample", "forward", "reverse")
        fastq_1 = Experiment[, "forward"]
        fastq_2 = Experiment[, "reverse"]
        .validate_STAR_fastq_samples(fastq_1, fastq_2)
        paired = TRUE
    }
    gzipped = all(grepl(paste0("\\", ".gz", "$"), fastq_1)) &&
        (!paired || all(grepl(paste0("\\", ".gz", "$"), fastq_2)))
    if(is_valid(trim_adaptor)) .validate_STAR_trim_sequence(trim_adaptor)
    
    # system2(command = "STAR", args = c(
        # "--genomeLoad", "LoadAndExit", "--genomeDir", STAR_ref_path
    # ))
    
    samples = unique(Experiment[, "sample"])
    SJ.files = NULL
    two_pass_genome = NULL
    for(pass in seq_len(ifelse(two_pass, 2, 1))) {
        if(two_pass && pass == 1) message("STAR - first pass")
        if(two_pass && pass == 2) message("STAR - second pass")
        for(i in seq_len(length(samples))) {
            sample = samples[i]
            Expr_sample = Experiment[Experiment[, "sample"] == sample,]
            if(paired) {
                fastq_1 = Expr_sample[, "forward"]
                fastq_2 = NULL        
            } else {
                fastq_1 = Experiment[, "forward"]
                fastq_2 = Experiment[, "reverse"]
            }
            ref = STAR_ref_path
            memory_mode = "LoadAndKeep"
            if(two_pass && pass == 1) {
                additional_args = c("--outSAMtype", "None")
            } else if(two_pass && pass == 2 && !is.null(SJ.files)) {
                additional_args = c("--sjdbFileChrStartEnd",
                    paste(SJ.files, collapse = " "),
                    "--sjdbInsertSave", "All"
                )
                two_pass_genome = file.path(BAM_output_path, sample, 
                    "_STARgenome")
                SJ.files = NULL
                memory_mode = "LoadAndRemove"
            } else if(two_pass && pass == 2 && !is.null(two_pass_genome)) {
                ref = two_pass_genome
                additional_args = NULL
            } else {
                additional_args = NULL
            }

            message(paste("Aligning", sample, "using STAR"))
            STAR_align_fastq(ref, 
                BAM_output_path = file.path(BAM_output_path, sample),
                fastq_1 = fastq_1, fastq_2 = fastq_2, 
                trim_adaptor = trim_adaptor,
                memory_mode = memory_mode,
                additional_args = additional_args,
                n_threads = n_threads)
        }
        if(two_pass && pass == 1) {
            SJ.files = Find_Samples(BAM_output_path, suffix = ".out.tab")
            if(nrow(SJ.files) == 0) {
                stop(paste("In STAR two-pass,",
                    "no SJ.out.tab files were found"
                ), call. = FALSE)
            }
        }
    }

    system2(command = "STAR", args = c(
        "--genomeLoad", "Remove", "--genomeDir", ref
    ))
}

#' @export
STAR_align_fastq <- function(STAR_ref_path, BAM_output_path,
        fastq_1 = c("./sample_1.fq"), fastq_2 = NULL,
        two_pass = FALSE,
        trim_adaptor = "AGATCGGAAG",
        memory_mode = "NoSharedMemory",
        additional_args = NULL,
        n_threads = 4) {
        
    .validate_STAR_version()
    STAR_ref_path = .validate_STAR_reference(STAR_ref_path)
    .validate_STAR_fastq_samples(fastq_1, fastq_2)
    
    paired = (length(fastq_1) == length(fastq_2))
    gzipped = all(grepl(paste0("\\", ".gz", "$"), fastq_1)) &&
        (!paired || all(grepl(paste0("\\", ".gz", "$"), fastq_2)))
    if(is_valid(trim_adaptor)) .validate_STAR_trim_sequence(trim_adaptor)

    BAM_output_path = .validate_path(BAM_output_path)
    # Load STAR reference
    
    args = c(
        "--genomeLoad", memory_mode,
        "--runThreadN", .validate_threads(n_threads, as_BPPARAM = FALSE),
        "--genomeDir", STAR_ref_path,

        "--outFileNamePrefix", paste0(BAM_output_path, "/"),
        "--outStd", "Log",      # Not Bam_Unsorted
        
        "--outSAMtype", "BAM", "Unsorted", 
        "--outSAMstrandField", "intronMotif",
        "--outSAMunmapped", "None",

        "--outFilterMultimapNmax", "1"
    )
    if(two_pass) {
        args = c(args, "--twopassMode", "Basic")
    }
    
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
    if(!is.null(additional_args) && all(is.character(additional_args))) {
        args = c(args, additional_args)
    }
    system2(command = "STAR", args = args)
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