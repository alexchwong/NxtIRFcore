#' A wrapper function to call NxtIRF/IRFinder
#'
#' This function calls IRFinder on one or more BAM files.
#' @param bamfiles A vector containing file paths of 1 or more BAM files
#' @param sample_names The sample names of the given BAM files. Must
#'   be a vector of the same length as `bamfiles`
#' @param reference_path The directory of the NxtIRF reference
#' @param output_path The directory where NxtIRF/IRFinder output
#'   should be stored
#' @param n_threads The number of threads to use. On Linux / Windows, this will
#'   use OpenMP from within the C++ subroutine. On Macs, BiocParallel
#'   MulticoreParam will be used on single-threaded NxtIRF/IRFinder
#' @param Use_OpenMP (default TRUE) Whether to use OpenMP to run IRFinder.
#'   Typically Windows and Linux will natively support OpenMP, whereas MacOS
#'   will not by default. If OpenMP is not available or `Use_OpenMP = FALSE`,
#'   then IRFinder multi-threading will use BiocParallel (SerialParam() for
#'   Windows, and MulticoreParam for Linux / MacOS).
#' @param overwrite (default FALSE) If IRFinder output files already exist,
#'   will not attempt to re-run.
#' @param run_featureCounts Whether this function will run 
#'   `Rsubread::featureCounts()` on the BAM files. If so, the output will be
#'   saved to "main.FC.Rds" in the output directory as a list object
#' @param verbose (default FALSE) Set to `TRUE` to allow IRFinder to output
#'   progress bars and messages
#' @return None. `IRFinder()` will save output to `output_path`. \cr\cr
#'   sample.txt.gz: The main IRFinder output file containing the quantitation
#'   of IR and splice junctions, as well as QC information\cr\cr
#'   sample.cov: Contains coverage information in compressed binary. This
#'   format is 5-10X faster than BigWig format (see [GetCoverage()])\cr\cr
#'   main.FC.Rds: A single file containing gene counts for the whole dataset
#'   (only if `run_featureCounts == TRUE`)
#' @examples
#' bams = NxtIRF_example_bams()
#' IRFinder(bams$path, bams$sample,
#'   reference_path = file.path(tempdir(), "Reference"),
#'   output_path = file.path(tempdir(), "IRFinder_output")
#' )
#' @md
#' @export
IRFinder <- function(
        bamfiles = "./Unsorted.bam", 
        sample_names = "sample1",
        reference_path = "./Reference",
        output_path = "./IRFinder_Output",
        n_threads = 1, Use_OpenMP = TRUE,
        overwrite = FALSE,
        run_featureCounts = FALSE,
        verbose = FALSE
) {
    # Check args
    if(length(bamfiles) != length(sample_names)) .log(paste("In IRFinder(),",
        "Number of BAM files and sample names must be the same"))
    if(!all(file.exists(bamfiles))) .log(paste("In IRFinder(),",
        "some BAMs in bamfiles do not exist"))
    if(!dir.exists(dirname(output_path))) .log(paste("In IRFinder(),",
        dirname(output_path), " - path does not exist"))
    if(!dir.exists(output_path)) dir.create(output_path)

    # Check which output already exists; prevent overwrite
    s_output = file.path(normalizePath(output_path), sample_names)
    if(!overwrite) {
        already_exist = (
            file.exists(paste(s_output, "txt.gz")) &
            file.exists(paste(s_output, ".cov"))
        )
    } else {
        already_exist = rep(FALSE, length(bamfiles))
    }

    # Call wrapper
    .run_IRFinder(
        reference_path = reference_path,
        bamfiles = bamfiles[!already_exist],
        output_files = s_output[!already_exist],
        max_threads = n_threads, Use_OpenMP = Use_OpenMP,
        run_featureCounts = run_featureCounts,
        overwrite_IRFinder_output = overwrite,
        verbose = verbose
    )
}

# IRFinder wrapper to R/C++. Handles whether OpenMP or BiocParallel is used
.run_IRFinder = function(
        reference_path = "./Reference", 
        bamfiles = "Unsorted.bam", 
        output_files = "./Sample",
        max_threads = max(parallel::detectCores() - 2, 1),
        Use_OpenMP = TRUE,
        run_featureCounts = FALSE,
        overwrite_IRFinder_output = FALSE,
        verbose = TRUE
    ) {
    .validate_reference(reference_path)     # Check valid NxtIRF reference
    s_bam = normalizePath(bamfiles)         # Clean path name for C/IRFinder
    s_ref = normalizePath(reference_path)   # Clean path name for C/IRFinder
    
    # Check args
    .irfinder_validate_args(s_bam, s_ref, max_threads, output_files)
    ref_file = normalizePath(file.path(s_ref, "IRFinder.ref.gz"))

    message("Running IRFinder ", appendLF = FALSE)
    n_threads = floor(max_threads)
    if(Has_OpenMP() > 0 & Use_OpenMP) {
        IRF_main_multi(ref_file, s_bam, output_files, n_threads, verbose)
    } else {
        # Use BiocParallel
        n_rounds = ceiling(length(s_bam) / floor(max_threads))
        n_threads = ceiling(length(s_bam) / n_rounds)

        BPPARAM_mod = .validate_threads(n_threads, as_BPPARAM = TRUE)

        row_starts = seq(1, by = n_threads, length.out = n_rounds)
        for(i in seq_len(n_rounds)) {
            selected_rows_subset = seq(row_starts[i], 
                min(length(s_bam), row_starts[i] + n_threads - 1)
            )
            BiocParallel::bplapply(selected_rows_subset,
                function(i, s_bam, reference_file, 
                        output_files, verbose, overwrite) {
                    .irfinder_run_single(s_bam[i], reference_file, 
                        output_files[i], verbose, overwrite)
                }, 
                s_bam = s_bam,
                reference_file = ref_file,
                output_files = output_files,
                verbose = verbose,
                overwrite = overwrite_IRFinder_output,
                BPPARAM = BPPARAM_mod
            )
        }
    }
    if(run_featureCounts == TRUE) {
        .irfinder_run_featureCounts(reference_path, output_files, 
            s_bam, n_threads)
    }
}

# Call C++/IRFinder on a single sample. Used for BiocParallel
.irfinder_run_single <- function(
    bam, ref, out, verbose, overwrite, n_threads = 1
) {
    file_gz = paste0(out, ".txt.gz")
    file_cov = paste0(out, ".cov")
    bam_short = file.path(basename(dirname(bam)), basename(bam))
    if(overwrite ||
        !(file.exists(file_gz) | file.exists(file_cov))) {
        ret = IRF_main(bam, ref, out, verbose, n_threads)
        # Check IRFinder returns all files successfully
        if(ret != 0) {
            .log(paste(
                "IRFinder exited with errors, see error messages above"))
        } else if(!file.exists(file_gz)) {
            .log(paste(
                "IRFinder failed to produce", file_gz))
        } else if(!file.exists(file_cov)) {
            .log(paste(
                "IRFinder failed to produce", file_cov))
        } else {
            .log(paste("IRFinder processed", bam_short), "message")
        }
    } else {
        .log(paste("IRFinder output for", bam_short, 
            "already exists, skipping..."), "message")
    }
}

# Runs featureCounts on given BAM files, intended to be run after IRFinder
# as the IRFinder determines the strandedness and paired-ness of the experiment
.irfinder_run_featureCounts <- function(reference_path, output_files, 
        s_bam, n_threads) {
    NxtIRF.CheckPackageInstalled("Rsubread", "2.4.0")
    gtf_file <- Get_GTF_file(reference_path)
    
    # determine paired-ness, strandedness, assume all BAMS are the same
    data.list = get_multi_DT_from_gz(
        normalizePath(paste0(output_files[1], ".txt.gz")), 
        c("BAM", "Directionality")
    )
    stats = data.list$BAM
    direct = data.list$Directionality

    paired = (stats$Value[3] == 0 & stats$Value[4] > 0) || 
        (stats$Value[3] > 0 && stats$Value[4] / stats$Value[3] / 1000)
    strand = direct$Value[9]
    if(strand == -1) strand = 2
    
    # Run FeatureCounts in bulk
    res = Rsubread::featureCounts(
        s_bam,
        annot.ext = gtf_file,
        isGTFAnnotationFile = TRUE,
        strandSpecific = strand,
        isPairedEnd = paired,
        requireBothEndsMapped = paired,
        nthreads = n_threads
    )
    
    # Append to existing main.FC.Rds if exists:
    outfile = file.path(dirname(output_files[1]), "main.FC.Rds")
    if(file.exists(outfile)) {
        res.old = readRDS(outfile)

        # Check annotation df and Status vectors are identical
        anno.old = res.old$annotation[, 
            c("GeneID", "Chr", "Start", "End", "Strand")]
        anno.new = res$annotation[, 
            c("GeneID", "Chr", "Start", "End", "Strand")]
        if(
            identical(anno.old, anno.new) & 
            identical(res.old$stat$Status, res$stat$Status)
        ) {
            new_samples = res$targets[!(res$targets %in% res.old$targets)]
            res$targets = c(res.old$targets, new_samples)
            res$stat = cbind(res.old$stat, res$stat[,new_samples])        
            res$counts = cbind(res.old$counts, res$counts[,new_samples])
        } else {
            .log(paste(
                "featureCounts output not compatible with previous output",
                "in main.FC.Rds. Overwriting previous output"
            ), "warning")
        }
    }
    if(all(c("counts", "annotation", "targets", "stat") %in% names(res))) {
        saveRDS(res, outfile)
    }
    .log(paste("featureCounts ran succesfully; saved to", outfile), "message")
}

# Validate arguments; return error if invalid
.irfinder_validate_args <- function(s_bam, s_ref, max_threads, output_files) {
    if(!is.numeric(max_threads)) max_threads = 1
    if(max_threads < 1) max_threads = 1
    max_threads = floor(max_threads)
    
    if(max_threads > 1 && max_threads > parallel::detectCores()) {
        .log(paste("In .run_IRFinder(), ",
            max_threads, " threads is not allowed for this system"))
    }
    
    if(!all(file.exists(s_bam))) {
        .log(paste("In .run_IRFinder(), ",
            paste(unique(s_bam[!file.exists(s_bam)]), collapse = ""),
            " - these BAM files were not found"))
    }    

    if(!all(dir.exists(dirname(output_files)))) {
        .log(paste("In .run_IRFinder(), ",
            paste(unique(dirname(
                    output_files[!dir.exists(dirname(output_files))])),
                collapse = ""),
            " - directories not found"))
    }

    if(!(length(s_bam) == length(output_files))) {
        .log(paste("In .run_IRFinder(), ",
            "Number of output files and bam files must be the same"))
    }
    return(TRUE)
}
