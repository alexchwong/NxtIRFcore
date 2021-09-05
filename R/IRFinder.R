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
        n_threads = 1,
        overwrite = FALSE,
        run_featureCounts = FALSE,
        verbose = FALSE
        ) {
    if(length(bamfiles) != length(sample_names)) {
        .log(paste("In IRFinder,",
            "Number of BAM files and sample names must be the same"))
    }
    if(!all(file.exists(bamfiles))) {
        .log(paste("In IRFinder,",
            "some BAMs in bamfiles do not exist"))
    }
    if(!dir.exists(dirname(output_path))) {
        .log(paste("In IRFinder,",
            dirname(output_path), " - path does not exist"))
    }
    if(!dir.exists(output_path)) dir.create(output_path)

    s_output = file.path(normalizePath(output_path), sample_names)

    if(!overwrite) {
        already_exist = (
            file.exists(paste(s_output, "txt.gz")) &
            file.exists(paste(s_output, ".cov"))
        )
    } else {
        already_exist = rep(FALSE, length(bamfiles))
    }

    .run_IRFinder(
        reference_path = reference_path,
        bamfiles = bamfiles[!already_exist],
        output_files = s_output[!already_exist],
        max_threads = n_threads,
        run_featureCounts = run_featureCounts,
        overwrite_IRFinder_output = overwrite,
        verbose = verbose
    )
}
