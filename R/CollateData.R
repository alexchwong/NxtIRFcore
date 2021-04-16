#' Searches a specified path for files of a particular file pattern
#'
#' This convenience function identifies files with the specified suffix.
#' The sample name is assumed to be the name of the file minus its suffix
#' (`use_subdir = FALSE`), but can also be changed to be the names 
#' of the parent directory (`use_subdir = TRUE`). See example below.
#'
#' @param sample_path The path in which to recursively search for files
#'   that match the given `suffix`
#' @param suffix A vector of or or more strings that specifies the file suffix 
#'   (e.g. '.bam' denotes BAM files, whereas ".txt.gz" denotes gzipped txt 
#'   files).
#' @param suffix_type A vector of string that determines the column
#'   names of the files retrieved by `suffix`. Must be the same length as 
#'   `suffix`
#' @param use_subdir Whether to assume the directory name containing
#'   the found files denote the sample name. If `FALSE`, the base name
#'   of the file is assumed to be the sample name. See below example.
#' @return A 2-column data frame with the first column containing
#'   the sample name, and the second column being the file path.
#' @examples
#' # Return all BAM files using file names as sample names
#' bam_path = dirname(example_bams())[1]
#' df = FindSamples(sample_path = bam_path, 
#'   suffix = ".bam", suffix_type = "bam_file", use_subdir = FALSE)
#' @md
#' @export
FindSamples <- function(sample_path, suffix = ".txt.gz", 
            suffix_type = "path", use_subdir = FALSE) {
    if(!dir.exists(sample_path)) {
        stop(paste("In FindSamples(),",
            sample_path, "- given path does not exist"
        ), call. = FALSE)
    }
    if(length(suffix) == 0 || length(suffix) != length(suffix_type)) {
        stop(paste("In FindSamples(),",
            "suffix must be of length greater than zero",
            "and same length as suffix_type"
        ), call. = FALSE)
    }
    DT.list = list()
    for(i in seq_len(length(suffix))) {
        files_found = list.files(pattern = paste0("\\", suffix[i], "$"),
        path = normalizePath(sample_path), full.names = TRUE, recursive = TRUE)
        if(length(files_found) > 0) {
            DT = data.table(sample = "", path = files_found)
            if(use_subdir) {
                DT$sample = basename(dirname(DT$path))
            } else {
                DT$sample = sub(suffix[i],"",basename(DT$path))
            }
            colnames(DT)[2] = suffix_type[i]
            DT.list[[i]] = DT
        } else {
            DT.list[[i]] = NULL
        }
    }
    if(length(DT.list) <= 1) return(as.data.frame(DT.list))
    final = c()
    if(length(suffix) > 1) {
        for(i in seq_len(length(suffix))) {
            if(!is.null(DT.list[[i]])) {
                if(is.null(final)) {
                    final = DT.list[[i]]
                } else {
                    final = DT.list[[i]][final, on = "sample"]
                }
            }
        }
    }
    cols = c("sample", suffix_type[suffix_type %in% colnames(final)])
    return(as.data.frame(final[, cols, with = FALSE]))
}

#' Convenience Function to find BAM files in a certain folder
#' 
#' Runs FindSamples to find BAM files. Assumes file names are names of samples.
#' 
#' @param sample_path The path in which to recursively search for BAM files
#' @param ... Additional parameters to pass into `FindSamples()`
#' @return A 2-column data frame with the first column containing
#'   the sample name, and the second column being the BAM file path.
#' @seealso [FindSamples]
#' @examples
#' bam_path = dirname(example_bams())[1]
#' bams = Find_Bams(bam_path)
#' @md
#' @export
Find_Bams <- function(sample_path, ...) {
    return(FindSamples(sample_path, ".bam", "BAM", ...))
}

#' Convenience Function to find IRFinder output and COV files
#' 
#' Runs FindSamples to find IRFinder .txt.gz and .cov files. 
#' Assumes file names are names of samples.
#' 
#' @param sample_path The path in which to recursively search for BAM files
#' @param ... Additional parameters to pass into `FindSamples()`
#' @return A 3-column data frame with the first column containing
#'   the sample name, and the second column being the IRFinder main output path,
#'   and the third column being the COV file path.
#' @seealso [FindSamples]
#' @examples
#' expr = Find_IRFinder_Output(file.path(tempdir(), "IRFinder_output"))
#' @md
#' @export
Find_IRFinder_Output <- function(sample_path, ...) {
    return(FindSamples(sample_path, 
        c(".txt.gz", ".cov"), 
        c("irf_file", "cov_file"), 
        ...)
    )
}

#' A wrapper function to call NxtIRF/IRFinder
#'
#' This function calls IRFinder on one or more BAM files.
#' @param bamfiles The file names of 1 or more BAM files
#' @param sample_names The sample names of the given BAM files. Must
#'   be a vector of the same length as `bamfiles`
#' @param reference_path The directory of the NxtIRF reference
#' @param output_path The directory where NxtIRF/IRFinder output
#'   should be stored
#' @param n_threads The number of threads to use. On Linux / Windows, this will
#'   use OpenMP from within the C++ subroutine. On Macs, BiocParallel
#'   MulticoreParam will be used on single-threaded NxtIRF/IRFinder
#' @param run_featureCounts Whether this function will run 
#'   `Rsubread::featureCounts()` on the BAM files. If so, the output will be
#'   saved to "main.FC.Rds" in the output directory as a list object
#' @return None. `IRFinder()` will save output to `output_path`. \cr\cr
#'   sample.txt.gz: The main IRFinder output file containing the quantitation
#'   of IR and splice junctions, as well as QC information\cr\cr
#'   sample.cov: Contains coverage information in compressed binary. This
#'   format is 5-10X faster than BigWig format (see [GetCoverage()])\cr\cr
#'   main.FC.Rds: A single file containing gene counts for the whole dataset
#'   (only if `run_featureCounts == TRUE`)
#' @examples
#' bams = NxtIRF_example_bams()
#' IRFinder(bams$BAM, bams$sample,
#'   reference_path = file.path(tempdir(), "Reference"),
#'   output_path = file.path(tempdir(), "IRFinder_output")
#' )
#' @md
#' @export
IRFinder <- function(
        bamfiles = "Unsorted.bam", 
        sample_names = "sample1",
        reference_path = "./Reference",
        output_path = "./IRFinder_Output",
        n_threads = 1,
        run_featureCounts = FALSE
        ) {
    if(length(bamfiles) != length(sample_names)) {
        stop(paste("In IRFinder,",
            "Number of BAM files and sample names must be the same"
        ), call. = FALSE)
    }
    if(!dir.exists(dirname(output_path))) {
        stop(paste("In IRFinder,",
            dirname(output_path), " - path does not exist"
        ), call. = FALSE)
    }
    if(!dir.exists(output_path)) dir.create(output_path)

    s_output = file.path(normalizePath(output_path), sample_names)
    run_IRFinder_multithreaded(
        reference_path = reference_path,
        bamfiles = bamfiles,
        output_files = s_output,
        max_threads = n_threads,
        run_featureCounts = run_featureCounts
    )
}

#' Processes data from IRFinder output
#'
#' CollateData unifies a list of IRFinder output files belonging to an 
#' experiment. It is assumed every sample is analysed using the same IRFinder
#' reference. The combination of junction counts and IR quantification from
#' IRFinder is used to calculate percentage spliced in (PSI) of alternative
#' splice events, and percent intron retention (PIR) of retained introns. Also,
#' QC information is extracted, and data is collated into fst files for fast
#' downstream access such as \code{MakeSE()}.
#'
#' @param Experiment A 2-column data frame (generated by \code{FindSamples()}),
#'   with the first column designating the sample names, and the 2nd column 
#'   containing the primary IRFinder output file (of type \code{.txt.gz}). A 
#'   third optional column can contain the coverage files of the corresponding 
#'   samples. NB: all other columns are ignored.\cr\cr
#'   We recommend using the output of \code{Find_IRFinder_Output()} here.
#' @param reference_path The path to the reference generated by BuildReference()
#' @param output_path The path for the output files to be generated by this
#'   function.
#' @param IRMode The algorithm to calculate 'splice abundance' in IR 
#'   quantification. The original algorithm by Middleton et al (2017) proposes
#'   \code{SpliceMax}, which calculates the number of mapped splice events
#'   that share the boundary coordinate of either the left or right flanking
#'   exon (SpliceLeft, SpliceRight) and defines splice abundance as the larger
#'   of the two values. NxtIRF proposes a new algorithm, \code{SpliceOverMax},
#'   to account for the possibility that the major isoform shares neither
#'   boundary, but arises from either of the flanking "exon islands". Exon
#'   islands are contiguous regions covered by exons from any transcript
#'   (except those designated as \code{retained_intron} or 
#'   \code{sense_intronic}), and are separated by
#'   obligate intronic regions (genomic regions that are introns for all 
#'   transcripts). Note for introns that are internal to a single exon island
#'   (i.e. akin to "known-exon" introns from IRFinder), \code{SpliceOverMax} 
#'   uses \code{GenomicRanges::findOverlaps()} to summate competing mapped
#'   splice reads.
#' @param low_memory_mode Use this mode in memory-limited systems with many
#'   samples (> 16). CollateData will write to file for every N samples as
#'   defined by \code{samples_per_block = N}. Memory usage is often a problem
#'   for large datasets or using multiple cores. If you experience crashes due
#'   to running out of memory, set this to true and make sure 
#'   \code{n_threads = 1}
#' @param samples_per_block How many samples to process per thread. Use in 
#'   conjunction with low_memory_mode to lower memory requirements
#' @param n_threads The number of threads to use.
#' @return None. \code{CollateData()} writes to the directory given by 
#'   \code{output_path}
#' @examples
#' bams = NxtIRF_example_bams()
#' IRFinder(bams$BAM, bams$sample,
#'   reference_path = file.path(tempdir(), "Reference"),
#'   output_path = file.path(tempdir(), "IRFinder_output")
#' )
#' expr = Find_IRFinder_Output(file.path(tempdir(), "IRFinder_output"))
#' CollateData(expr, 
#'   reference_path = file.path(tempdir(), "Reference"),
#'   output_path = file.path(tempdir(), "NxtIRF_output")
#' )
#' @seealso "NxtIRF: 2 - Building an Experiment" 
#'   (accessed via `browseVignettes("NxtIRF")`)
#' @export
CollateData <- function(Experiment, reference_path, output_path,
        IRMode = c("SpliceOverMax", "SpliceMax"), 
        low_memory_mode = FALSE, samples_per_block = 16, n_threads = 1) {

    IRMode = match.arg(IRMode)
    if(IRMode == "") {
        stop(paste("In CollateData(),",
            "IRMode must be either 'SpliceOverMax' (default) or 'SpliceMax'"
        ), call. = FALSE)
    }

    BPPARAM_mod = .validate_threads(n_threads)
    norm_output_path = .collateData_validate(Experiment, 
        reference_path, output_path)   
    coverage_files = .collateData_COV(Experiment)
    
    df.internal <- .collateData_expr(Experiment)
    # jobs <- .collateData_jobs(nrow(df.internal), BPPARAM_mod, samples_per_block)
    jobs = NxtIRF.SplitVector(seq_len(nrow(df.internal)),
        ceiling(nrow(df.internal) / samples_per_block))
    n_jobs = length(jobs)
    N <- 7
    dash_progress("Compiling Sample Stats", N)
    message("Compiling Sample Stats")
    df.internal = .collateData_stats(df.internal, jobs, BPPARAM_mod)
    if(any(df.internal$strand == 0)) {
        runStranded = FALSE
    } else {
        runStranded = TRUE
    }    
    
    dash_progress("Compiling Junction List", N)
    message("Compiling Junction List")       
    junc.common <- .collateData_junc_merge(df.internal, jobs, BPPARAM_mod, 
        output_path)
    gc()

    dash_progress("Compiling Intron Retention List", N)
    message("Compiling Intron Retention List")
    irf.common <- .collateData_irf_merge(df.internal, jobs, BPPARAM_mod, 
        output_path, runStranded)
    gc()

    irf.common[, start := start + 1]
    junc.common[, start := start + 1]

# Reassign +/- based on junctions.fst annotation
    # Annotate junctions
    dash_progress("Tidying up splice junctions and intron retentions", N)
    message("Tidying up splice junctions and intron retentions...")
    .collateData_annotate(reference_path, norm_output_path,
        junc.common, irf.common, runStranded)
    message("done\n")
    gc()

    dash_progress("Generating NxtIRF FST files", N)
    message("Generating NxtIRF FST files")
    
    item.todo = c("Included", "Excluded", "Depth", "Coverage", "minDepth", 
        "Up_Inc", "Down_Inc", "Up_Exc", "Down_Exc", "junc_PSI", "junc_counts")

    agg.list <- suppressWarnings(BiocParallel::bplapply(seq_len(n_jobs),
        .collateData_compile_agglist, 
        jobs = jobs, df.internal = df.internal, 
        norm_output_path = norm_output_path, IRMode = IRMode,
        low_memory_mode = low_memory_mode,
        BPPARAM = BPPARAM_mod
    ))
    gc()

    dash_progress("Building Final SummarizedExperiment Object", N)
    message("Building Final SummarizedExperiment Object")
    assays <- .collateData_compile_assays(agg.list, df.internal,
        norm_output_path, low_memory_mode,
        item.todo, jobs, n_jobs)
    # .collateData_fix_juncnames(norm_output_path)
    .collateData_write_stats(df.internal, norm_output_path)
    .collateData_write_colData(df.internal, coverage_files, norm_output_path)
    cov_data <- prepare_covplot_data(reference_path)
    saveRDS(cov_data, file.path(norm_output_path, "Annotation", "cov_data.Rds"))
    
    # NEW compile HDF5SummarizedExperiment:
    colData.Rds <- readRDS(file.path(norm_output_path, "colData.Rds"))
    colData <- .makeSE_colData_clean(
        colData.Rds$df.anno)
    se <- .makeSE_initialise_HDF5(norm_output_path, colData)
    
    HDF5Array::saveHDF5SummarizedExperiment(se, 
        dir = norm_output_path, prefix = "NxtSE_")
    
    dash_progress("NxtIRF Collation Finished", N)
    message("NxtIRF Collation Finished")
}

################################################################################
# CollateData helper functions:

.collateData_validate <- function(Experiment, reference_path, output_path) {
    if(!is(Experiment, "data.frame")) {
        stop(paste("In CollateData(),",
            "Experiment object needs to be a data frame"
        ), call. = FALSE)
    }
    if(ncol(Experiment) < 2) {
        stop(paste("In CollateData(),",
            "Experiment needs to contain at least two columns,",
            "with the first 2 columns containing",
            "(1) sample name and (2) IRFinder output"
        ), call. = FALSE)
    }
    .validate_reference(reference_path)
    if(!dir.exists(dirname(output_path))) {
        stop(paste("In CollateData(),",
            "Parent directory of output path:", 
            dirname(output_path), "needs to exist"
        ), call. = FALSE)
    }
    base_output_path = normalizePath(dirname(output_path)) 
    norm_output_path = file.path(base_output_path, basename(output_path))
    if(!dir.exists(norm_output_path)) {
        dir.create(norm_output_path)
    }
    temp_output_path = file.path(norm_output_path, "temp")
    if(!dir.exists(temp_output_path)) {
        dir.create(temp_output_path)
    }        
    if(!dir.exists(file.path(norm_output_path, "samples"))) {
        dir.create(file.path(norm_output_path, "samples"))
    }
    if(!dir.exists(file.path(norm_output_path, "annotation"))) {
        dir.create(file.path(norm_output_path, "annotation"))
    }
    return(norm_output_path)
}

.collateData_COV <- function(Experiment) {
    coverage_files = ""
    Experiment = as.data.frame(Experiment)
    if(ncol(Experiment) > 2 && all(file.exists(Experiment[,3]))) {
        coverage_files = Experiment[,3]
    }
    if(!IsCOV(coverage_files)){
        message("Some coverage files do not exist or are corrupted")
        coverage_files = ""
    }
    if(length(coverage_files) != nrow(Experiment)) {
        message("The number of coverage files must equal the number of samples")
        coverage_files = ""
    }
    return(coverage_files)
}

.collateData_expr <- function(Experiment) {
    Experiment = as.data.frame(Experiment)
    colnames(Experiment)[c(1,2)] = c("sample", "path")
    if(!all(vapply(Experiment$path, file.exists, logical(1)))) {
        stop(paste("In CollateData(),",
            "Some files in Experiment do not exist"
        ), call. = FALSE)
    }

    df.internal = as.data.table(Experiment[,c(1,2)])
    df.internal$paired = FALSE
    df.internal$strand = 0
    df.internal$depth = 0
    df.internal$mean_frag_size = 0
    df.internal$directionality_strength = 0
    df.internal$Intergenic_Fraction = 0
    df.internal$rRNA_Fraction = 0
    df.internal$NonPolyA_Fraction = 0
    df.internal$Mitochondrial_Fraction = 0
    df.internal$Unanno_Jn_Fraction = 0
    df.internal$NMD_Jn_Fraction = 0
    df.internal$Fraction_Splice_Reads = 0
    df.internal$Fraction_Span_Reads = 0

    df.internal$IRBurden_clean = 0
    df.internal$IRBurden_exitrons = 0
    df.internal$IRBurden_clean_unstranded = 0
    df.internal$IRBurden_exitrons_unstranded = 0
    df.internal$IRBurden_antisense = 0
    return(df.internal)
}

.collateData_jobs <- function(n_expr, BPPARAM_mod, samples_per_block) {
    n_jobs = min(ceiling(n_expr / samples_per_block), 
        BPPARAM_mod$workers)
    jobs = NxtIRF.SplitVector(seq_len(n_expr), n_jobs)  
    n_jobs = length(jobs)
    return(jobs)
}



################################################################################
# Sub

.collateData_stats <- function(df.internal, jobs, BPPARAM_mod) {
    n_jobs = length(jobs)
    df.internal = suppressWarnings(rbindlist(
        BiocParallel::bplapply(
            seq_len(n_jobs),
            function(x, jobs, df.internal) {
                suppressPackageStartupMessages({
                    requireNamespace("data.table")
                    requireNamespace("stats")
                })
                work = jobs[[x]]
                block = df.internal[work]
                for(i in seq_len(length(work))) {
                    if(file.exists(paste0(block$path[i], ".sqlite"))) {
                        sqldb <- .open_sqlite(paste0(block$path[i], ".sqlite"))
                        stats = .from_sqlite(sqldb, "stats", DT = TRUE)
                        direct = .from_sqlite(sqldb, "direction", DT = TRUE)
                        QC = .from_sqlite(sqldb, "QC", DT = TRUE)
                        .close_sqlite(sqldb)
                    } else {
                        data.list = get_multi_DT_from_gz(
                            normalizePath(block$path[i]), 
                            c("BAM", "Directionality", "QC")) 
                        stats = data.list$BAM
                        direct = data.list$Directionality
                        QC = data.list$QC                    
                    }
                    block <- .collateData_stats_reads(block, i, stats, direct)
                    block <- .collateData_stats_QC(block, i, QC, direct)
                }
                return(block)
            }, jobs = jobs, df.internal = df.internal, BPPARAM = BPPARAM_mod
        )
    ))
    return(df.internal)
}

################################################################################

.collateData_stats_reads <- function(block, i, stats, direct) {
    if(stats$Value[3] == 0 & stats$Value[4] > 0) {
        block$paired[i] = TRUE
        block$depth[i] = stats$Value[4]
        block$mean_frag_size[i] = stats$Value[2] / 
            stats$Value[4]
    } else if(stats$Value[3] > 0 && 
            stats$Value[4] / stats$Value[3] / 1000) {
        block$paired[i] = TRUE
        block$depth[i] = stats$Value[4]
        block$mean_frag_size[i] = stats$Value[2] / 
            stats$Value[4]
    } else {
        block$paired[i] = FALSE
        block$depth[i] = stats$Value[3]
        block$mean_frag_size[i] = stats$Value[2] / 
            stats$Value[3]
    }
    block$strand[i] = direct$Value[9]
    return(block)
}

.collateData_stats_QC <- function(block, i, QC, direct) {
    block$directionality_strength[i] = direct$Value[8]
    block$Intergenic_Fraction[i] =
        QC$Value[QC$QC == "Intergenic Reads"] / 
            block$depth[i]
    block$rRNA_Fraction[i] =    
        QC$Value[QC$QC == "rRNA Reads"] / 
            block$depth[i]
    block$NonPolyA_Fraction[i] =
        QC$Value[QC$QC == "NonPolyA Reads"] / 
            block$depth[i]
    block$Mitochondrial_Fraction[i] =
        QC$Value[QC$QC == "Mitochondrial Reads"] / 
            block$depth[i]
    block$Unanno_Jn_Fraction[i] =
        QC$Value[QC$QC == "Unannotated Junctions"] / 
        (QC$Value[QC$QC == "Unannotated Junctions"] +
        QC$Value[QC$QC == "Annotated Junctions"])
    block$NMD_Jn_Fraction[i] =
        QC$Value[QC$QC == "NMD Junctions"] / 
        QC$Value[QC$QC == "Annotated Junctions"]
    block$Fraction_Splice_Reads[i] =
        QC$Value[QC$QC == "Annotated Junctions"] / 
        block$depth[i]
    block$Fraction_Span_Reads[i] =
        QC$Value[QC$QC == "Spans Reads"] / 
            block$depth[i]
# IRBurden calculations                            
    block$IRBurden_clean_unstranded[i] =
        QC$Value[QC$QC == "Non-Directional Clean IntronDepth Sum"] / 
        (QC$Value[QC$QC == "Non-Directional Clean IntronDepth Sum"] +
        QC$Value[QC$QC == "Annotated Junctions"])
    block$IRBurden_exitrons_unstranded[i] =
        QC$Value[QC$QC == "Non-Directional Known-Exon IntronDepth Sum"] / 
        (QC$Value[QC$QC == "Non-Directional Known-Exon IntronDepth Sum"] +
        QC$Value[QC$QC == "Annotated Junctions"])
    block$IRBurden_antisense[i] =
        QC$Value[QC$QC == "Non-Directional Anti-Sense IntronDepth Sum"] / 
        (QC$Value[QC$QC == "Non-Directional Anti-Sense IntronDepth Sum"] +
        QC$Value[QC$QC == "Annotated Junctions"])
    if(block$strand[i] != 0) {
        block$IRBurden_clean[i] =
            QC$Value[QC$QC == "Directional Clean IntronDepth Sum"] / 
            (QC$Value[QC$QC == "Directional Clean IntronDepth Sum"] +
            QC$Value[QC$QC == "Annotated Junctions"])
        block$IRBurden_exitrons[i] =
            QC$Value[QC$QC == "Directional Known-Exon IntronDepth Sum"] / 
            (QC$Value[QC$QC == "Directional Known-Exon IntronDepth Sum"] +
            QC$Value[QC$QC == "Annotated Junctions"])
    }
    return(block)    
}

################################################################################
# Sub

.collateData_junc_merge <- function(df.internal, jobs, BPPARAM_mod, 
        output_path) {
    temp_output_path = file.path(output_path, "temp")
    n_jobs = length(jobs)
    # Compile junc.common via merge
    junc.list = suppressWarnings(BiocParallel::bplapply(
        seq_len(n_jobs),
        function(x, jobs, df.internal, temp_output_path) {
            suppressPackageStartupMessages({
                requireNamespace("data.table")
                requireNamespace("stats")
            })
            work = jobs[[x]]
            block = df.internal[work]
            junc.segment = NULL
            for(i in seq_len(length(work))) {
                if(file.exists(paste0(block$path[i], ".sqlite"))) {
                    sqldb <- .open_sqlite(paste0(block$path[i], ".sqlite"))
                    junc <- .from_sqlite(sqldb, "JC", DT = TRUE)
                    .close_sqlite(sqldb)
                } else {
                    junc = suppressWarnings(as.data.table(
                        fread(block$path[i], skip = "JC_seqname")))
                }
                setnames(junc, "JC_seqname", "seqnames")
                if(is.null(junc.segment)) {
                    junc.segment = junc[,seq_len(4), with = FALSE]
                } else {
                    junc.segment = merge(junc.segment, 
                        junc[,seq_len(4), with = FALSE], 
                        all = TRUE)
                }
                # Write temp file
                fst::write.fst(as.data.frame(junc), 
                    file.path(temp_output_path, paste(block$sample[i], 
                    "junc.fst.tmp", sep=".")))        
            }
            return(junc.segment)
        }, jobs = jobs, df.internal = df.internal, 
            temp_output_path = temp_output_path, BPPARAM = BPPARAM_mod
    ))
    junc.common = NULL
    for(i in seq_len(length(junc.list))) {
        if(is.null(junc.common)) {
            junc.common = junc.list[[i]]
        } else {
            junc.common = merge(junc.common, junc.list[[i]], 
                all = TRUE, by = colnames(junc.common))
        }
    }
    return(junc.common)
}

.collateData_irf_merge <- function(df.internal, jobs, BPPARAM_mod, 
        output_path, runStranded) {
    # Check MD5 of all introns are the same in all samples
    temp_output_path = file.path(output_path, "temp")
    n_jobs = length(jobs)
    irf.list = suppressWarnings(BiocParallel::bplapply(seq_len(n_jobs),
        function(x, jobs, df.internal, temp_output_path, runStranded) {
            suppressPackageStartupMessages({
                requireNamespace("data.table")
                requireNamespace("stats")
                requireNamespace("openssl")
            })
            work = jobs[[x]]
            block = df.internal[work]
            irf.md5 = c()
            phrase = ifelse(runStranded, "Dir_Chr", "Nondir_Chr")
            for(i in seq_len(length(work))) {
                if(file.exists(paste0(block$path[i], ".sqlite"))) {
                    sqldb <- .open_sqlite(paste0(block$path[i], ".sqlite"))
                    if(runStranded) {
                        irf = .from_sqlite(sqldb, "Dir", DT = TRUE)
                    } else {
                        irf = .from_sqlite(sqldb, "Nondir", DT = TRUE)
                    }
                    .close_sqlite(sqldb)
                } else {
                    irf = suppressWarnings(
                        fread(block$path[i], skip = phrase)
                    )
                }
                setnames(irf, c(phrase, "Start", "End", "Strand"), 
                    c("seqnames","start","end", "strand"))
                irf.md5 = unique(c(irf.md5, 
                    as.character(openssl::md5(paste(irf$Name, collapse=" ")))))
                fst::write.fst(as.data.frame(irf), 
                    file.path(temp_output_path, 
                        paste(block$sample[i], "irf.fst.tmp", sep=".")))
            }
            return(irf.md5)
        },  jobs = jobs, df.internal = df.internal, 
            temp_output_path = temp_output_path, 
            runStranded = runStranded, BPPARAM = BPPARAM_mod
    ))
    irf.md5.check = unique(unlist(irf.list))
    if(length(irf.md5.check) > 1) {
        stop(paste(
            "MD5 check of IRFinder introns are not the same.",
            "Perhaps some samples were processed by a different reference.",
            "NxtIRF needs all samples to be processed by the same reference"
        ), call. = FALSE)
    }
    irf = fst::read.fst(file.path(temp_output_path, 
        paste(df.internal$sample[1], "irf.fst.tmp", sep=".")),
        as.data.table = TRUE)
    irf.common = irf[,seq_len(6), with = FALSE]
    return(irf.common)
}

################################################################################
# Sub

.collateData_annotate <- function(reference_path, norm_output_path,
        junc.common, irf.common, runStranded) {
    
    message("...annotating splice junctions")
    junc.common = .collateData_junc_annotate(junc.common, reference_path)
    message("...grouping splice junctions")
    junc.common = .collateData_junc_group(junc.common, reference_path)
    irf.common = .collateData_irf_group(irf.common, reference_path, runStranded)
    irf.common[, c("EventRegion") := 
        paste0(get("seqnames"), ":", get("start"), "-", 
            get("end"), "/", get("strand"))]
    message("...loading splice events")
    Splice.Anno = .collateData_splice_anno(reference_path)
    message("...saving annotations")
    # Save annotation
    write.fst(as.data.frame(junc.common), 
        file.path(norm_output_path, "annotation", "Junc.fst"))
    write.fst(as.data.frame(irf.common), 
        file.path(norm_output_path, "annotation", "IR.fst"))
    write.fst(as.data.frame(Splice.Anno), 
        file.path(norm_output_path, "annotation", "Splice.fst"))      
    .collateData_rowEvent(irf.common, Splice.Anno, 
        norm_output_path, reference_path)
    # Write junc_PSI index
    junc_PSI = junc.common[, c("seqnames", "start", "end", "strand")]
    write.fst(junc_PSI, file.path(norm_output_path, "junc_PSI_index.fst"))
    
    message("...compiling rowEvents")
    .collateData_rowEvent(irf.common, Splice.Anno, 
        norm_output_path, reference_path)
}

################################################################################
.collateData_junc_annotate <- function(junc.common, reference_path) {
    junc.strand = read.fst(
        path = file.path(reference_path, "fst", "junctions.fst"),
        columns = c("seqnames", "start", "end", "strand", "splice_motif"),
        as.data.table = TRUE
    )
    junc.strand = unique(junc.strand)
    junc.common[, c("strand") := NULL]
    junc.common = unique(junc.common)
    junc.common = merge(junc.common, junc.strand, 
        all = TRUE, by = c("seqnames", "start", "end"))
    junc.common[is.na(get("strand")), c("strand") := "*"]
    junc.common.unanno = junc.common[get("strand") == "*"]
    junc.common.anno = junc.common[get("strand") != "*"]
    
    left.gr = with(junc.common.unanno, 
        GRanges(seqnames = seqnames, 
        ranges = IRanges(start = start, end = start + 1), 
        strand = "+"))
    right.gr = with(junc.common.unanno, 
        GRanges(seqnames = seqnames, 
        ranges = IRanges(start = end - 1, end = end), 
        strand = "+"))
    
    genome = Get_Genome(reference_path)
    junc.common.unanno[, c("splice_motif") := paste0(
        as.character(getSeq(genome, left.gr)), 
        as.character(getSeq(genome, right.gr))
    )]
    splice_motifs = data.frame(
        seqs = c("GTAG", "GCAG", "ATAC", "ATAG"),
        seqs_r = c("CTAC", "CTGC", "GTAT", "CTAT")
    )
    junc.common.unanno[ get("splice_motif") %in% splice_motifs$seqs,
        c("strand") := "+"]
    junc.common.unanno[ get("splice_motif") %in% splice_motifs$seqs_r,
        c("strand") := "-"]
    # Do not accept un-annotated GTACs - too confusing
    # Exclude unannotated non-splice motifs
    junc.common.unanno = junc.common.unanno[
        get("strand") != "*"]
    junc.common.unanno[get("strand") == "-",
        c("splice_motif") := splice_motifs$seqs[match(
            get("splice_motif"), splice_motifs$seqs_r)]]
    junc.final = rbindlist(list(junc.common.anno, junc.common.unanno))
    # Assign region names to junctions:
    junc.final[, c("Event") := paste0(get("seqnames"), ":", 
        get("start"), "-", get("end"), "/", get("strand"))]
    return(junc.final)
}

.collateData_junc_group <- function(junc.common, reference_path) {
    # Use Exon Groups file to designate exon groups to all junctions
    Exon.Groups = as.data.table(
        read.fst(file.path(reference_path, "fst", "Exons.Group.fst"))
    )
    # Always calculate stranded for junctions
    Exon.Groups.S = Exon.Groups[get("strand") != "*"]    
    
    junc.common <- .process_introns_group_overlap(
        junc.common, Exon.Groups.S,
        c("gene_group_up", "exon_group_up",
            "gene_group_down", "exon_group_down"),
        c("gene_group", "exon_group", "gene_group", "exon_group")
    )
    junc.common <- .process_introns_group_fix_RI(
        junc.common, Exon.Groups.S,
        c("gene_group_up", "exon_group_up",
            "gene_group_down", "exon_group_down"),
        c("gene_group", "intron_number", "gene_group", "intron_number")
    )
    junc.common[, 
        c("JG_up", "JG_down") := list(
            paste(get("gene_group_up"), get("exon_group_up"), sep="_"),
            paste(get("gene_group_down"), get("exon_group_down"), sep="_")
        )
    ]

    junc.common$gene_group_up = NULL
    junc.common$gene_group_down = NULL
    junc.common$exon_group_up = NULL
    junc.common$exon_group_down = NULL
    return(junc.common)
}

.collateData_irf_group <- function(
        irf.common, reference_path, stranded = TRUE) {
    # Use Exon Groups file to designate exon groups to all junctions
    Exon.Groups = as.data.table(
        read.fst(file.path(reference_path, "fst", "Exons.Group.fst"))
    )

    # Always calculate stranded for junctions
    Exon.Groups.S = Exon.Groups[get("strand") != "*"]    
    
    irf.common <- .process_introns_group_overlap(
        irf.common, Exon.Groups.S,
        c("gene_group_up", "exon_group_up",
            "gene_group_down", "exon_group_down"),
        c("gene_group", "exon_group", "gene_group", "exon_group")
    )
    irf.common <- .process_introns_group_fix_RI(
        irf.common, Exon.Groups.S,
        c("gene_group_up", "exon_group_up",
            "gene_group_down", "exon_group_down"),
        c("gene_group", "intron_number", "gene_group", "intron_number")
    )
    irf.common[, 
        c("JG_up", "JG_down") := list(
            paste(get("gene_group_up"), get("exon_group_up"), sep="_"),
            paste(get("gene_group_down"), get("exon_group_down"), sep="_")
        )
    ]

    irf.common$gene_group_up = NULL
    irf.common$gene_group_down = NULL
    irf.common$exon_group_up = NULL
    irf.common$exon_group_down = NULL

    if(!stranded) {
        Exon.Groups.S = Exon.Groups[get("strand") == "*"]
    } else {
        Exon.Groups.S = Exon.Groups[get("strand") != "*"]    
    }

    irf.common <- .process_introns_group_overlap(
        irf.common, Exon.Groups.S,
        c("gene_group_up", "exon_group_up",
            "gene_group_down", "exon_group_down"),
        c("gene_group", "exon_group", "gene_group", "exon_group")
    )
    irf.common <- .process_introns_group_fix_RI(
        irf.common, Exon.Groups.S,
        c("gene_group_up", "exon_group_up",
            "gene_group_down", "exon_group_down"),
        c("gene_group", "intron_number", "gene_group", "intron_number")
    )
    irf.common[, 
        c("IRG_up", "IRG_down") := list(
            paste(get("gene_group_up"), get("exon_group_up"), sep="_"),
            paste(get("gene_group_down"), get("exon_group_down"), sep="_")
        )
    ]

    return(irf.common)
}

.collateData_splice_anno <- function(reference_path) {
    candidate.introns = as.data.table(
        read.fst(file.path(reference_path, "fst", "junctions.fst")))
    candidate.introns[, c("transcript_biotype_2") := get("transcript_biotype")]
    candidate.introns[
        !(get("transcript_biotype") %in% 
            c("protein_coding", "processed_transcript",
            "lincRNA", "antisense", "nonsense_mediated_decay")), 
        c("transcript_biotype_2") := "other"]
    candidate.introns[, c("transcript_biotype_2") := 
        factor(get("transcript_biotype_2"), 
            c("protein_coding", "processed_transcript",
            "lincRNA", "antisense", "other", "nonsense_mediated_decay"), 
        ordered = TRUE)]
    if("transcript_support_level" %in% colnames(candidate.introns)) {
        setorderv(candidate.introns, 
            c("transcript_biotype_2", "transcript_support_level"))
    } else {
        setorder(candidate.introns, "transcript_biotype_2")    
    }
    candidate.introns[, c("Event1a") := get("Event")]
    candidate.introns[, c("Event2a") := get("Event")]

    Splice.Anno = read.fst(file.path(reference_path, "fst", "Splice.fst"),
        as.data.table = TRUE)
    Splice.Anno[candidate.introns, on = "Event1a", 
        c("up_1a") := paste(get("i.gene_group_stranded"), 
        get("i.exon_group_stranded_upstream"), sep="_")]
    Splice.Anno[candidate.introns, on = "Event1a", 
        c("down_1a") := paste(get("i.gene_group_stranded"), 
        get("i.exon_group_stranded_downstream"), sep="_")]
    Splice.Anno[candidate.introns, on = "Event2a", 
        c("down_2a") := paste(get("i.gene_group_stranded"), 
        get("i.exon_group_stranded_downstream"), sep="_")]
    
    Splice.Anno[get("EventType") %in% c("MXE", "SE", "ALE", "A3SS"),
        c("JG_up") := get("up_1a")]
    Splice.Anno[get("EventType") %in% c("SE", "AFE", "A5SS"),
        c("JG_down") := get("down_1a")]
    Splice.Anno[get("EventType") %in% c("MXE"),
        c("JG_down") := get("down_2a")]
    
    Splice.Anno$up_1a = NULL
    Splice.Anno$down_1a = NULL
    Splice.Anno$down_2a = NULL
    Splice.Anno[, c("strand") := 
        tstrsplit(get("Event1a"), split="/")[[2]]]
    return(Splice.Anno)
}

.collateData_rowEvent <- function(irf.common, Splice.Anno, 
        norm_output_path, reference_path) {
    .collateData_rowEvent_brief(irf.common, Splice.Anno, 
        norm_output_path)
    Splice.Options.Summary <- .collateData_rowEvent_splice_option(
        reference_path, Splice.Anno)
    .collateData_rowEvent_full(Splice.Options.Summary, Splice.Anno,
        norm_output_path, reference_path)
}

.collateData_rowEvent_brief <- function(irf.common, Splice.Anno, 
        norm_output_path) {
    # make rowEvent brief here
    irf.anno.brief = irf.common[, c("Name", "EventRegion")]
    setnames(irf.anno.brief, "Name", "EventName")
    irf.anno.brief[, c("EventType") := "IR"]
    irf.anno.brief = irf.anno.brief[, 
        c("EventName", "EventType", "EventRegion")]
    splice.anno.brief = Splice.Anno[, 
        c("EventName", "EventType", "EventRegion")]
    
    rowEvent = rbind(irf.anno.brief, splice.anno.brief)    
    write.fst(rowEvent, file.path(norm_output_path, "rowEvent.brief.fst"))
}
.collateData_rowEvent_splice_option <- function(reference_path,
        Splice.Anno) {
    Splice.Options = as.data.table(read.fst(
        file.path(reference_path, "fst", "Splice.options.fst")))
    Transcripts = as.data.table(read.fst(
        file.path(reference_path, "fst", "Transcripts.fst")))
    Splice.Options[Splice.Anno, on = "EventID", 
        c("EventName") := get("i.EventName")]
    Splice.Options[Transcripts, on = "transcript_id", 
        c("transcript_biotype") := get("i.transcript_biotype")]

    Splice.Options.Summary = copy(Splice.Options)
    Splice.Options.Summary[, 
        c("tsl_min") := min(get("transcript_support_level")), 
        by = c("EventID", "isoform")]
    Splice.Options.Summary[, 
        c("any_is_PC") := any(get("is_protein_coding")), 
        by = c("EventID", "isoform")]  
    Splice.Options.Summary[, 
        c("all_is_NMD") := all(grepl("decay", get("transcript_biotype"))), 
        by = c("EventID", "isoform")]
    return(Splice.Options.Summary)
}

.collateData_rowEvent_full <- function(Splice.Options.Summary, Splice.Anno,
        norm_output_path, reference_path) {
    rowEvent.Extended = read.fst(
        file.path(norm_output_path, "rowEvent.brief.fst"),
        as.data.table = TRUE)
    IR_NMD = read.fst(file.path(reference_path, "fst", "IR.NMD.fst"),
        as.data.table = TRUE)
    candidate.introns = as.data.table(
        read.fst(file.path(reference_path, "fst", "junctions.fst")))
        
    rowEvent.Extended[get("EventType") == "IR", 
        c("intron_id") := tstrsplit(get("EventName"), split="/")[[2]]]
    rowEvent.Extended[, c("Inc_Is_Protein_Coding") := FALSE]
    rowEvent.Extended[, c("Exc_Is_Protein_Coding") := FALSE]
    rowEvent.Extended[IR_NMD, on = "intron_id", 
        c("Exc_Is_Protein_Coding") := TRUE]
    rowEvent.Extended[IR_NMD, on = "intron_id", 
        c("Inc_Is_Protein_Coding") := (get("i.intron_type") == "CDS")]

    rowEvent.Extended[Splice.Options.Summary[get("isoform") == "A"], 
        on = "EventName", c("Inc_Is_Protein_Coding") := get("i.any_is_PC")]
    rowEvent.Extended[Splice.Options.Summary[get("isoform") == "B"], 
        on = "EventName", c("Exc_Is_Protein_Coding") := get("i.any_is_PC")]
    rowEvent.Extended[, c("Inc_Is_NMD", "Exc_Is_NMD") := list(FALSE, FALSE)]
    rowEvent.Extended[IR_NMD[!is.na(get("splice_is_NMD"))], on = "intron_id", 
        c("Exc_Is_NMD") := get("i.splice_is_NMD")]
    rowEvent.Extended[IR_NMD, on = "intron_id", 
        c("Inc_Is_NMD") := get("i.IRT_is_NMD")]
    rowEvent.Extended[get("EventType") == "IR" & 
        get("Exc_Is_Protein_Coding") == FALSE, c("Exc_Is_NMD") := NA]
    rowEvent.Extended[get("EventType") == "IR" & 
        get("Inc_Is_Protein_Coding") == FALSE, c("Inc_Is_NMD") := NA]

    rowEvent.Extended[Splice.Options.Summary[get("isoform") == "A"], 
        on = "EventName", c("Inc_Is_NMD") := get("i.all_is_NMD")]
    rowEvent.Extended[Splice.Options.Summary[get("isoform") == "B"], 
        on = "EventName", c("Exc_Is_NMD") := get("i.all_is_NMD")]
    rowEvent.Extended[candidate.introns, on = "intron_id", 
        c("Inc_TSL") := get("i.transcript_support_level")]
    rowEvent.Extended[candidate.introns, on = "intron_id", 
        c("Exc_TSL") := get("i.transcript_support_level")]
    rowEvent.Extended[Splice.Options.Summary[get("isoform") == "A"], 
        on = "EventName", c("Inc_TSL") := get("i.tsl_min")]
    rowEvent.Extended[Splice.Options.Summary[get("isoform") == "B"], 
        on = "EventName", c("Exc_TSL") := get("i.tsl_min")]
    # define Event1 / Event2
    rowEvent.Extended[get("EventType") == "IR", 
        c("Event1a") := get("EventRegion")]
    rowEvent.Extended[Splice.Anno, on = "EventName",
        c("Event1a", "Event2a", "Event1b", "Event2b") := 
        list(get("i.Event1a"), get("i.Event2a"), 
            get("i.Event1b"), get("i.Event2b"))]
    write.fst(rowEvent.Extended, file.path(norm_output_path, "rowEvent.fst"))
}

################################################################################
# Sub

.collateData_compile_agglist <- function(x, jobs, df.internal, 
        norm_output_path, IRMode, low_memory_mode) {
    # suppressPackageStartupMessages({
        # requireNamespace("data.table")
        # requireNamespace("stats")
    # })
    rowEvent = as.data.table(read.fst(
        file.path(norm_output_path, "rowEvent.brief.fst")))
    junc.common = as.data.table(read.fst(
        file.path(norm_output_path, "annotation", "Junc.fst")))
    irf.common = as.data.table(read.fst(
        file.path(norm_output_path, "annotation","IR.fst")))
    Splice.Anno = as.data.table(read.fst(
        file.path(norm_output_path, "annotation","Splice.fst")))
    junc_PSI = as.data.table(read.fst(
        file.path(norm_output_path, "junc_PSI_index.fst")
    ))
    work = jobs[[x]]
    block = df.internal[work]
    assays <- .collateData_assays_init(rowEvent, junc_PSI)

    for(i in seq_len(length(work))) {
        junc <- .collateData_process_junc(
            block$sample[i], block$strand[i], 
                junc.common, norm_output_path)
        splice <- .collateData_process_splice(
            junc, Splice.Anno)
        irf <- .collateData_process_irf(
            block$sample[i], block$strand[i], junc,
                irf.common, norm_output_path)                
        splice <- .collateData_process_splice_depth(
            splice, irf)
        assays <- .collateData_process_assays(assays,
            block$sample[i], junc, irf, splice, IRMode)
        file.remove(file.path(norm_output_path, "temp", 
            paste(block$sample[i], "junc.fst.tmp", sep=".")))
        file.remove(file.path(norm_output_path, "temp", 
            paste(block$sample[i], "irf.fst.tmp", sep=".")))               
    } # end FOR loop
    
    # if(low_memory_mode == TRUE) {
        # .collateData_save_assays_lowmem(assays, norm_output_path, x)
        # .collateData_save_assays_lowmem_fst(assays, norm_output_path, x)
        # return(NULL)
    # } else {

        # return(final)
    # }

    final = list(
        Included = DelayedArray(
            assays[["Included"]][, -seq_len(3), with = FALSE]),
        Excluded = DelayedArray(
            assays[["Excluded"]][, -seq_len(3), with = FALSE]),
        Depth = DelayedArray(
            assays[["Depth"]][, -seq_len(3), with = FALSE]),
        Coverage = DelayedArray(
            assays[["Coverage"]][, -seq_len(3), with = FALSE]),
        minDepth = DelayedArray(
            assays[["minDepth"]][, -seq_len(3), with = FALSE]),
        Up_Inc = DelayedArray(
            assays[["Up_Inc"]][, -seq_len(3), with = FALSE]),
        Down_Inc = DelayedArray(
            assays[["Down_Inc"]][, -seq_len(3), with = FALSE]),
        Up_Exc = DelayedArray(
            assays[["Up_Exc"]][, -seq_len(3), with = FALSE]),
        Down_Exc = DelayedArray(
            assays[["Down_Exc"]][, -seq_len(3), with = FALSE]),
        junc_PSI = DelayedArray(
            assays[["junc_PSI"]][, -seq_len(4), with = FALSE]),
        junc_counts = DelayedArray(
            assays[["junc_counts"]][, -seq_len(4), with = FALSE])
    )
    return(final)
}



.collateData_assays_init <- function(rowEvent, junc_PSI) {
    assays = list(
        Included = copy(rowEvent),
        Excluded = copy(rowEvent),
        Depth = copy(rowEvent),
        Coverage = copy(rowEvent),
        minDepth = copy(rowEvent),
        Up_Inc = rowEvent[get("EventType") %in% c("IR", "MXE", "SE")],
        Down_Inc = rowEvent[get("EventType") %in% c("IR", "MXE", "SE")],
        Up_Exc = rowEvent[get("EventType") %in% c("MXE")],
        Down_Exc = rowEvent[get("EventType") %in% c("MXE")],
        junc_PSI = copy(junc_PSI),
        junc_counts = copy(junc_PSI)
    )
    return(assays)
}

.collateData_process_junc <- function(sample, strand, 
        junc.common, norm_output_path) {
    junc = as.data.table(
        read.fst(file.path(norm_output_path, "temp", 
            paste(sample, "junc.fst.tmp", sep=".")))
    )
    junc[, c("start") := get("start") + 1]
    junc$strand = NULL

    junc = junc[junc.common, on = colnames(junc.common)[c(1,2,3)]]
    if(strand == 0) {
        junc$count = junc$total    
    } else if(strand == -1) {
        junc$count = 0
        junc[get("strand") == "+", c("count") := get("neg")]
        junc[get("strand") == "-", c("count") := get("pos")]
        junc[get("strand") == "*", c("count") := get("total")]
    } else {
        junc$count = 0
        junc[get("strand") == "+", c("count") := get("pos")]
        junc[get("strand") == "-", c("count") := get("neg")]
        junc[get("strand") == "*", c("count") := get("total")]    
    }
    junc[is.na(get("count")), c("count") := 0]
    junc = junc[,c("seqnames", "start", "end", "strand", "Event", "count")]
    junc = cbind(junc, junc.common[, c("JG_up", "JG_down")])
    junc[, c("SO_L") := 0]
    junc[, c("SO_R") := 0]
    junc[, c("SO_I") := 0]

    # SpliceLeft and SpliceRight calculations
    junc[, c("SL") := sum(get("count")), 
        by = c("seqnames", "start", "strand")]
    junc[, c("SR") := sum(get("count")), 
        by = c("seqnames", "end", "strand")]

    # first overlap any junction that has non-same-island junctions
    junc[get("JG_up") != get("JG_down") & 
            get("JG_up") != "" & get("strand") == "+", 
        c("SO_L") := sum(get("count")), by = "JG_up"]
    junc[get("JG_up") != get("JG_down") & 
            get("JG_down") != "" & get("strand") == "+", 
        c("SO_R") := sum(get("count")), by = "JG_down"]
    junc[get("JG_up") != get("JG_down") & 
            get("JG_up") != "" & get("strand") == "-",
        c("SO_R") := sum(get("count")), by = "JG_up"]
    junc[get("JG_up") != get("JG_down") & 
            get("JG_down") != "" & get("strand") == "-", 
        c("SO_L") := sum(get("count")), by = "JG_down"]

    # message("Calculating SpliceOver for annotated IR events")
            
    # Then use a simple overlap method to account for the remainder
    junc.subset = junc[get("JG_up") == get("JG_down") & 
        get("JG_up") != "" & get("JG_down") != ""]
    junc.from = .grDT(junc.subset)
    junc.to = .grDT(junc)
    
    OL = findOverlaps(junc.from, junc.to)

    splice.overlaps.DT = data.table(from = from(OL), to = to(OL))
    splice.overlaps.DT[, 
        c("count") := junc$count[to(OL)]]
    splice.overlaps.DT[, 
        c("count_sum") := sum(get("count")), by = "from"]
    splice.summa = unique(
        splice.overlaps.DT[, c("from", "count_sum")])        

    junc.subset[splice.summa$from, 
        c("SO_I") := splice.summa$count_sum]

    junc[junc.subset, on = c("Event"), c("SO_I") := get("i.SO_I")]

    # For annotated junctions, take SpliceOver as max of 
    # SpliceLeft, SpliceRight, or SpliceOver
    junc[get("SO_L") < get("SO_I"), c("SO_L") := get("SO_I")]
    junc[get("SO_R") < get("SO_I"), c("SO_R") := get("SO_I")]
    # Finally, for extreme cases, make SO_L = SL if underestimates
    junc[get("SO_L") < get("SL"), c("SO_L") := get("SL")]
    junc[get("SO_R") < get("SR"), c("SO_R") := get("SR")]
    junc[, c("SO_I") := NULL]
    
    return(junc)
}

.collateData_process_splice <- function(junc, 
        Splice.Anno) {
    splice = copy(Splice.Anno)
    
    splice[, c("count_Event1a", "count_Event2a",
        "count_Event1b", "count_Event2b") := list(0,0,0,0)]
    splice[!is.na(get("Event1a")), 
        c("count_Event1a") := junc$count[match(get("Event1a"), junc$Event)]]
    splice[is.na(get("count_Event1a")), 
        c("count_Event1a") := 0]
    splice[!is.na(get("Event2a")), 
        c("count_Event2a") := junc$count[match(get("Event2a"), junc$Event)]]
    splice[is.na(get("count_Event2a")), 
        c("count_Event2a") := 0]
    splice[!is.na(get("Event1b")), 
        c("count_Event1b") := junc$count[match(get("Event1b"), junc$Event)]]
    splice[is.na(get("count_Event1b")), 
        c("count_Event1b") := 0]
    splice[!is.na(get("Event2b")), 
        c("count_Event2b") := junc$count[match(get("Event2b"), junc$Event)]]
    splice[is.na(get("count_Event2b")), 
        c("count_Event2b") := 0]

    splice[, c("count_JG_up", "count_JG_down") := list(0,0)]
    splice[!is.na(get("JG_up")) & get("strand") == "+", 
        c("count_JG_up") := junc$SO_L[match(get("JG_up"), junc$JG_up)]]
    splice[!is.na(get("JG_up")) & get("strand") == "-", 
        c("count_JG_up") := junc$SO_R[match(get("JG_up"), junc$JG_up)]]
    splice[is.na(get("count_JG_up")), c("count_JG_up") := 0]
    splice[!is.na(get("JG_down")) & get("strand") == "-", 
        c("count_JG_down") := junc$SO_L[match(get("JG_down"), junc$JG_down)]]
    splice[!is.na(get("JG_down")) & get("strand") == "+", 
        c("count_JG_down") := junc$SO_R[match(get("JG_down"), junc$JG_down)]]
    splice[is.na(get("count_JG_down")), c("count_JG_down") := 0]

    # Splice participation: sum of two events compared to JG_up / JG_down
    splice[, c("partic_up", "partic_down") := list(0, 0)]
    splice[get("EventType") %in% c("MXE", "SE", "ALE", "A3SS"), 
        c("partic_up") := get("count_Event1a") + get("count_Event1b")]
    splice[get("EventType") %in% c("MXE"), 
        c("partic_down") := get("count_Event2a") + get("count_Event2b")]
    splice[get("EventType") %in% c("SE"), 
        c("partic_down") := get("count_Event2a") + get("count_Event1b")]
    splice[get("EventType") %in% c("AFE", "A5SS"), 
        c("partic_down") := get("count_Event1a") + get("count_Event1b")]

    # Splice coverage = participation / max_JG
    splice[, c("cov_up", "cov_down") := list(0, 0)]
    splice[get("count_JG_up") > 0, 
        c("cov_up") := get("partic_up") / get("count_JG_up")]
    splice[get("count_JG_down") > 0, 
        c("cov_down") := get("partic_down") / get("count_JG_down")]
    splice[get("EventType") %in% c("MXE", "SE") & 
        get("cov_up") < get("cov_down"), 
        c("coverage") := get("cov_up")]
    splice[get("EventType") %in% c("MXE", "SE") & 
        get("cov_up") >= get("cov_down"), 
        c("coverage") := get("cov_down")]
    splice[get("EventType") %in% c("ALE", "A3SS"), 
        c("coverage") := get("cov_up")]
    splice[get("EventType") %in% c("AFE", "A5SS"), 
        c("coverage") := get("cov_down")]        
        
    return(splice)
}

.collateData_process_irf <- function(sample, strand, junc,
        irf.common, norm_output_path) {
    irf = as.data.table(
            read.fst(file.path(norm_output_path, "temp", 
                paste(sample, "irf.fst.tmp", sep=".")))
    )

    irf[, c("start") := get("start") + 1]
    irf = irf[irf.common, on = colnames(irf.common)[seq_len(6)], 
        c("EventRegion") := get("i.EventRegion")]
    
    # Extra statistics:
    irf[, c("SpliceMax") := 0]
    irf[get("SpliceLeft") >= get("SpliceRight"), 
        c("SpliceMax") := get("SpliceLeft")]
    irf[get("SpliceLeft") < get("SpliceRight"), 
        c("SpliceMax") := get("SpliceRight")]

    irf[junc, on = c("seqnames", "start", "end", "strand"), 
        c("SpliceOverLeft") := get("SO_L")]
    irf[junc, on = c("seqnames", "start", "end", "strand"), 
        c("SpliceOverRight") := get("SO_R")]
    irf[get("SpliceOverLeft") >= get("SpliceOverRight"), 
        c("SpliceOverMax") := get("SpliceOverLeft")]
    irf[get("SpliceOverLeft") < get("SpliceOverRight"), 
        c("SpliceOverMax") := get("SpliceOverRight")]
    
    irf[, c("IROratio") := 0]
    irf[get("IntronDepth") < 1 & get("IntronDepth") > 0 & 
        (get("Coverage") + get("SpliceOverMax")) > 0, 
        c("IROratio") := get("Coverage") / (
            get("Coverage") + get("SpliceOverMax"))]
    irf[get("IntronDepth") >= 1, 
        c("IROratio") := get("IntronDepth") / 
            (get("IntronDepth") + get("SpliceOverMax"))]

    irf[, c("TotalDepth") := get("IntronDepth") + get("SpliceOverMax")]
    setnames(irf, "Name", "EventName")
    return(irf)
}

.collateData_process_splice_depth <- function(splice, irf) {
    splice.no_region = splice[!(get("EventRegion") %in% irf$EventRegion)]
    splice.no_region[, c("Depth1a") := 
        irf$TotalDepth[match(get("Event1a"), irf$EventRegion)]]
    splice.no_region[, c("Depth2a") := 
        irf$TotalDepth[match(get("Event2a"), irf$EventRegion)]]
    splice.no_region[, c("Depth1b") := 
        irf$TotalDepth[match(get("Event1b"), irf$EventRegion)]]
    splice.no_region[, c("Depth2b") := 
        irf$TotalDepth[match(get("Event2b"), irf$EventRegion)]]
    splice.no_region[, c("Depth") := 0]
    splice.no_region[get("count_JG_up") > get("count_JG_down"), 
        c("Depth") := get("count_JG_up")]
    splice.no_region[get("count_JG_up") <= get("count_JG_down"), 
        c("Depth") := get("count_JG_down")]
    splice.no_region[is.na(get("Depth1a")), c("Depth1a") := 0]
    splice.no_region[is.na(get("Depth1b")), c("Depth1b") := 0]
    splice.no_region[is.na(get("Depth2a")), c("Depth2a") := 0]
    splice.no_region[is.na(get("Depth2b")), c("Depth2b") := 0]
    splice.no_region[get("Depth1a") > get("Depth2a"), 
        c("DepthA") := get("Depth1a")]
    splice.no_region[get("Depth1b") > get("Depth2b"), 
        c("DepthB") := get("Depth1b")]
    splice.no_region[get("Depth1a") <= get("Depth2a"), 
        c("DepthA") := get("Depth2a")]
    splice.no_region[get("Depth1b") <= get("Depth2b"), 
        c("DepthB") := get("Depth2b")]
    splice.no_region[get("DepthA") > get("DepthB"), 
        c("Depth") := get("DepthA")]
    splice.no_region[get("DepthA") <= get("DepthB"), 
        c("Depth") := get("DepthB")]

    splice[, c("TotalDepth") := 0]
    splice[irf, on = "EventRegion", 
        c("TotalDepth") := get("i.TotalDepth")]
    splice[splice.no_region, on = "EventName", 
        c("TotalDepth") := get("i.Depth")]
    return(splice)
}

.collateData_process_assays <- function(assays, 
        sample, junc, irf, splice, IRMode) {
    assays[["Included"]][, c(sample) := c(
        irf$IntronDepth, 
        0.5 * (splice$count_Event1a[splice$EventType %in% c("SE", "MXE")] + 
            splice$count_Event2a[splice$EventType %in% c("SE", "MXE")]),
        splice$count_Event1a[!splice$EventType %in% c("SE", "MXE")]
    )]
    if(IRMode == "SpliceOverMax") {
        assays[["Excluded"]][, c(sample) := c(
            irf$SpliceOverMax,
            0.5 * (splice$count_Event1b[splice$EventType %in% c("MXE")] + 
                splice$count_Event2b[splice$EventType %in% c("MXE")]),
            splice$count_Event1b[!splice$EventType %in% c("MXE")]
        )]
    } else {
        assays[["Excluded"]][, c(sample) := c(
            irf$SpliceMax,
            0.5 * (splice$count_Event1b[splice$EventType %in% c("MXE")] + 
                splice$count_Event2b[splice$EventType %in% c("MXE")]),
            splice$count_Event1b[!splice$EventType %in% c("MXE")]
        )]      
    }
    # Validity checking for IR, MXE, SE
    irf[get("strand") == "+", c("Up_Inc") := get("ExonToIntronReadsLeft")]
    irf[get("strand") == "-", c("Up_Inc") := get("ExonToIntronReadsRight")]
    irf[get("strand") == "+", c("Down_Inc") := get("ExonToIntronReadsRight")]
    irf[get("strand") == "-", c("Down_Inc") := get("ExonToIntronReadsLeft")]   
    assays[["Up_Inc"]][, c(sample) := 
        c(irf$Up_Inc, 
            splice$count_Event1a[splice$EventType %in% c("MXE", "SE")])]
    assays[["Down_Inc"]][, c(sample) := 
        c(irf$Down_Inc, 
            splice$count_Event2a[splice$EventType %in% c("MXE", "SE")])]
    assays[["Up_Exc"]][, c(sample) := 
        splice$count_Event1b[splice$EventType %in% c("MXE")]]
    assays[["Down_Exc"]][, c(sample) := 
        splice$count_Event2b[splice$EventType %in% c("MXE")]]
    
    assays[["Depth"]][, c(sample) := 
        c(irf$TotalDepth, splice$TotalDepth)]
    assays[["Coverage"]][, c(sample) := 
        c(irf$Coverage, splice$coverage)]
    
    splice[get("EventType") %in% c("MXE", "SE") & 
        get("cov_up") < get("cov_down"), 
        c("minDepth") := get("count_JG_up")]
    splice[get("EventType") %in% c("MXE", "SE") & 
        get("cov_up") >= get("cov_down"), 
        c("minDepth") := get("count_JG_down")]
    splice[get("EventType") %in% c("ALE", "A3SS"), 
        c("minDepth") := get("count_JG_up")]
    splice[get("EventType") %in% c("AFE", "A5SS"), 
        c("minDepth") := get("count_JG_down")]   
    assays[["minDepth"]][, c(sample) := c(
        irf$IntronDepth,
        splice$minDepth)]

    junc[get("count") == 0, c("PSI") := 0]
    junc[get("SO_L") > get("SO_R"), 
        c("PSI") := get("count") / get("SO_L")]
    junc[get("SO_R") >= get("SO_L") & get("SO_R") > 0, 
        c("PSI") := get("count") / get("SO_R")]               
    assays[["junc_PSI"]][junc, on = c("seqnames", "start", "end", "strand"),
        c(sample) := get("i.PSI")]
    assays[["junc_counts"]][junc, on = c("seqnames", "start", "end", "strand"),
        c(sample) := get("i.count")]
        
    return(assays)
}

.collateData_save_assays_lowmem <- function(assays, norm_output_path, x) {
    value = t(as.matrix(assays[["Included"]][, -c(1,2,3)]))
    fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
        paste("Included", as.character(x), "txt.gz", sep=".")), 
        col.names = FALSE, row.names = FALSE)
    value = t(as.matrix(assays[["Excluded"]][, -c(1,2,3)]))
    fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
        paste("Excluded", as.character(x), "txt.gz", sep=".")), 
        col.names = FALSE, row.names = FALSE)
    value = t(as.matrix(assays[["Depth"]][, -c(1,2,3)]))
    fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
        paste("Depth", as.character(x), "txt.gz", sep=".")), 
        col.names = FALSE, row.names = FALSE)
    value = t(as.matrix(assays[["Coverage"]][, -c(1,2,3)]))
    fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
        paste("Coverage", as.character(x), "txt.gz", sep=".")), 
        col.names = FALSE, row.names = FALSE)
    value = t(as.matrix(assays[["minDepth"]][, -c(1,2,3)]))
    fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
        paste("minDepth", as.character(x), "txt.gz", sep=".")), 
        col.names = FALSE, row.names = FALSE)
    value = t(as.matrix(assays[["Up_Inc"]][, -c(1,2,3)]))
    fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
        paste("Up_Inc", as.character(x), "txt.gz", sep=".")), 
        col.names = FALSE, row.names = FALSE)
    value = t(as.matrix(assays[["Down_Inc"]][, -c(1,2,3)]))
    fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
        paste("Down_Inc", as.character(x), "txt.gz", sep=".")), 
        col.names = FALSE, row.names = FALSE)
    value = t(as.matrix(assays[["Up_Exc"]][, -c(1,2,3)]))
    fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
        paste("Up_Exc", as.character(x), "txt.gz", sep=".")), 
        col.names = FALSE, row.names = FALSE)
    value = t(as.matrix(assays[["Down_Exc"]][, -c(1,2,3)]))
    fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
        paste("Down_Exc", as.character(x), "txt.gz", sep=".")), 
        col.names = FALSE, row.names = FALSE)
    value = t(as.matrix(assays[["junc_PSI"]][, -c(1,2,3,4)]))
    fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
        paste("junc_PSI", as.character(x), "txt.gz", sep=".")), 
        col.names = FALSE, row.names = FALSE)
    value = t(as.matrix(assays[["junc_counts"]][, -c(1,2,3,4)]))
    fwrite(as.data.frame(value), file.path(norm_output_path, "temp", 
        paste("junc_counts", as.character(x), "txt.gz", sep=".")), 
        col.names = FALSE, row.names = FALSE)
}

.collateData_save_assays_lowmem_fst <- function(assays, norm_output_path, x) {
    assaynames = c("Included", "Excluded", "Depth", "Coverage", "minDepth",
        "Up_Inc", "Down_Inc", "Up_Exc", "Down_Exc", 
        "junc_PSI", "junc_counts")
    for(assayname in assaynames) {
        if(grepl("junc", assayname)) {
            .collateData_save_assays_lowmem_fst_indiv(
                assays[[assayname]][, -seq_len(4), with = FALSE],
                assayname, norm_output_path, x)
        } else {
            .collateData_save_assays_lowmem_fst_indiv(
                assays[[assayname]][, -seq_len(3), with = FALSE],
                assayname, norm_output_path, x)
        }
    }
}

.collateData_save_assays_lowmem_fst_indiv <- function(mat, assayname, 
        norm_output_path, x) {
    mat = as.data.frame(mat)
    colnames(mat) = as.character(seq_len(ncol(mat)))
    fst::write.fst(mat, file.path(norm_output_path, "temp", 
        paste(assayname, as.character(x), "fst", sep=".")))
}

.collateData_compile_assays <- function(agg.list, df.internal,
        norm_output_path, low_memory_mode,
        item.todo, jobs, n_jobs) {
    # if(low_memory_mode) {
        # for(item in item.todo) {
            # file.DT = data.table(file = list.files(pattern = item, 
                # path = file.path(norm_output_path, "temp")))
            # file.DT[, c("index") := 
                # as.numeric(tstrsplit(file.DT, split=".", fixed = TRUE)[[2]])]
            # setorder(file.DT, "index")
            # mat = NULL
            # for(x in seq_len(n_jobs)) {
                # temp = t(fread(file.path(file.path(norm_output_path, "temp"), 
                    # file.DT$file[x]), data.table = FALSE))
                # temp <- fst::read.fst(file.path(norm_output_path, "temp", 
                    # file.DT$file[x]))
                # colnames(temp) = df.internal$sample[jobs[[x]]]
                # if(!is.null(mat)) {
                    # mat <- cbind(mat, temp)
                # } else mat <- temp
                # file.remove(file.path(file.path(norm_output_path, "temp"), 
                    # file.DT$file[x]))
            # }
            # rownames(mat) <- seq_len(nrow(mat))
            # outfile = file.path(norm_output_path, paste(item, "fst", sep="."))
            # write.fst(as.data.frame(mat), outfile)
        # }
    # } else {
        # item.DTList = list()
        # for(item in item.todo) {
            # for(x in seq_len(n_jobs)) {
                # if(x == 1) {
                    # item.DTList[[item]] = agg.list[[x]][[item]]
                # } else {
                    # item.DTList[[item]] = cbind(item.DTList[[item]], 
                        # agg.list[[x]][[item]])
                # }
            # }
            # outfile = file.path(norm_output_path, paste(item, "fst", sep="."))
            # write.fst(as.data.frame(item.DTList[[item]]), outfile)
        # }
    # }

    junc_index <- fst::read.fst(file.path(
        norm_output_path, "junc_PSI_index.fst"
    ))
    junc_rownames <- with(junc_index, 
        paste0(seqnames, ":", start, "-", end, "/", strand))
    # Realize all the DelayedArrays as H5:
    assays = list()
    agg_t = purrr::transpose(agg.list)
    outfile = file.path(norm_output_path, paste("data", "h5", sep="."))
    if(file.exists(outfile)) file.remove(outfile)
    item.DTList = list()
    for(item in item.todo) {
        # for(x in seq_len(n_jobs)) {
            # if(x == 1) {
                # item.DTList[[item]] = agg.list[[x]][[item]]
            # } else {
                # item.DTList[[item]] = cbind(item.DTList[[item]], 
                    # agg.list[[x]][[item]])
            # }
        # }
        item.DTList[[item]] = do.call(cbind, agg_t[[item]])
        
        if(grepl("junc", item)) {
            rownames(item.DTList[[item]]) <- junc_rownames
        }
        
        assays[[item]] <- HDF5Array::writeHDF5Array(item.DTList[[item]], 
            filepath = outfile,
            name = item, with.dimnames = TRUE
        )
    }
    return(assays)
}

.collateData_fix_juncnames <- function(norm_output_path) {
    # Rewrite junc_PSI and junc_count by adding in NxtIRF rownames 
    #   (as first column)
    junc_index = fst::read.fst(file.path(
        norm_output_path, "junc_PSI_index.fst"
    ))
    junc_PSI = fst::read.fst(file.path(
        norm_output_path, "junc_PSI.fst"
    ))  
    junc_counts = fst::read.fst(file.path(
        norm_output_path, "junc_counts.fst"
    ))
    junc_PSI$rownames = with(junc_index, 
        paste0(seqnames, ":", start, "-", end, "/", strand))
    junc_counts$rownames = junc_PSI$rownames
    fst::write.fst(cbind(junc_PSI[,ncol(junc_PSI),drop=FALSE],
        junc_PSI[,-ncol(junc_PSI)]),file.path(
        norm_output_path, "junc_PSI.fst"
    ))
    fst::write.fst(cbind(junc_counts[,ncol(junc_counts),drop=FALSE],
        junc_counts[,-ncol(junc_counts)]),file.path(
        norm_output_path, "junc_counts.fst"
    ))
}

.collateData_write_stats <- function(df.internal, norm_output_path) {
    outfile = file.path(norm_output_path, paste("stats", "fst", sep="."))
    write.fst(as.data.frame(df.internal), outfile)
}

.collateData_write_colData <- function(df.internal, coverage_files,
        norm_output_path) {
    # Create barebones colData.Rds - save coverage files as well
    if(length(coverage_files) == nrow(df.internal) & IsCOV(coverage_files)) {
        df.files = data.table(
            sample = df.internal$sample,
            bam_file = "",
            irf_file = df.internal$path,
            cov_file = coverage_files
        )
    } else {
        df.files = data.table(
            sample = df.internal$sample,
            bam_file = "",
            irf_file = df.internal$path,
            cov_file = ""
        )
    }
    # 
    colData = list(
        df.files = df.files,
        df.anno = data.table(sample = df.internal$sample)
    )
    saveRDS(colData, file.path(norm_output_path, "colData.Rds"))
}


################################################################################

#' Constructs a SummarizedExperiment object from the collated data
#'
#' MakeSE() creates a SummarizedExperiment object from the data collated
#' from IRFinder output using CollateData().
#'
#' @param collate_path The output path given to CollateData() pointing to the
#'   collated data
#' @param colData A data frame containing the sample annotation information.
#'   Note that the first column must contain the sample names. If the names of 
#'   only a subset of samples are given, then `MakeSE()` will construct the SE 
#'   object based only on the samples given. Omit `colData` to generate an SE 
#'   object based on the whole dataset. The colData can be set later using 
#'   `colData()`
#' @param RemoveOverlapping (default = TRUE) Whether to filter out overlapping 
#'   introns of IR events belonging to minor isoforms. MakeSE will try to 
#'   identify which junctions belong to major isoforms, then select the 
#'   junctions from non-overlapping minor isoforms in an iterative approach, 
#'   until no non-overlapping introns remain. This is important
#'   to make sure IR events are not 'double-counted'
#'
#' @return A NxtIRF SummarizedExperiment (`NxtSE`) object 
#'
#' @examples
#' se = MakeSE(collate_path = file.path(tempdir(), "NxtIRF_output"))
#' @md
#' @export
MakeSE = function(collate_path, colData, RemoveOverlapping = TRUE) {
    # Includes iterative filtering for IR events with highest mean PSI
        # To annotate IR events of major isoforms

    colData <- .makeSE_validate_args(collate_path, colData)
    colData <- .makeSE_colData_clean(colData)

    N <- 8
    dash_progress("Loading NxtSE object from file...", N)
    message("Loading NxtSE object from file...", appendLF = FALSE)
    # se <- .makeSE_initialise_HDF5(collate_path, colData)
    se <- HDF5Array::loadHDF5SummarizedExperiment(
        dir = collate_path, prefix = "NxtSE_")
    # Encapsulate as NxtSE object
    se = se[, colData$sample]
    if(ncol(colData) > 1) {
        colData_use <- colData[, -1, drop = FALSE]
        rownames(colData_use) <- colData$sample
        colData(se) <- as(colData_use, "DataFrame")    
    }
    se = as(se, "NxtSE")
    message("done\n")
    
    if(RemoveOverlapping == TRUE) {
        # Iterative filtering of IR
        # tryCatch({
            # se <- .makeSE_iterate_IR(se, collate_path)
        # }, error = function(e) {
            # message(paste(
                # "Iterative filtering of IR appears to have",'
                # "run into an error.",
                # "Using RemoveOverlapping = FALSE"))
        # })
        dash_progress("Removing overlapping introns...", N)
        se <- .makeSE_iterate_IR(se, collate_path)
    }
    return(se)
}

################################################################################
# helpers

.makeSE_validate_args <- function(collate_path, colData) {
    item.todo = c("rowEvent", "Included", "Excluded", "Depth", "Coverage", 
        "minDepth", "Up_Inc", "Down_Inc", "Up_Exc", "Down_Exc")
    # files.todo = file.path(normalizePath(collate_path), 
        # paste(item.todo, "fst", sep="."))
    # if(!all(file.exists(files.todo))) {
        # stop(paste(
            # "FST File generation appears incomplete.",
            # "Suggest run CollateData() again"
        # ), call. = FALSE)
    # }
    if(!file.exists(file.path(collate_path, "colData.Rds"))) {
        stop(paste("In MakeSE():",
            file.path(collate_path, "colData.Rds"),
            "was not found"
        ), call. = FALSE)
    }
    colData.Rds = readRDS(file.path(collate_path, "colData.Rds"))
    if(!("df.anno" %in% names(colData.Rds))) {
        stop(paste("In MakeSE():",
            file.path(collate_path, "colData.Rds"),
            "must contain df.anno containing annotations"
        ), call. = FALSE)
    }
    if(missing(colData)) {    
        colData = colData.Rds$df.anno
    } else {
        if(!("sample" %in% colnames(colData))) {
            stop(paste("In MakeSE():",
                "'sample' must be the name of the first column",
                "in colData, containing sample names"
            ), call. = FALSE)
        }
        if(!all(colData$sample %in% colData.Rds$df.anno$sample)) {
            stop(paste("In MakeSE():",
                "some samples in colData were not found in given path"
            ), call. = FALSE)
        }
    }
    colData = as.data.frame(colData)
    return(colData)
}

.makeSE_colData_clean <- function(colData) {
    remove_na = NULL
    if(ncol(colData) > 1) {
        for(i in seq(2, ncol(colData))) {
            if(is(colData[,i], "character")) {
                colData[,i] = factor(unlist(colData[,i]))      
            } else if(is(colData[,i], "logical")) {
                colData[,i] <- factor(unlist(
                    ifelse(colData[,i], "TRUE","FALSE")))                
            } else if(all(is.na(unlist(colData[,i])))) {
                remove_na = append(remove_na, i)
            }
        }
    }
    if(!is.null(remove_na)) {
        colData = colData[,-remove_na]
    }
    return(colData)
}

.makeSE_initialise <- function(collate_path, colData) {
    item.todo = c("rowEvent", "Included", "Excluded", "Depth", "Coverage", 
        "minDepth", "Up_Inc", "Down_Inc", "Up_Exc", "Down_Exc")
    files.todo = file.path(normalizePath(collate_path), 
        paste(item.todo, "fst", sep="."))
    colData.Rds = readRDS(file.path(collate_path, "colData.Rds"))
    rowData = read.fst(files.todo[1])
    Included = as.matrix(read.fst(files.todo[2], columns = colData$sample))
    Excluded = as.matrix(read.fst(files.todo[3], columns = colData$sample))
    Depth = as.matrix(read.fst(files.todo[4], columns = colData$sample))
    Coverage = as.matrix(read.fst(files.todo[5], columns = colData$sample))
    minDepth = as.matrix(read.fst(files.todo[6], columns = colData$sample))
    Up_Inc = as.matrix(read.fst(files.todo[7], columns = colData$sample))
    Down_Inc = as.matrix(read.fst(files.todo[8], columns = colData$sample))
    Up_Exc = as.matrix(read.fst(files.todo[9], columns = colData$sample))
    Down_Exc = as.matrix(read.fst(files.todo[10], columns = colData$sample))
    rownames(Up_Inc) = rowData$EventName[
        rowData$EventType %in% c("IR", "MXE", "SE")]
    rownames(Down_Inc) = rowData$EventName[
        rowData$EventType %in% c("IR", "MXE", "SE")]
    rownames(Up_Exc) = rowData$EventName[
        rowData$EventType %in% c("MXE")]
    rownames(Down_Exc) = rowData$EventName[
        rowData$EventType %in% c("MXE")]
    
    # Annotate NMD direction
    rowData = as.data.table(rowData)
    rowData[, c("NMD_direction") := 0]
    rowData[get("Inc_Is_NMD") & !get("Exc_Is_NMD"), c("NMD_direction") := 1]
    rowData[!get("Inc_Is_NMD") & get("Exc_Is_NMD"), c("NMD_direction") := -1]
    rowData = as.data.frame(rowData)

    se = SummarizedExperiment(
        assays = SimpleList(
            Included = Included, Excluded = Excluded, 
            Depth = Depth, Coverage = Coverage, minDepth = minDepth
        ),
        rowData = rowData, 
        colData = as.data.frame(
            colData[, -1, drop=FALSE], row.names = colData$sample)
    )
    rownames(se) = rowData(se)$EventName

    metadata(se)$Up_Inc = Up_Inc
    metadata(se)$Down_Inc = Down_Inc
    metadata(se)$Up_Exc = Up_Exc
    metadata(se)$Down_Exc = Down_Exc
    if("df.files" %in% names(colData.Rds) &&
        "cov_file" %in% colnames(colData.Rds$df.files)) {
        metadata(se)$cov_file = colData.Rds$df.files$cov_file
        # names(metadata(se)$cov_file) = colData.Rds$df.files$sample       
    }
    metadata(se)$ref = readRDS(file.path(
        collate_path, "Annotation", "cov_data.Rds"
    ))
    return(se)
}

.makeSE_initialise_HDF5 <- function(collate_path, colData) {
    item.todo = c("rowEvent", "Included", "Excluded", "Depth", "Coverage", 
        "minDepth", "Up_Inc", "Down_Inc", "Up_Exc", "Down_Exc")
    files.todo = file.path(normalizePath(collate_path), 
        paste(item.todo, "fst", sep="."))
    colData.Rds = readRDS(file.path(collate_path, "colData.Rds"))
    h5file = file.path(collate_path, "data.h5")
    rowData = read.fst(file.path(collate_path, "rowEvent.fst"))
    Included = HDF5Array(h5file, item.todo[2])[,colData$sample]
    Excluded = HDF5Array(h5file, item.todo[3])[,colData$sample]
    Depth = HDF5Array(h5file, item.todo[4])[,colData$sample]
    Coverage = HDF5Array(h5file, item.todo[5])[,colData$sample]
    minDepth = HDF5Array(h5file, item.todo[6])[,colData$sample]
    Up_Inc = HDF5Array(h5file, item.todo[7])[,colData$sample]
    Down_Inc = HDF5Array(h5file, item.todo[8])[,colData$sample]
    Up_Exc = HDF5Array(h5file, item.todo[9])[,colData$sample]
    Down_Exc = HDF5Array(h5file, item.todo[10])[,colData$sample]
    rownames(Up_Inc) = rowData$EventName[
        rowData$EventType %in% c("IR", "MXE", "SE")]
    rownames(Down_Inc) = rowData$EventName[
        rowData$EventType %in% c("IR", "MXE", "SE")]
    rownames(Up_Exc) = rowData$EventName[
        rowData$EventType %in% c("MXE")]
    rownames(Down_Exc) = rowData$EventName[
        rowData$EventType %in% c("MXE")]
    
    # Annotate NMD direction
    rowData = as.data.table(rowData)
    rowData[, c("NMD_direction") := 0]
    rowData[get("Inc_Is_NMD") & !get("Exc_Is_NMD"), c("NMD_direction") := 1]
    rowData[!get("Inc_Is_NMD") & get("Exc_Is_NMD"), c("NMD_direction") := -1]
    rowData = as.data.frame(rowData)
    colData = as.data.frame(colData)
    
    colData_use = as.data.frame(
            colData[, -1, drop=FALSE], row.names = colData$sample)
    se = SummarizedExperiment(
        assays = SimpleList(
            Included = Included, Excluded = Excluded, 
            Depth = Depth, Coverage = Coverage, minDepth = minDepth
        ),
        rowData = rowData, 
        colData = colData_use
    )
    rownames(se) = rowData(se)$EventName

    metadata(se)$Up_Inc = Up_Inc
    metadata(se)$Down_Inc = Down_Inc
    metadata(se)$Up_Exc = Up_Exc
    metadata(se)$Down_Exc = Down_Exc
    if("df.files" %in% names(colData.Rds) &&
        "cov_file" %in% colnames(colData.Rds$df.files)) {
        metadata(se)$cov_file = colData.Rds$df.files$cov_file
        # names(metadata(se)$cov_file) = colData.Rds$df.files$sample       
    }
    metadata(se)$ref = readRDS(file.path(
        collate_path, "Annotation", "cov_data.Rds"
    ))
    return(se)
}

.makeSE_iterate_IR <- function(se, collate_path) {
    # junc_PSI = fst::read.fst(file.path(
        # normalizePath(collate_path), "junc_PSI.fst"
    # ))
    # rownames(junc_PSI) = junc_PSI$rownames
    # junc_PSI = junc_PSI[,-1,drop=FALSE]
    junc_PSI <- HDF5Array(file.path(normalizePath(collate_path), 
        "data.h5"), "junc_PSI")[, colnames(se)]

    se.IR = se[rowData(se)$EventType == "IR",,drop = FALSE]
    se.coords = rowData(se.IR)$EventRegion[
        rowData(se.IR)$EventRegion %in% rownames(junc_PSI)]
    
    if(length(se.coords) > 0) {
        message(paste(
            "Iterating through IR events to determine introns",
            "of main isoforms"))
        include <- .makeSE_iterate_IR_select_events(se.coords, junc_PSI)
        se.coords.final = se.coords[include]
        se.coords.excluded = se.coords[!include]

        # Iteration to find events not overlapping with se.IR.final
        include <- .makeSE_iterate_IR_retrieve_excluded_introns(
            se.coords.final, se.coords.excluded)
        iteration = 0
        while(length(include) > 0 & length(se.coords.final) > 0) {
            iteration = iteration + 1
            message(paste("Iteration", iteration))
            dash_progress(paste("Iteration", iteration), 8)
            se.coords.excluded = se.coords.excluded[include]

            include <- .makeSE_iterate_IR_select_events(
                    se.coords.excluded, junc_PSI)

            if(length(include) > 0) {
                se.coords.final = c(se.coords.final, 
                    se.coords.excluded[include])
                se.coords.excluded = 
                    se.coords.excluded[!include]
                include <- .makeSE_iterate_IR_retrieve_excluded_introns(
                    se.coords.final, se.coords.excluded)
            } else {
                # final.gr = c()
                include = c()
            }
        }
        # se = rbind(
            # se.IR[rownames(se.IR) %in% rownames(se.IR.final),,drop = FALSE],
            # se[rowData(se)$EventType != "IR",,drop = FALSE]
        # )
        se = se[c(
            which(rowData(se.IR)$EventRegion %in% se.coords.final), 
            which(rowData(se)$EventType != "IR")
        ),]
    }
    return(se)
}

.makeSE_iterate_IR_select_events <- function(se.coords, junc_PSI) {
    gr = NxtIRF.CoordToGR(se.coords)
    gr.reduced = reduce(gr)

    OL = findOverlaps(gr, gr.reduced)
    junc_PSI.group = as.data.table(junc_PSI[se.coords,, drop = FALSE])
    junc_PSI.group$means = rowMeans(junc_PSI.group)
    junc_PSI.group$group = to(OL)
    junc_PSI.group[, c("max_means") := max(get("means")), 
        by = "group"]
    return(junc_PSI.group$means == junc_PSI.group$max_means)
}

.makeSE_iterate_IR_retrieve_excluded_introns <- function(
        se.coords.final, se.coords.excluded) {
    if(length(se.coords.excluded) > 0) {
        final.gr = NxtIRF.CoordToGR(se.coords.final)
        excluded.gr = NxtIRF.CoordToGR(se.coords.excluded)

        OL = findOverlaps(excluded.gr, final.gr)
        include = which(!(
            seq_len(length(excluded.gr))) %in% sort(unique(from(OL))))
    } else {
        include = c()
    }
    return(include)
}

################################################################################



#' Filtering for IR and Alternative Splicing Events
#'
#' This function implements filtering of IR or AS events based on customisable
#' criteria
#' 
#' @details
#'   \strong{Annotation Filters}\cr\cr 
#'     \strong{Protein_Coding}: Filters for alternative splicing or IR events 
#'       within protein reading frames. No additional parameters required.\cr\cr
#'     \strong{NMD_Switching}: Filters for events in which one isoform is a
#'       predicted NMD substrate.\cr\cr 
#'     \strong{Transcript_Support_Level}: filters for events in which both
#'       isoforms have a TSL level below or equal to filterVars$minimum\cr\cr 
#'   \strong{Data Filters}\cr\cr 
#'     \strong{Depth}: Filters IR or alternative splicing events of transcripts
#'       that are "expressed" with adequate \code{Depth} as calculated by the
#'       sum of all splicing and IR reads spanning the event. Events with 
#'       \code{Depth} below filterVars$minimum are excluded\cr\cr
#'     \strong{Coverage}: Coverage means different things to IR and alternative
#'       splicing.\cr\cr 
#'     For \emph{IR}, Coverage refers to the percentage of the measured intron
#'       covered with reads. Introns of samples with an IntronDepth above 
#'       \code{filterVars$minDepth} are assessed, with introns with coverage 
#'       below \code{filterVars$minimum} are excluded.\cr\cr 
#'     For \emph{Alternative Splicing}, Coverage refers to the percentage of all
#'       splicing events observed across the genomic region that is compatible
#'       with either the included or excluded event. This prevents NxtIRF from 
#'       doing differential analysis between two minor isoforms. Instead of
#'       IntronDepth, in AS events NxtIRF considers events where the spliced
#'       reads from both exonic regions exceed \code{filterVars$minDepth}.
#'       Then, events with a splicing coverage below \code{filterVars$minimum}
#'       are excluded. We recommend testing IR events for > 90% coverage and AS
#'       events for > 60% coverage as given in the default filters which can be
#'       accessed using \code{\link{get_default_filters}}\cr\cr  
#'     \strong{Consistency}: Skipped exons (SE) and mutually exclusive exons
#'       (MXE) comprise reads of two contiguous splice junctions (for the
#'       included casette exon). Summating counts from both junctions is
#'       misleading as there may be overlapping events (e.g. alternate first 
#'       / last exons) that only rely on one splice event. To ensure the SE /
#'       MXE is the dominant event, we require both splice junctions to have
#'       comparable counts.\cr\cr
#'     Events are excluded if either of the upstream or downstream
#'       event is lower than total splicing events by a log-2 magnitude 
#'       above filterVars$maximum. For example, if 
#'       \code{filterVars$maximum = 2}, we require both upstream and downstream
#'       events to represent at least 1/(2^2) = 1/4 of the sum of upstream
#'       and downstream event.
#'       This is considered for each isoform of each event, as long as the
#'       total counts belonging to the considered isoform is above 
#'       \code{filterVars$minDepth}.
#'
#'   We highly recommend using the default filters, which can be acquired 
#'     using \code{\link{get_default_filters}}
#' 
#' @param filterClass One of either \code{Annotation} or \code{Data}
#' @param filterType For filterClass \code{Annotation}, either one of 
#'   \code{Protein_Coding}, \code{NMD_Switching}, 
#'   \code{Transcript_Support_Level}.
#'   For filterClass \code{Data}, either one of
#'   \code{Depth}, \code{Coverage}, \code{Consistency}.
#' @param filterVars A list of parameters, as explained below
#' @param filterObject the SummarizedExperiment to filter
#'
#' @return A vector of type \code{logical} designating which events to retain 
#'   \code{TRUE} and which to remove \code{FALSE}.
#' @examples
#' # see ?MakeSE on example code of generating this NxtSE object
#' se = NxtIRF_example_NxtSE()
#' 
#' # Get NxtIRF recommended filters
#' filters = get_default_filters()
#' 
#' # Filter the SummarizedExperiment using the first default filter ("Depth")
#' se.depthfilter = se[runFilter(
#'         filterClass = filters[[1]]$filterClass,
#'         filterType = filters[[1]]$filterType,
#'         filterVars = filters[[1]]$filterVars,
#'         filterObject = se
#'     ), ]
#' @seealso \code{\link{get_default_filters}}, 
#' \code{\link{apply_filters}}
#' @export
runFilter <- function(filterClass, filterType, filterVars, filterObject) {
    # filterClass: can be one of 'Annotation', 'Data', 'Runtime'
    # filterType:
    # - Annotation:
    # - Data:
    # -     Depth: 1-minimum, 2-minCond, 3-pcTRUE
    # -     Coverage: 1-minimum, 1a-minDepth, 2-minCond, 3-pcTRUE
    # -     Consistency: 1-maximum, 1a-minDepth, 2-minCond, 3-pcTRUE
    # - for Consistency, maximum is the max(abs(log2_delta)) between comparison
    #       and calculated value
    if(filterClass == "Data") {
        if(filterType == "Depth") {
            message("Running Depth filter")
            return(.runFilter_data_depth(filterObject, filterVars))
        } else if(filterType == "Coverage") {
            message("Running Coverage filter")
            return(.runFilter_data_coverage(filterObject, filterVars))
        } else if(filterType == "Consistency") {    # requires: 
            message("Running Consistency filter")
            return(.runFilter_data_consistency(filterObject, filterVars))
        }
    } else if(filterClass == "Annotation") {
        if(filterType == "Protein_Coding") {
            # returns if any of included or excluded is protein_coding
            rowSelected = as.data.table(rowData)
            rowSelected = rowSelected[
                get("Inc_Is_Protein_Coding") == TRUE | 
                get("Exc_Is_Protein_Coding") == TRUE]
            rowSelected = rowSelected[get("EventType") != "IR" | 
                get("Inc_Is_Protein_Coding") == TRUE] # filter for CDS introns
            res = rowData$EventName %in% rowSelected$EventName
        } else if(filterType == "NMD_Switching") {
            rowSelected = as.data.table(rowData)
            rowSelected = rowSelected[!is.na(get("Inc_Is_NMD")) & 
                !is.na(get("Exc_Is_NMD"))]
            rowSelected = rowSelected[get("Inc_Is_NMD") != get("Exc_Is_NMD")]
            res = rowData$EventName %in% rowSelected$EventName
        } else if(filterType == "Transcript_Support_Level") {
            if(!("minimum" %in% names(filterVars))) {
                minimum = 1
            } else {
                minimum = filterVars$minimum
            }
            rowSelected = as.data.table(rowData)
            rowSelected = rowSelected[get("Inc_TSL") != "NA" & 
                get("Exc_TSL") != "NA"]
            rowSelected[, c("Inc_TSL") := as.numeric(get("Inc_TSL"))]
            rowSelected[, c("Exc_TSL") := as.numeric(get("Exc_TSL"))]
            rowSelected = rowSelected[get("Inc_TSL") <= minimum & 
                get("Exc_TSL") <= minimum]
            res = rowData$EventName %in% rowSelected$EventName
        }
        if("EventTypes" %in% names(filterVars)) {
            res[!(rowData$EventType %in% filterVars$EventTypes)] = TRUE
        }
        return(res)
    } else {
        return(rep(TRUE, nrow(filterObject)))
    }
}

################################################################################
# Individual functions:

.runFilter_data_depth <- function(se, filterVars) {
    usePC = ifelse("pcTRUE" %in% names(filterVars),
        filterVars$pcTRUE, 100)
    use_cond = ifelse(
            !is.null(names(filterVars)) && 
            "condition" %in% names(filterVars) && 
            filterVars$condition %in% colnames(colData), 
        TRUE, FALSE)

    colData = as.data.frame(colData(se))
    rowData = as.data.frame(rowData(se))

    minimum = ifelse("minimum" %in% names(filterVars),
        filterVars$minimum, 20)
        
    if(use_cond == TRUE) {
        cond_vec = unlist(
            colData[, which(colnames(colData) == filterVars$condition)])
        cond_vars = unique(cond_vec)
    }
    depth = as.matrix(assay(se, "Depth"))
    sum_res = rep(0, nrow(se))
    if(use_cond == TRUE) {
        for(cond in cond_vars) {
            depth.subset = depth[, which(cond_vec == cond)]
            sum = rowSums(depth.subset >= minimum)
            sum_res = sum_res + 
                ifelse(sum * 100 / ncol(depth.subset) >= usePC, 1, 0)
        }
        n_TRUE = ifelse(
            !is.na(suppressWarnings(as.numeric(filterVars$minCond))),
            as.numeric(filterVars$minCond), -1
        )
        if(n_TRUE == -1) n_TRUE = length(cond_vars)
        res = (sum_res >= n_TRUE)
    } else {
        sum = rowSums(depth >= minimum)
        res = ifelse(sum * 100 / ncol(depth) >= usePC, TRUE, FALSE)
    }
    if("EventTypes" %in% names(filterVars)) {
        res[!(rowData$EventType %in% filterVars$EventTypes)] = TRUE
    }
    return(res)
}

.runFilter_data_coverage <- function(se, filterVars) {
    usePC = ifelse("pcTRUE" %in% names(filterVars),
        filterVars$pcTRUE, 100)
    use_cond = ifelse(
            !is.null(names(filterVars)) && 
            "condition" %in% names(filterVars) && 
            filterVars$condition %in% colnames(colData), 
        TRUE, FALSE)
    minDepth = ifelse("minDepth" %in% names(filterVars),
        filterVars$minDepth, 0)

    colData = as.data.frame(colData(se))
    rowData = as.data.frame(rowData(se))

    minimum = ifelse("minimum" %in% names(filterVars),
        filterVars$minimum, 20)
    if(use_cond == TRUE) {
        cond_vec = unlist(
            colData[, which(colnames(colData) == filterVars$condition)])
        cond_vars = unique(cond_vec)
    }
    cov = as.matrix(assay(se, "Coverage"))
    depth = as.matrix(assay(se, "minDepth"))
    
    # do not test if depth below threshold
    cov[depth < minDepth] = 1    

    sum_res = rep(0, nrow(se))
    if(use_cond == TRUE) {
        for(cond in cond_vars) {
            cov.subset = cov[, which(cond_vec == cond)]
            sum = rowSums(cov.subset >= minimum / 100)
            sum_res = sum_res + 
                ifelse(sum * 100 / ncol(cov.subset) >= usePC, 1, 0)
        }
        n_TRUE = ifelse(
            !is.na(suppressWarnings(as.numeric(filterVars$minCond))), 
            as.numeric(filterVars$minCond), -1)
        if(n_TRUE == -1) n_TRUE = length(cond_vars)
        res = (sum_res >= n_TRUE)
    } else {
        sum = rowSums(cov >= minimum / 100)
        res = ifelse(sum * 100 / ncol(cov) >= usePC, TRUE, FALSE)
    }
    if("EventTypes" %in% names(filterVars)) {
        res[!(rowData$EventType %in% filterVars$EventTypes)] = TRUE
    }
    return(res)
}

.runFilter_data_consistency <- function(se, filterVars) {
    usePC = ifelse("pcTRUE" %in% names(filterVars),
        filterVars$pcTRUE, 100)
    use_cond = ifelse(
            !is.null(names(filterVars)) && 
            "condition" %in% names(filterVars) && 
            filterVars$condition %in% colnames(colData), 
        TRUE, FALSE)
    minDepth = ifelse("minDepth" %in% names(filterVars),
        filterVars$minDepth, 0)

    colData = as.data.frame(colData(se))
    rowData = as.data.frame(rowData(se))

    maximum = ifelse("maximum" %in% names(filterVars),
        filterVars$maximum, 1)

    if(use_cond == TRUE) {
        cond_vec = unlist(colData[, 
            which(colnames(colData) == filterVars$condition)])
        cond_vars = unique(cond_vec)
    }
    Up_Inc = as.matrix(S4Vectors::metadata(se)$Up_Inc)[
        rowData(se)$EventName[
            rowData(se)$EventType %in% c("IR", "MXE", "SE")
        ],
    ]
    Down_Inc = as.matrix(S4Vectors::metadata(se)$Down_Inc)[
        rowData(se)$EventName[
            rowData(se)$EventType %in% c("IR", "MXE", "SE")
        ],
    ]
    IntronDepth = as.matrix(assay(se, "Included"))
    IntronDepth = IntronDepth[rowData$EventType %in% c("IR", "MXE", "SE"),]
    minDepth.Inc = Up_Inc + Down_Inc
    # do not test if depth below threshold
    Up_Inc[minDepth.Inc < minDepth] = IntronDepth[minDepth.Inc < minDepth]
    Down_Inc[minDepth.Inc < minDepth] = IntronDepth[minDepth.Inc < minDepth]

    Excluded = as.matrix(assay(se, "Excluded"))
    Excluded = Excluded[rowData$EventType %in% c("MXE"),]
    Up_Exc = as.matrix(S4Vectors::metadata(se)$Up_Exc)
    Down_Exc = as.matrix(S4Vectors::metadata(se)$Down_Exc)
    minDepth.Exc = Up_Exc + Down_Exc
    # do not test if depth below threshold
    Up_Exc[minDepth.Exc < minDepth] = Excluded[minDepth.Exc < minDepth]    
    Down_Exc[minDepth.Exc < minDepth] = Excluded[minDepth.Exc < minDepth]

    sum_res = rep(0, nrow(se))
    if(use_cond == TRUE) {
        for(cond in cond_vars) {
            Up_Inc.subset = Up_Inc[, which(cond_vec == cond)]
            Down_Inc.subset = Down_Inc[, which(cond_vec == cond)]
            IntronDepth.subset = IntronDepth[, which(cond_vec == cond)]
            Up_Exc.subset = Up_Exc[, which(cond_vec == cond)]
            Down_Exc.subset = Down_Exc[, which(cond_vec == cond)]
            Excluded.subset = Excluded[, which(cond_vec == cond)]

            sum_inc = rowSums(
                abs(log2(Up_Inc.subset + 1) - log2(IntronDepth.subset + 1)) 
                    < maximum &
                abs(log2(Down_Inc.subset + 1) - log2(IntronDepth.subset + 1)) 
                    < maximum
            )
            sum_exc = rowSums(
                abs(log2(Up_Exc.subset + 1) - log2(Excluded.subset + 1)) 
                    < maximum &
                abs(log2(Down_Exc.subset + 1) - log2(Excluded.subset + 1)) 
                    < maximum
            )
            sum_inc = c(sum_inc, rep(ncol(Up_Inc.subset), 
                sum(!(rowData$EventType %in% c("IR", "MXE", "SE")))))
            sum_exc = c(rep(ncol(Up_Inc.subset), 
                sum(rowData$EventType == "IR")),
            sum_exc, rep(ncol(Up_Inc.subset), 
                sum(!(rowData$EventType %in% c("IR", "MXE")))))
            sum = 0.5 * (sum_inc + sum_exc)
            sum_res = sum_res + 
                ifelse(sum * 100 / ncol(Up_Inc.subset) >= usePC, 1, 0)
        }
        n_TRUE = ifelse(
            !is.na(suppressWarnings(as.numeric(filterVars$minCond))), 
            as.numeric(filterVars$minCond), -1)
        if(n_TRUE == -1) n_TRUE = length(cond_vars)
        res = (sum_res >= n_TRUE)
    } else {
        sum_inc = rowSums(
            abs(log2(Up_Inc + 1) - log2(IntronDepth +1)) < filterVars$maximum &
            abs(log2(Down_Inc + 1) - log2(IntronDepth +1)) < filterVars$maximum
        )
        sum_exc = rowSums(
            abs(log2(Up_Exc + 1) - log2(Excluded +1)) < filterVars$maximum &
            abs(log2(Down_Exc + 1) - log2(Excluded +1)) < filterVars$maximum
        )
        sum_inc = c(sum_inc, rep(ncol(Up_Inc), 
            sum(!(rowData$EventType %in% c("IR", "MXE", "SE")))))
        sum_exc = c(
            rep(ncol(Up_Inc), sum(rowData$EventType == "IR")),
            sum_exc, 
            rep(
                ncol(Up_Inc), 
                sum(!(rowData$EventType %in% c("IR", "MXE")))
            )
        )
        sum = 0.5 * (sum_inc + sum_exc)
        res = ifelse(sum * 100 / ncol(Up_Inc) >= usePC, TRUE, FALSE)
    }
    if("EventTypes" %in% names(filterVars)) {
        res[!(rowData(se)$EventType %in% filterVars$EventTypes)] = TRUE
    }
    return(res)
}

################################################################################



#' Convenience function to apply a list of filters to a SummarizedExperiment
#'   object
#'
#' See [runFilter()] for details regarding filters
#' 
#' @param se A SummarizedExperiment object created by `MakeSE()`
#' @param filters A list of filters to apply. Each filter must contain the
#'   elements `filterClass`, `filterType` and `filterVars`. 
#'   See `?runFilter` for details
#' @return A vector of logicals, with `TRUE` indicating events to be retained,
#'   and `FALSE` for events to be filtered out
#' @examples
#' # see ?MakeSE on example code of generating this NxtSE object
#' se = NxtIRF_example_NxtSE()
#' 
#' # Get NxtIRF recommended filters
#' filters = get_default_filters()
#' 
#' # Apply all recommended filters:
#' se.filtered = apply_filters(
#'     se = se, 
#'     filters = filters
#' )
#' @seealso [get_default_filters()], [runFilter()]
#' @md
#' @export
apply_filters <- function(se, filters) {
    # filters are a list of filters to apply on se
    # returns a vector of TRUE / FALSE
    # a filtered se can be made using:
    #       se.filtered = se[apply_filters(se, filters),]
    if(!is(filters, "list")) {
        stop(paste("In apply_filters(),",
            "filters must be a list"
        ), call. = FALSE)
    }
    for(i in seq_len(length(filters))) {
        if(!("filterVars" %in% names(filters[[i]]))) {
            stop(paste("In apply_filters(),",
                "filterVars is missing from filters @ index #", i
            ), call. = FALSE)
        }
        if(!("filterClass" %in% names(filters[[i]]))) {
            stop(paste("In apply_filters(),",
                "filterClass is missing from filters @ index #", i
            ), call. = FALSE)
        }
        if(!("filterType" %in% names(filters[[i]]))) {
            stop(paste("In apply_filters(),",
                "filterType is missing from filters @ index #", i
            ), call. = FALSE)
        }
    }
    if(!is(se, "NxtSE")) {
        stop(paste("In apply_filters(),",
            "se must be a NxtSE object"
        ), call. = FALSE)
    }
    filterSummary = rep(TRUE, nrow(se))
    for(i in seq_len(length(filters))) {
        filterSummary = filterSummary & runFilter(
            filters[[i]]$filterClass,
            filters[[i]]$filterType,
            filters[[i]]$filterVars,
            se
        )
    }
    
    return(filterSummary)
}