# wrappers to R/C++

run_IRFinder_multithreaded = function(
        reference_path = "./Reference", 
        bamfiles = "Unsorted.bam", 
        output_files = "./Sample",
        max_threads = max(parallel::detectCores() - 2, 1),
        Use_OpenMP = TRUE,
        run_featureCounts = FALSE
    ) {
    .validate_reference(reference_path)

    s_bam = normalizePath(bamfiles)
    s_ref = normalizePath(reference_path)
    
    .irfinder_validate_args(s_bam, s_ref, max_threads, output_files)
    
    ref_file = normalizePath(file.path(s_ref, "IRFinder.ref.gz"))


    if(Has_OpenMP() > 0 & Use_OpenMP) {
        n_threads = floor(max_threads)
        n_threads = min(n_threads, length(s_bam))
        IRF_main_multithreaded(ref_file, s_bam, output_files, n_threads)
    } else {
        # Use BiocParallel
        n_rounds = ceiling(length(s_bam) / floor(max_threads))
        n_threads = ceiling(length(s_bam) / n_rounds)

        BPPARAM = BiocParallel::bpparam()
        if(Sys.info()["sysname"] == "Windows") {
            BPPARAM_mod = BiocParallel::SnowParam(n_threads)
            message(paste("Using SnowParam", BPPARAM_mod$workers, "threads"))
        } else {
            BPPARAM_mod = BiocParallel::MulticoreParam(n_threads)
            message(paste("Using MulticoreParam", BPPARAM_mod$workers, 
                "threads"))
        }

        row_starts = seq(1, by = n_threads,
            length.out = n_rounds)
            
        for(i in seq_len(n_rounds)) {
            selected_rows_subset = seq(row_starts[i], 
                min(length(s_bam), row_starts[i] + n_threads - 1)
            )
            BiocParallel::bplapply(selected_rows_subset,
                function(i, IRF_main, s_bam, reference_file, output_files) {
                    IRF_main(s_bam[i], reference_file, output_files[i])
                }, 

                s_bam = s_bam,
                output_files = output_files,
                IRF_main = IRF_main, 
                reference_file = ref_file,
                BPPARAM = BPPARAM_mod
            )
        }
    }
    
    if(run_featureCounts == TRUE) {

        NxtIRF.CheckPackageInstalled("Rsubread", "2.4.0")
        gtf_file <- Get_GTF_file(reference_path)
        
        # determine paired_ness, strandedness, assume all BAMS are the same
        data.list = get_multi_DT_from_gz(
            normalizePath(paste0(output_files[1], ".txt.gz")), 
            c("BAM", "Directionality")
        )
        stats = data.list$BAM
        direct = data.list$Directionality

        if(stats$Value[3] == 0 & stats$Value[4] > 0) {
            paired = TRUE
        } else if(stats$Value[3] > 0 && 
                stats$Value[4] / stats$Value[3] / 1000) {
            paired = TRUE
        } else {
            paired = FALSE
        }
        strand = direct$Value[9]
        if(strand == -1) strand = 2
        
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
        
        if(file.exists(file.path(dirname(output_files[1]), "main.FC.Rds"))) {
            res.old = readRDS(
                file.path(dirname(output_files[1]), "main.FC.Rds"))

            # Check md5 of annotation to show same reference was used
            md5.old = openssl::md5(paste(
                res.old$annotation$GeneID, res.old$annotation$Chr,
                res.old$annotation$Start, res.old$annotation$End, 
                res.old$annotation$Strand, collapse=" "
                ))
            md5 = openssl::md5(paste(
                res$annotation$GeneID, res$annotation$Chr,
                res$annotation$Start, res$annotation$End, 
                res$annotation$Strand, collapse=" "
                ))
            md5.old.stat = openssl::md5(paste(
                res.old$stat$Status, collapse=" "
                ))
            md5.stat = openssl::md5(paste(
                res$stat$Status, collapse=" "
                ))
            if(md5 == md5.old & md5.stat == md5.old.stat) {
                # cbind stats
                new_samples = res$targets[!(res$targets %in% res.old$targets)]
                res$targets = c(res.old$targets, new_samples)

                res$stat = cbind(res.old$stat, res$stat[,new_samples])
                # cbind counts            
                res$counts = cbind(res.old$counts, res$counts[,new_samples])
            }
        }
        
        if(all(c("counts", "annotation", "targets", "stat") %in% names(res))) {
            saveRDS(res, file.path(dirname(output_files[1]), "main.FC.Rds"))
        }
        
    }   

}    



.irfinder_validate_args <- function(s_bam, s_ref, max_threads, output_files) {
    if(max_threads != 1 && max_threads > parallel::detectCores() - 1) {
        stop(paste("In run_IRFinder_multithreaded(), ",
            max_threads, " threads is not allowed for this system"
        ), call. = FALSE)
    }
    if(!all(file.exists(s_bam))) {
        stop(paste("In run_IRFinder_multithreaded(), ",
            paste(
                unique(s_bam[!file.exists(s_bam)]),
                collapse = ""
            ),
            " - files not found"
        ), call. = FALSE)
    }    

    if(!all(dir.exists(dirname(output_files)))) {
        stop(paste("In run_IRFinder_multithreaded(), ",
            paste(
                unique(dirname(
                    output_files[!dir.exists(dirname(output_files))])),
                collapse = ""
            ),
            " - directories not found"
        ), call. = FALSE)
    }

    if(!(length(s_bam) == length(output_files))) {
        stop(paste("In run_IRFinder_multithreaded(), ",
            "Number of output files and bam files must be the same"
        ), call. = FALSE)
    }
    return(TRUE)
}

run_IRFinder_GenerateMapReads = function(genome.fa = "", out.fa, 
    read_len = 70, read_stride = 10, error_pos = 35) {
    return(
        IRF_GenerateMappabilityReads(normalizePath(genome.fa), 
            file.path(normalizePath(dirname(out.fa)), basename(out.fa)),
            read_len = read_len, 
            read_stride = read_stride, 
            error_pos = error_pos)
    )
}

run_IRFinder_MapExclusionRegions = function(bamfile = "", output_file, 
        threshold = 4, includeCov = FALSE) {
    s_bam = normalizePath(bamfile)
    if(!file.exists(s_bam)) {
        stop(paste("In run_IRFinder_MapExclusionRegions(),",
            s_bam, "does not exist"
        ), call. = FALSE)
    }
    return(
        IRF_GenerateMappabilityRegions(s_bam, 
            output_file,
            threshold = threshold,
            includeCov = includeCov)
    )
}

run_Gunzip = function(infile = "", outfile) {
    file_to_read = normalizePath(infile)
    if(!file.exists(file_to_read)) {
        stop(paste("In run_Gunzip(),",
            file_to_read, "does not exist"
        ), call. = FALSE)
    }
    if(!dir.exists(dirname(outfile))) {
        stop(paste("In run_Gunzip(),",
            dirname(outfile), "does not exist"
        ), call. = FALSE)
    }
    IRF_gunzip(file_to_read, outfile)
}

get_multi_DT_from_gz = function(infile = "", 
        block_headers = c("Header1", "Header2")) {
    file_to_read = normalizePath(infile)
    if(!file.exists(file_to_read)) {
        stop(paste("In get_multi_DT_from_gz(),",
            file_to_read, "does not exist"
        ))
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

#' Calls NxtIRF's C++ function to retrieve coverage
#'
#' This function returns an RLE or RLEList containing coverage data from the
#' given COV file
#' @param file The file name of the COV file
#' @param seqname Either blank, or a character string denoting the chromosome 
#'  name
#' @param start The 0-based start coordinate 
#' @param end The 0-based end coordinate
#' @param strand An integer denoting ths strand: "+" = 0, "-" = 1, "*" = 2
#' @return If seqname is left as "", returns an RLEList of the whole BAM file.
#'   If seqname and coordinates are given, returns an RLE containing the
#'   chromosome coordinate. Coordinates outside the given range will be set to 0
#' @examples
#' se <- NxtIRF_example_NxtSE()
#' 
#' cov_file <- covfile(se)[1]
#'
#' cov <- GetCoverage(cov_file, seqname = "chrZ", 
#'   start = 10000, end = 20000,
#'   strand = 2
#' )
#' @export
GetCoverage <- function(file, seqname = "", start = 0, end = 0, strand = 2) {
    if(!(as.numeric(strand) %in% c(0,1,2))) {
        stop(paste("In GetCoverage(),",
            "Invalid strand. Must be either 0 (+), 1 (-) or 2(*)"
        ), call. = FALSE)
    }
    if(!is.numeric(start) || !is.numeric(end) || 
            (as.numeric(start) > as.numeric(end)) || 
            as.numeric(end) == 0) {
        stop(paste("In GetCoverage(),",
            "Null or negative regions not allowed"
        ), call. = FALSE)
    }

    if(seqname == "") {
        raw_list = IRF_RLEList_From_Cov(normalizePath(file), strand)
        final_list = list()
        if(length(raw_list) > 0) {
            for(i in seq_len(length(raw_list))) {
                final_list[[i]] = S4Vectors::Rle(
                    raw_list[[i]]$values, raw_list[[i]]$length
                )
            }
        } else {
            return(NULL)
        }
        final_RLE = as(final_list, "RleList")
        names(final_RLE) = names(raw_list)
        return(final_RLE)
    } else if(end == 0) {
        raw_RLE = IRF_RLE_From_Cov(
            normalizePath(file), as.character(seqname), 
            0,0, as.numeric(strand)
        )
        final_RLE = S4Vectors::Rle(raw_RLE$values, raw_RLE$length)
    } else {
        raw_RLE = IRF_RLE_From_Cov(
            normalizePath(file), as.character(seqname), 
            round(as.numeric(start)), round(as.numeric(end)), 
            as.numeric(strand)
        )
        final_RLE = S4Vectors::Rle(raw_RLE$values, raw_RLE$length)
    }
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