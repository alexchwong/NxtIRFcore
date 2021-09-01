#' NxtIRF package for IRFinder-based differential Alternative Splicing
#' and Intron Retention analysis
#' 
#' The NxtIRF package analyses RNA-seq datasets (aligned as BAM files using 
#' splice aware genome aligners such as STAR). It fully incorporates the
#' IRFinder algorithm, porting its C++ code into the R platform using RCPP.
#' @details
#' Enhancements include:
#'
#' \itemize{
#' \item One-step reference generation from user-supplied or AnnotationHub 
#'   genome and gene annotations;
#' \item Multi-threaded support for core IRFinder algorithm to simultaneously
#'   analyse multiple BAM files;
#' \item Streamlined integration of IRFinder results from large datasets;
#' \item Easy generation of experimental designs and data handling using
#'   NxtSE which is a SummarizedExperiment object extended for NxtIRF-based
#'   analysis;
#' \item Simplified workflow for filtering low-abundance splicing events and
#'   limma- or DESeq2- based differential alternative splicing events analysis;
#' \item Interactive workflow and visualisation tools using shinydashboard;
#' \item Advanced RNA-seq coverage visualisation, including the ability to 
#'   combine RNA-seq coverage of multiple samples using advanced normalisation 
#'   methods across samples grouped by conditions;
#' }
#' 
#' The main functions are:
#'
#' \itemize{
#' \item \code{\link{BuildReference}} - Prepares genome and gene annotation
#'   references (FASTA/TwoBit, GTF files), and synthesises the NxtIRF reference
#'   for the IRFinder engine and NxtIRF-based analysis
#' \item \code{\link{IRFinder}} - Runs the IRFinder C++ routine to analyse
#'   single or multiple BAM files using the NxtIRF/IRFinder reference.
#' \item \code{\link{CollateData}} - Combines multiple IRFinder outputs
#'   into one unified data structure
#' \item \code{\link{MakeSE}} - Constructs a \code{NxtSE} 
#'   object (see \code{\link{NxtSE-methods}}) to contain IRFinder output 
#'   produced by \code{\link{CollateData}}
#' \item \code{\link{apply_filters}} - Easily specify and apply a list of
#'   filters to exclude low-abundance splicing events from downstream analysis
#' \item \code{\link{limma_ASE}}, \code{\link{DESeq_ASE}} - perform
#'   differential alternate splice event (ASE) analysis on a filtered NxtSE 
#    object using limma or DESeq2
#' \item \code{\link{make_diagonal}}, \code{\link{make_matrix}}: Generate
#'   data to produce scatter plots or heatmaps of ASE events
#' \item \code{\link{Plot_Coverage}}: Generate RNA-seq coverage plots of
#'   individual samples or across samples grouped by user-specified conditions
#' }
#' 
#' See \code{vignette("NxtIRFcore")} for worked examples on how to use NxtIRF
#' 
#' @author Alex Wong, Ulf Schmitz, William Ritchie
#' 
#' @docType package
#' @name NxtIRFcore-package
#' @aliases NxtIRFcore-package
#' @keywords package
#' @md
NULL

#' Builds reference files used by IRFinder / NxtIRF.
#'
#' @description
#' These function builds the reference required by the IRFinder engine, as well
#' as access-ready refined splice annotation data for NxtIRF. See details below.
#' @details
#' `GetReferenceResource()` processes the files, downloads resources from
#' web links or from `AnnotationHub()`, and saves a local copy in the "resource"
#' subdirectory within the given `reference_path`\cr\cr
#' `BuildReference()` runs `GetReferenceResource()` if resources are not 
#' saved locally (i.e. `GetReferenceResource()` is not already run). Then,
#' it creates the NxtIRF / IRFinder references.\cr\cr
#' NB: the parameters `fasta` and `gtf` can be omitted in `BuildReference()` if
#' `GetReferenceResource()` is already run. See examples below.\cr\cr
#' The NxtIRF reference can be created using either:\cr\cr
#' 1. User-supplied FASTA and GTF file. This can be a file path, or a web link
#'   (e.g. 'http://', 'https://' or 'ftp://'). Use `fasta_file` and `gtf_file`
#'    to specify the files or web paths to use.\cr\cr
#' 2. AnnotationHub genome and gene annotation (Ensembl): supply the names of
#'    the genome sequence and gene annotations  \cr\cr
#'
#' @param reference_path The directory to store the reference files
#' @param fasta The file path or web link to the user-supplied genome
#'   FASTA file. Alternatively, the name of the AnnotationHub record containing
#'   the genome resource (for Ensembl, this is a TwoBit file).
#' @param gtf The file path or web link  to the user-supplied transcript 
#'   GTF file (or gzipped GTF file). Alternatively, the name of the
#'   AnnotationHub record containing the transcript GTF file
#' @param generate_mappability_reads (default FALSE) Whether reads should be
#'   generated for the purpose of Mappability calculations. See 
#'   \code{\link{Mappability-methods}}
#' @param convert_chromosome_names (Optional) A 2-column data frame containing 
#'   chromosome name conversions. The first column lists the chromosome names 
#'   of the source reference files, and the second column gives the desired 
#'   chromosome names). See example below, or refer to
#'   <https://github.com/dpryan79/ChromosomeMappings> for a list of
#'   chromosome conversion resources. If chromosome conversion occurs, only
#'   annotations that match the source chromosome names will be included in
#'   the final NxtIRF reference (i.e. records with non-matching chromosome
#'   names will be discarded)
#' @param overwrite_resource (default FALSE) `BuildReference()` will first 
#'   run `GetReferenceResource()` before generating the NxtIRF reference.
#'   If the genome TwoBit and gene annotation GTF files (generated by 
#'   `GetReferenceResource()`) are found in the "resource" subdirectory, 
#'   `BuildReference()` will use these resources (and not overwrite these),
#'   unless `overwrite_resource` is set to TRUE.
#' @param ... In GetReferenceResource(), if `generate_mappability_reads` is set
#'   to `TRUE`, additional arguments to parse to `Mappability_GenReads()`
#'   See \link{Mappability-methods}.
#' @param genome_type Allows `BuildReference()` to select default 
#'   `nonPolyARef` and `MappabilityRef` for selected genomes. Allowed options 
#'   are: 'hg38', 'hg19', 'mm9', 'mm10'.
#'   Note that the user can still set `nonPolyARef` and `MappabilityRef` values
#'   to override the default files. Alternatively, use `GetNonPolyARef()` and
#'   `GetMappabilityRef()` to retrieve NxtIRF-supplied default files
#'   for the supported chromosomes (see examples below).
#' @param nonPolyARef A BED file (3 unnamed columns containing chromosome, 
#'   start and end coordinates) of regions defining known non-polyadenylated 
#'   transcripts. This file is used for QC analysis of IRFinder-processed files 
#'   to measure Poly-A enrichment quality of samples. Leave blank to not use a 
#'   `nonPolyARef` file (or to use default - see `genome_type`).
#' @param MappabilityRef A BED file (3 unnamed columns containing chromosome, 
#'   start and end coordinates) of poorly-mapped regions due to repeat elements 
#'   in the genome. We recommend using the default Mappability files supplied 
#'   (see `genome_type`). Alternately, this reference can be generated by 
#'   running `Mappability_GenReads()` on the genome sequence, followed by 
#'   alignment of the produced fasta file to an aligner of choice (e.g. STAR, 
#'   HISAT2). The aligned sequences (as BAM file) should then be analysed using 
#'   `Mappability_CalculateExclusions()`, which will provide the Mappability 
#'   file to be used here.
#' @param BlacklistRef A BED file (3 unnamed columns containing chromosome, 
#'   start and end coordinates) of regions to be otherwise excluded from IR 
#'   analysis. Leave blank to not use a `BlacklistRef` file.
#' @param UseExtendedTranscripts (default TRUE) Should non-protein-coding 
#'   transcripts such as anti-sense and lincRNAs be included in searching
#'   for IR / AS events? Setting `FALSE` (vanilla IRFinder) will exclude 
#'   transcripts other than `protein_coding` and 
#'   `processed_transcript` transcripts from IR analysis.
#' @param n_threads The number of threads used to generate the STAR reference
#'   and mappability calculations. Multi-threading is not used for NxtIRF
#'   reference generation (but multiple cores are utilised by data.table
#'   and fst packages automatically, where available). See [STAR-methods]
#' @return Nothing. The created reference will be written to the given 
#'   directory. 
#'   This includes:
#' * `reference_path`/settings.Rds: An RDS file containing parameters used
#'   to generate the NxtIRF reference
#' * `reference_path`/IRFinder.ref.gz: A gzipped text file containing collated 
#'   IRFinder references to be used as input for the IRFinder analysis
#' * `reference_path`/fst/: Contains fst files for subsequent easy access to 
#'   NxtIRF generated references
#' * `reference_path`/resource/genome.2bit: Contains a TwoBitFile generated 
#'   by this function for easy subsequent access to the genome. This step is
#'   skipped if an AnnotationHub resource is used and convert_chromosome_names
#'   is not used, as the TwoBit file is identical to the cached AnnotationHub
#'   resource.
#' * `reference_path`/resource/transcripts.gtf.gz: Contains a copy of the GTF
#'   file used by this reference. This is used by featureCounts option in
#'   `IRFinder()` 
#' @examples
#' # Reference generation from NxtIRF's example mock genome
#' 
#' GetReferenceResource(
#'     reference_path = file.path(tempdir(), "Reference"),
#'     fasta = mock_genome(), gtf = mock_gtf()
#' )
#' BuildReference(
#'     reference_path = file.path(tempdir(), "Reference")
#' )
#' 
#' # Gets path to the Non-PolyA BED file for hg19
#'
#' GetNonPolyARef("hg19")
#'
#' # Gets path to the Mappability Exclusion BED file for mm10 (from a tempdir())
#'
#' GetMappabilityRef("mm10")
#'
#' \dontrun{
#' # Reference generation from user supplied FASTA and GTF files
#' 
#' GetReferenceResource(
#'     reference_path = "./Reference_user",
#'     fasta = "genome.fa", gtf = "transcripts.gtf"
#' )
#' BuildReference(
#'     reference_path = "./Reference_user",
#'     genome_type = "hg38" 
#' )
#'
#' # Reference generation from Ensembl's FTP links:
#' 
#' FTP = "ftp://ftp.ensembl.org/pub/release-94/"
#' GetReferenceResource(
#'     reference_path = "./Reference_FTP",
#'     fasta = paste0(FTP, "fasta/homo_sapiens/dna/",
#'         "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"), 
#'     gtf = paste0(FTP, "gtf/homo_sapiens/",
#'         "Homo_sapiens.GRCh38.94.chr.gtf.gz")
#' )
#' BuildReference(reference_path = "./Reference_FTP",
#'     genome_type = "hg38" 
#' )
#' 
#' # Get AnnotationHub record names for Ensembl release-94:
#' 
#' ah = AnnotationHub::AnnotationHub()
#' ah_r94 = AnnotationHub::query(ah, c("Homo Sapiens", "release-94"))
#' ah_r94
#' # AnnotationHub with 9 records
#' # # snapshotDate(): 2021-04-12
#' # # $dataprovider: Ensembl
#' # # $species: Homo sapiens
#' # # $rdataclass: TwoBitFile, GRanges
#' # # additional mcols(): taxonomyid, genome, description, coordinate_1_based, 
#' # # maintainer, rdatadateadded,
#' # #   preparerclass, tags, rdatapath, sourceurl, sourcetype
#' # # retrieve records with, e.g., 'object[["AH64628"]]'
#' # 
#' #             title                                           
#' #   AH64628 | Homo_sapiens.GRCh38.94.abinitio.gtf
#' #   AH64629 | Homo_sapiens.GRCh38.94.chr.gtf
#' #   AH64630 | Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf
#' #   AH64631 | Homo_sapiens.GRCh38.94.gtf
#' #   AH65744 | Homo_sapiens.GRCh38.cdna.all.2bit
#' #   AH65745 | Homo_sapiens.GRCh38.dna.primary_assembly.2bit
#' #   AH65746 | Homo_sapiens.GRCh38.dna_rm.primary_assembly.2bit
#' #   AH65747 | Homo_sapiens.GRCh38.dna_sm.primary_assembly.2bit
#' #   AH65748 | Homo_sapiens.GRCh38.ncrna.2bit
#' 
#' # Reference generation from AnnotationHub's Ensembl release-94:
#' 
#' GetReferenceResource(
#'     reference_path = "./Reference_AH",
#'     fasta = "AH65745", 
#'     gtf = "AH64631"
#' )
#' BuildReference(
#'     reference_path = "./Reference_AH",
#'     genome_type = "hg38" 
#' )
#' 
#' # AnnotationHub's Ensembl release-94, converting Ensembl chromosomes
#' # to UCSC style:
#' 
#' chrom.df = GenomeInfoDb::genomeStyles()$Homo_sapiens
#' 
#' GetReferenceResource(
#'     reference_path = "./Reference_UCSC",
#'     fasta = "AH65745", 
#'     gtf = "AH64631",
#'     convert_chromosome_names = chrom.df[, c("Ensembl", "UCSC")]
#' )
#' BuildReference(
#'     reference_path = "./Reference_UCSC",
#'     genome_type = "hg38", 
#'     convert_chromosome_names = chrom.df[, c("Ensembl", "UCSC")]
#' )
#'
#' # One-step generation of NxtIRF and STAR references, using 4 threads.
#' # NB: requires a linux-based system with STAR installed.
#'
#' BuildReference_Full(
#'     reference_path = "./Reference_with_STAR",
#'     fasta = "genome.fa", gtf = "transcripts.gtf"
#'     genome_type = "hg38",
#'     n_threads = 4
#' )
#'
#' }
#' @seealso 
#' [Mappability-methods]\cr\cr
#' [STAR-methods]\cr\cr
#' \link[AnnotationHub]{AnnotationHub}\cr\cr
#' @name BuildReference
#' @md
NULL

#' @describeIn BuildReference Processes / downloads a copy of the genome
#' and gene annotations and stores this in the "resource" subdirectory
#' of the given reference path
#' @export
GetReferenceResource <- function(
        reference_path = "./Reference",
        fasta, gtf,
        generate_mappability_reads = FALSE,
        convert_chromosome_names = NULL,
        overwrite_resource = FALSE,
        ...
) {
    .validate_path(reference_path, subdirs = "fst")

    chromosomes = .convert_chromosomes(convert_chromosome_names)
    reference_data = .get_reference_data(reference_path = reference_path,
        fasta = fasta, gtf = gtf, 
        chromosomes = chromosomes, 
        overwrite_resource = overwrite_resource
    )
    map_reads = file.path(normalizePath(reference_path), 
        "Mappability", "Reads.fa")
    if(generate_mappability_reads &&
            (overwrite_resource || !file.exists(map_reads))) {
        Mappability_GenReads(reference_path, ...)
    }
}

#' @describeIn BuildReference First calls \code{GetReferenceResource()}
#' (if required). Afterwards creates the NxtIRF reference in the
#' given reference path
#' @export
BuildReference <- function(
        reference_path = "./Reference",
        fasta, gtf,
        convert_chromosome_names = NULL,
        overwrite_resource = FALSE, genome_type, 
        nonPolyARef = "", MappabilityRef = "", BlacklistRef = "", 
        UseExtendedTranscripts = TRUE
    ) {
    if(missing(genome_type)) genome_type = ""
    .validate_path(reference_path, subdirs = "fst")
    extra_files <- .fetch_genome_defaults(reference_path,
        genome_type, nonPolyARef, MappabilityRef, BlacklistRef
    )
    N <- 8
    dash_progress("Reading Reference Files", N)
    chromosomes = .convert_chromosomes(convert_chromosome_names)
    reference_data = .get_reference_data(reference_path = reference_path,
        fasta = fasta, gtf = gtf, 
        chromosomes = chromosomes, 
        overwrite_resource = overwrite_resource
    )
    dash_progress("Processing gtf file", N)
    reference_data$gtf_gr = .validate_gtf(reference_data$genome, 
        reference_data$gtf_gr)
    reference_data$gtf_gr = .fix_gtf(reference_data$gtf_gr)
    .process_gtf(reference_data$gtf_gr, reference_path)
    reference_data$gtf_gr = NULL # To save memory, remove original gtf
    gc()

    dash_progress("Processing introns", N)
    .process_introns(reference_path, reference_data$genome, 
        UseExtendedTranscripts)
    dash_progress("Generating IRFinder Reference", N)
    .gen_irf(reference_path, extra_files, reference_data$genome, chromosomes)
    gc()
    # Annotate IR-NMD
    dash_progress("Annotating IR-NMD", N)
    dash_withProgress(message = 'Determining NMD Transcripts', value = 0, {
        .gen_nmd(reference_path, reference_data$genome)
    })
    # Annotating Alternative Splicing Events
    dash_progress("Annotating Splice Events", N)
    .gen_splice(reference_path, reference_data$genome)
    if(file.exists(file.path(reference_path, "fst", "Splice.fst"))) {
        dash_progress("Translating AS Peptides", N)
        .gen_splice_proteins(reference_path, reference_data$genome)    
        message("Splice Annotations finished\n")
    } else {
        dash_progress("No alternate splicing events detected in this reference",
            N)
    }
    message("Reference build finished")
    dash_progress("Reference build finished", N)

    # Update settings.Rds only after everything is finalised
    settings.list = readRDS(file.path(reference_path, "settings.Rds"))
    settings.list$genome_type = genome_type
    settings.list$nonPolyARef = nonPolyARef
    settings.list$MappabilityRef = MappabilityRef
    settings.list$BlacklistRef = BlacklistRef
    settings.list$UseExtendedTranscripts = UseExtendedTranscripts
    
    settings.list$BuildVersion = buildref_version
    
    saveRDS(settings.list, file.path(reference_path, "settings.Rds"))
}

#' @describeIn BuildReference Returns the path to the BED file containing
#'   coordinates of known non-polyadenylated transcripts for genomes 
#'   \code{hg38}, \code{hg19}, \code{mm10} and \code{mm9}, 
#' @export
GetNonPolyARef <- function(genome_type) {
    if (genome_type == "hg38") {
        nonPolyAFile <- system.file(
            "extra-input-files/Human_hg38_nonPolyA_ROI.bed",
            package = "NxtIRFcore"
        )
    } else if (genome_type == "hg19") {
        nonPolyAFile <- system.file(
            "extra-input-files/Human_hg19_nonPolyA_ROI.bed",
            package = "NxtIRFcore"
        )
    } else if (genome_type == "mm10") {
        nonPolyAFile <- system.file(
            "extra-input-files/Mouse_mm10_nonPolyA_ROI.bed",
            package = "NxtIRFcore"
        )
    } else if (genome_type == "mm9") {
        nonPolyAFile <- system.file(
            "extra-input-files/Mouse_mm9_nonPolyA_ROI.bed",
            package = "NxtIRFcore"
        )
    } else {
        nonPolyAFile <- ""
    }
    return(nonPolyAFile)
}

#' @describeIn BuildReference Returns the path to the NxtIRF-provided
#'   Mappability Exclusion coordinates of lowly-mappable regions for genomes 
#'   \code{hg38}, \code{hg19}, \code{mm10} and \code{mm9}, 
#' @export
GetMappabilityRef <- function(genome_type) {
    if (!(genome_type %in% c("hg38", "hg19", "mm9", "mm10"))) {
        mapfile <- ""
    } else {
        mapfile <- get_mappability_exclusion(genome_type)
    }
    return(mapfile)
}

Get_Genome <- function(reference_path, validate = TRUE) {
    if(validate) .validate_reference(reference_path)
    if(file.exists(file.path(reference_path, "resource", "genome.2bit"))) {
        return(rtracklayer::TwoBitFile(
            file.path(reference_path, "resource", "genome.2bit")))
    } else if(file.exists(file.path(reference_path, "settings.Rds"))){
        settings = readRDS(file.path(reference_path, "settings.Rds"))
        genome = .fetch_AH(settings$ah_genome, 
            rdataclass = "TwoBitFile")
    } else {
        .log("In Get_Genome, invalid reference_path supplied")
    }
    return(genome)
}

Get_GTF_file <- function(reference_path) {
    .validate_reference(reference_path)
    if(file.exists(file.path(reference_path, 
            "resource", "transcripts.gtf.gz"))) {
        return(file.path(reference_path, "resource", "transcripts.gtf.gz"))
    } else {
        .log("In Get_GTF_file, invalid reference_path supplied")
    }
}

################################################################################
# Validation functions

.validate_genome_type <- function(genome_type) {
    if(genome_type != "") return(TRUE)
    .log(paste("In BuildReference(),",
        "genome_type not specified.",
        "This should be either one of 'hg38', 'hg19', 'mm10', 'mm9', or",
        "'other'. If 'other', please provide a nonPolyARef file or leave",
        "blank to omit polyA profiling."
    ))
}

.validate_path <- function(reference_path, subdirs = NULL) {
    if({
        reference_path != "" &&
        tryCatch(
            ifelse(
                normalizePath(dirname(reference_path)) != "", TRUE, TRUE
            ),
            error = function(e) FALSE
        )
    }) {
        # continue
    } else {
        .log(paste("Error in 'reference_path',",
            paste0("base path of '", reference_path, "' does not exist")
        ))
    }

    base <- normalizePath(dirname(reference_path))
    if (!dir.exists(file.path(base, basename(reference_path)))) {
        dir.create(file.path(base, basename(reference_path)))
    }
    if(!is.null(subdirs)) {
        for(subdir in subdirs) {
            if (!dir.exists(file.path(base, basename(reference_path), subdirs))) {
                dir.create(file.path(base, basename(reference_path), subdirs))
            }
        }
    }

    return(file.path(base, basename(reference_path)))
}

.validate_reference_resource <- function(reference_path, from = "") {
    ref <- normalizePath(reference_path)
    from_str = ifelse(from == "", "", 
        paste("In function", from, ":"))
    if (!dir.exists(ref)) {
        .log(paste(from_str,
            "in reference_path =", reference_path,
            ": this path does not exist"))
    }
    if (!file.exists(file.path(ref, "settings.Rds"))) {
        .log(paste(from_str,
            "in reference_path =", reference_path,
            ": settings.Rds not found"))
    }    
    settings.list <- readRDS(file.path(ref, "settings.Rds"))
    if (!("BuildVersion" %in% names(settings.list)) ||
            settings.list[["BuildVersion"]] < buildref_version) {
        .log(paste(from_str,
            "in reference_path =", reference_path,
            "NxtIRF reference is earlier than current version",
            buildref_version))
    }     
}

.validate_reference <- function(reference_path, from = "") {
    ref <- normalizePath(reference_path)
    from_str = ifelse(from == "", "", 
        paste("In function", from, ":"))
    if (!dir.exists(ref)) {
        .log(paste(from_str,
            "in reference_path =", reference_path,
            ": this path does not exist"))
    }
    if (!file.exists(file.path(ref, "settings.Rds"))) {
        .log(paste(from_str,
            "in reference_path =", reference_path,
            ": settings.Rds not found"))
    }    
    if (!file.exists(file.path(ref, "IRFinder.ref.gz"))) {
        .log(paste(from_str,
            "in reference_path =", reference_path,
            ": IRFinder.ref.gz not found"))
    }
    settings.list <- readRDS(file.path(ref, "settings.Rds"))
    if (!("BuildVersion" %in% names(settings.list)) ||
            settings.list[["BuildVersion"]] < buildref_version) {
        .log(paste(from_str,
            "in reference_path =", reference_path,
            "NxtIRF reference is earlier than current version",
            buildref_version))
    }     
}

.fetch_genome_defaults <- function(reference_path, 
        genome_type, nonPolyARef = "", 
        MappabilityRef = "", BlacklistRef = "") {
    if(!is_valid(nonPolyARef)) {
        nonPolyAFile <- GetNonPolyARef(genome_type)
        nonPolyAFile <- .parse_valid_file(nonPolyAFile, "non-polyA reference")        
    } else {
        nonPolyAFile <- .parse_valid_file(nonPolyARef, "non-polyA reference") 
    }
    map_file = file.path(normalizePath(reference_path), "Mappability",
        "MappabilityExclusion.bed.txt")
    if(is_valid(MappabilityRef)) {
        MappabilityFile <- .parse_valid_file(MappabilityRef)
    } else if(file.exists(map_file)) {
        MappabilityFile <- .parse_valid_file(map_file)
    } else if(genome_type %in% c("hg38", "hg19", "mm9", "mm10")) {
        MappabilityFile <- .parse_valid_file(GetMappabilityRef(genome_type))
    } else {
        MappabilityFile <-
            .parse_valid_file(MappabilityRef, "Mappability reference")        
    }
    BlacklistFile <-
        .parse_valid_file(BlacklistRef, "Blacklist exclusion")
        
    final <- list(
        nonPolyAFile = nonPolyAFile, MappabilityFile = MappabilityFile,
        BlacklistFile = BlacklistFile
    )
    return(final)
}
.convert_chromosomes <- function(convert_chromosome_names) {
    if(is.null(convert_chromosome_names)) return(NULL)
    df = as.data.frame(convert_chromosome_names)
    df = df[!duplicated(df[,1]),]
    df = df[!duplicated(df[,2]),]
    df = df[as.vector(df[,1]) != "",]
    df = df[as.vector(df[,2]) != "",]
    df = as.data.table(df)
    colnames(df) = c("Original", "New")
    return(df)
}

################################################################################
# Sub

.get_reference_data <- function(reference_path, fasta, gtf, 
        chromosomes = NULL, overwrite_resource = FALSE) {
    # Validate arguments:
    .get_reference_data_validate(fasta, gtf, reference_path, 
        chromosomes = NULL, overwrite_resource = FALSE)
    # If resources already exist in 'resource', then recall these:
    if(     !( is_valid(fasta) | is_valid(gtf) ) & 
            !overwrite_resource) {
        resource_path <- file.path(reference_path, "resource")
        # settings.file <- file.path(reference_path, "settings.Rds")
        # settings.list <- readRDS(settings.file)
        genome <- .fetch_fasta(reference_path,
            convert_chromosome_names = chromosomes,
            overwrite = overwrite_resource)
        gtf_gr <- .fetch_gtf(
            gtf = file.path(resource_path, "transcripts.gtf.gz"),
            reference_path = reference_path, 
            convert_chromosome_names = chromosomes,
            overwrite = overwrite_resource
        )
        # settings.list$chromosomes = chromosomes
        # settings.list$BuildVersion = buildref_version
        # saveRDS(settings.list, file.path(reference_path, "settings.Rds"))
    } else if(!( is_valid(fasta) | is_valid(gtf) )){
        .log("Both fasta and gtf are required if overwrite_resource = TRUE")
    } else {
        fasta_use <- gtf_use <- ""
        ah_genome_use <- ah_gtf_use <- ""
        if(.is_AH_pattern(fasta)) {
            ah_genome_use <- fasta
        } else {
            fasta_use <- fasta
        }
        if(.is_AH_pattern(gtf)) {
            ah_gtf_use <- gtf
        } else {
            gtf_use <- gtf
        }
        # Check web links are valid
        test_urls <- c(fasta, gtf)
        for(url in test_urls) {
            if(any(startsWith(url, c("http", "ftp")))) {
                ret = .check_if_url_exists(url)
                if(!ret) {
                    .log(paste(url, "is not accessible at this time.",
                    "Please try again later"))
                }
            }
        }
        genome <- .fetch_fasta(
            reference_path = reference_path,
            fasta = fasta_use, ah_genome = ah_genome_use,
            convert_chromosome_names = chromosomes,
            overwrite = overwrite_resource
        )
        gtf_gr <- .fetch_gtf(
            gtf = gtf_use, ah_transcriptome = ah_gtf_use,
            reference_path = reference_path, 
            convert_chromosome_names = chromosomes,
            overwrite = overwrite_resource
        )
        # Save Resource details to settings.Rds:
        settings.list <- list(fasta_file = fasta_use, gtf_file = gtf_use,
            ah_genome = ah_genome_use, ah_transcriptome = ah_gtf_use,
            chromosomes = chromosomes, reference_path = reference_path
        )
        settings.list$BuildVersion = buildref_version
        saveRDS(settings.list, file.path(reference_path, "settings.Rds"))
    }
    settings.list = readRDS(file.path(reference_path, "settings.Rds"))
    final = list(
        genome = genome, gtf_gr = gtf_gr
    )
    return(final)
}

.is_AH_pattern <- function(word) {
    if(substr(word, 1, 2) == "AH" && !file.exists(word)) return(TRUE)
    return(FALSE)
}

################################################################################

.get_reference_data_validate <- function(
        fasta, gtf, reference_path, 
        chromosomes = NULL, overwrite_resource = FALSE) {

    if(     !( is_valid(fasta) | is_valid(gtf) ) & 
            !overwrite_resource) {
        resource_path <- file.path(reference_path, "resource")
        settings.file <- file.path(reference_path, "settings.Rds")
        if(!file.exists(settings.file)) {
            .log(paste("Invalid reference path:", reference_path))     
        }
        settings.list <- readRDS(settings.file)
        if(!file.exists(file.path(resource_path, "genome.2bit"))) {
            .log(paste("Genome could not be found inside", reference_path))
        }
        if(!file.exists(file.path(resource_path, "transcripts.gtf.gz"))) {
            .log(paste("Gene annotations (GTF) could not be found inside", 
                reference_path))
        }
    } else {
        if(!is_valid(fasta) | !is_valid(gtf)) {
            .log("Both fasta and gtf are required")      
        }
    }
}

################################################################################

.fetch_fasta <- function(reference_path = "./Reference",
        fasta = "", ah_genome = "", verbose = TRUE,
        convert_chromosome_names = NULL,
        exclude_scaffolds = TRUE,
        overwrite = FALSE) {
    if (ah_genome != "") {
        fasta_file = ""
        genome <- .fetch_fasta_ah(ah_genome, 
            exclude_scaffolds = exclude_scaffolds,
            verbose = TRUE)
    } else if(fasta == "") {
        file.2bit = file.path(reference_path, "resource", "genome.2bit")
        if(!file.exists(file.2bit)) {
            .log(paste("Could not find", file.2bit))
        }
        genome <- TwoBitFile(file.2bit)
        return(genome)
    } else {
        fasta_file <- .fetch_fasta_file_validate(fasta)
        genome <- .fetch_fasta_file(fasta_file)
        if(exclude_scaffolds) {
            genome = genome[!grepl("scaffold", names(genome))]
        }
    }
    # Convert chromosome names to appropriate
    genome <- .fetch_fasta_convert_chrom(genome, convert_chromosome_names)

    # Save local copy of FASTA
    # .fetch_fasta_save_fasta(genome, reference_path, overwrite)
    .fetch_fasta_save_2bit(genome, reference_path, overwrite)    
    gc()
    message("Connecting to genome TwoBitFile...", appendLF = FALSE)
        genome_2bit <- Get_Genome(reference_path, validate = FALSE)
    message("done\n")
    return(genome_2bit)
}

.fetch_fasta_ah <- function(ah_genome, exclude_scaffolds = TRUE, 
        verbose = verbose) {
    if(substr(ah_genome, 1, 2) != "AH") {
        .log("Given genome AnnotationHub reference is incorrect")
    }
    genome <- .fetch_AH(ah_genome, verbose = verbose, 
        rdataclass = "TwoBitFile", 
        as_DNAStringSet = exclude_scaffolds, 
        exclude_scaffolds = exclude_scaffolds)

}

.fetch_fasta_file_validate <- function(fasta) {
    fasta_file = .parse_valid_file(fasta)
    if(!file.exists(fasta_file)) {
        .log(paste("Given genome fasta file", fasta_file, "not found"))
    }
    return(fasta_file)
}

.fetch_fasta_file <- function(fasta_file) {
    message("Importing genome into memory...", appendLF = FALSE)
        genome <- Biostrings::readDNAStringSet(fasta_file)
    message("done")
    return(genome)
}

.fetch_fasta_convert_chrom <- function(genome, convert_chromosome_names) {
    names(genome) = tstrsplit(names(genome), split = " ")[[1]]
    if(!is.null(convert_chromosome_names)) {
        converter = data.table(Original = 
                tstrsplit(names(genome), split=" ")[[1]])
        converter[convert_chromosome_names,
            on = "Original", c("NewName") := get("i.New")]
        # converter[is.na(get("NewName")), c("NewName") := get("Original")]
        converter = converter[!is.na(get("NewName"))]
        chrOrder = converter$NewName
        genome = genome[converter$Original]
        names(genome) <- chrOrder
    }
    return(genome)
}

.fetch_fasta_convert_chrom_species <- function(genome, species) {
    db = GenomeInfoDb::genomeStyles()
    use_species = NULL
    for(test_species in names(db)) {
        if(grepl(sub("_", " ", test_species), species, ignore.case = TRUE)) 
            use_species = test_species
    }
    if(is.null(use_species)) return(genome)

    .log(paste("Species detected", use_species), type = "message")

    infer_style <- GenomeInfoDb::seqlevelsStyle(genome)
    if(length(infer_style) >= 1) {
        infer_style = infer_style[1]
        .log(paste("Chrom style detected", infer_style), type = "message")
        df = db[[use_species]]
        return(.fetch_fasta_convert_chrom(genome,
            .convert_chromosomes(data.frame(
                old = df[, infer_style],
                new = df[, infer_style],
                stringsAsFactors = FALSE
            ))
        ))
    } else {
        return(genome)
    }
}

.fetch_fasta_save_fasta <- function(genome, reference_path, overwrite) {
    if (!dir.exists(file.path(reference_path, "resource"))) {
        dir.create(file.path(reference_path, "resource"))
    }
    genome.fa = file.path(reference_path, "resource", "genome.fa")
    if(overwrite || !file.exists(paste0(genome.fa, ".gz"))) {
        message("Saving local copy as FASTA...", appendLF = FALSE)
            if(overwrite && file.exists(paste0(genome.fa, ".gz"))) {
                file.remove(file.exists(paste0(genome.fa, ".gz")))
            }
            rtracklayer::export(genome, genome.fa, "fasta")
        message("Compressing FASTA", appendLF = FALSE)
            R.utils::gzip(genome.fa)
        message("done")
    }
}

.fetch_fasta_save_2bit <- function(genome, reference_path, overwrite) {
    if (!dir.exists(file.path(reference_path, "resource"))) {
        dir.create(file.path(reference_path, "resource"))
    }
    genome.2bit = file.path(reference_path, "resource", "genome.2bit")
    if(overwrite || !file.exists(paste0(genome.2bit, ".2bit"))) {
    # Convert to local 2bit for better memory management
        message("Saving genome as TwoBitFile...", appendLF = FALSE)
            if(overwrite && file.exists(genome.2bit)) {
                file.remove(file.exists(genome.2bit))
            }
            rtracklayer::export(genome, genome.2bit, "2bit")
        message("done\n")
    }
}


################################################################################

.fetch_gtf <- function(reference_path = "./Reference",
        gtf = "", ah_transcriptome = "",  verbose = TRUE,
        convert_chromosome_names = NULL,
        overwrite = FALSE) {
    gtf_gr <- NULL
    if (ah_transcriptome != "") {
        if(substr(ah_transcriptome, 1, 2) != "AH") {
            .log(paste("In .fetch_gtf(),",
                "Given transcriptome AnnotationHub reference is incorrect"))
        }
        gtf_gr <- .fetch_AH(ah_transcriptome, verbose = verbose)
    } else {
        gtf_file = .parse_valid_file(gtf)
        if(!file.exists(gtf_file)) {
            .log(paste("In .fetch_gtf(),",
                "Given transcriptome gtf file", gtf, "not found"))
        }
        message("Reading source GTF file...", appendLF = FALSE)
        gtf_gr <- rtracklayer::import(gtf_file, "gtf")
        message("done\n")
    }
    if(!is.null(convert_chromosome_names)) {
        gtf_gr = .gr_convert_seqnames(gtf_gr, convert_chromosome_names)
    }
    gc()
    # Convert to local gtf.gz for later FeatureCounts
    if (!dir.exists(file.path(reference_path, "resource"))) {
        dir.create(file.path(reference_path, "resource"))
    }
    r_path = file.path(reference_path, "resource")
    if(overwrite || !file.exists(file.path(r_path, "transcripts.gtf.gz"))) {
        message("Saving local copy of GTF file...", appendLF = FALSE)
        if(file.exists(file.path(r_path, "transcripts.gtf"))) {
            file.remove(file.path(r_path, "transcripts.gtf"))
        }
        rtracklayer::export(gtf_gr, file.path(r_path, "transcripts.gtf"), "gtf")
        if(!file.exists(file.path(r_path, "transcripts.gtf"))) {
            .log(paste("In .fetch_gtf(),",
                "Unable to save local copy of gene annotations"))
        }
        message("Compressing GTF file...", appendLF = FALSE)
        if(file.exists(file.path(r_path, "transcripts.gtf.gz"))) {
            file.remove(file.path(r_path, "transcripts.gtf.gz"))
        }
        gzip(filename = file.path(r_path, "transcripts.gtf"),
            destname = file.path(r_path, "transcripts.gtf.gz")
        )
        if(file.exists(file.path(r_path, "transcripts.gtf.gz")) &
                file.exists(file.path(r_path, "transcripts.gtf"))) {
            file.remove(file.path(r_path, "transcripts.gtf"))
        }
    }
    message("done\n")    
    return(gtf_gr)
}

.validate_gtf <- function(genome, gtf_gr) {
    chrOrder <- names(seqinfo(genome))
    if(!any(as.character(GenomicRanges::seqnames(gtf_gr)) %in% chrOrder)) {
        .log(paste("In .validate_gtf(),",
            "Chromosomes in genome and gene annotation does not match",
            "likely incompatible FASTA and GTF file"))
    }
    seqlevels(gtf_gr, pruning.mode = "tidy") <- chrOrder
    return(gtf_gr)
}

.gr_convert_seqnames <- function(gr, convert_chromosome_names) {
    if(!is.null(convert_chromosome_names)) {
        converter = data.table(Original = GenomeInfoDb::seqlevels(gr))
        converter[convert_chromosome_names,
            on = "Original",
            c("NewName") := get("i.New")]
        converter[is.na(get("NewName")),
            c("NewName") := get("Original")]
        chrOrder <- converter$NewName
        GenomeInfoDb::seqlevels(gr) <- chrOrder
    }
    gr
}

################################################################################

.fetch_AH <- function(ah_record_name, rdataclass = c("GRanges", "TwoBitFile"),
        localHub = FALSE, ah = AnnotationHub(localHub = localHub), 
        as_DNAStringSet = TRUE, exclude_scaffolds = TRUE, 
        verbose = FALSE) {
    rdataclass = match.arg(rdataclass)
    if(!substr(ah_record_name, 1, 2) == "AH") {
        .log(paste(ah_record_name,
            "does not appear to be a valid AnnotationHub record name"))
    }
    if(!(ah_record_name %in% names(ah))) {
        .log(paste(ah_record_name,
            "is not found in AnnotationHub index.",
            "Perhaps check online connection or record name"))
    }
    ah_record <- ah[names(ah) == ah_record_name]
    if(ah_record$rdataclass != rdataclass) {
        .log(paste(ah_record_name,
            "is of type", ah_record$rdataclass,
            "and not of expected:", rdataclass))
    }
    if (verbose) {
        message(paste("Downloading", rdataclass,
            "from AnnotationHub, if required..."),
            appendLF = FALSE
        )
    }
    cache_loc <- AnnotationHub::cache(ah_record)
    if (verbose) message("done")
    if(!file.exists(cache_loc)) {
        .log("AnnotationHub cache error - asset not found")
    }   
    if (ah_record$rdataclass == "GRanges") {
        if (verbose) {
            message("Importing to memory as GRanges object...",
                appendLF = FALSE
            )
        }
        gtf <- rtracklayer::import(cache_loc, "gtf")
        if (verbose) message("done\n")
        return(gtf)
    } else if (ah_record$rdataclass == "TwoBitFile") {
        if (verbose) {
            message("Importing to memory as TwoBitFile object...",
                appendLF = FALSE
            )
        }
        twobit <- rtracklayer::TwoBitFile(cache_loc)
        if (verbose) message("done\n")
        if(as_DNAStringSet) {
            message("Importing genome into memory...", appendLF = FALSE)
                genome = rtracklayer::import(twobit)
            message("done")
            if(exclude_scaffolds) {
                genome <- .fetch_fasta_convert_chrom_species(
                    genome, ah_record$species)
            }
            return(genome)
        } else {
            return(twobit)
        }
    }
}

.parse_valid_file <- function(file, msg = "") {
    if (!is_valid(file)) {
        .log(paste("Reference generated without", msg), type = "message")
        return("")
    } else if ( any(startsWith(file, c("http", "ftp")))) {
        url <- file
        # BiocFileCache this and return file path
        cache <- tools::R_user_dir(package = "NxtIRFcore", which="cache")
        bfc <- BiocFileCache::BiocFileCache(cache, ask = FALSE)
        path <- BiocFileCache::bfcrpath(bfc, url)
        return(unname(path))
    } else if (!file.exists(file)) {
        .log(paste(file, "not found.", "Reference generated without", msg),
            type = "message")
        return("")
    } else if (file.exists(file)) {
        return(file)
    } else {
        .log(paste("Reference generated without", msg), type = "message")
        return("")
    }
}

.fix_gtf <- function(gtf_gr) {
    # fix gene / transcript names with '/' (which breaks IRFinder code)
    # 22-08-21: also fix missing gene_name and transcript_name 
    #           in newer Ensembl refs
    if("gene_name" %in% names(S4Vectors::mcols(gtf_gr))) {
        gtf_gr$gene_name[is.na(gtf_gr$gene_name)] = 
            gtf_gr$gene_id[is.na(gtf_gr$gene_name)]
        gtf_gr$gene_name <- gsub("/", "_", gtf_gr$gene_name)
    } else {
        gtf_gr$gene_name <- gtf_gr$gene_id
    }
    
    if("transcript_name" %in% names(S4Vectors::mcols(gtf_gr))) {
        gtf_gr$transcript_name[is.na(gtf_gr$transcript_name)] = 
            gtf_gr$transcript_id[is.na(gtf_gr$transcript_name)]
        gtf_gr$transcript_name <- gsub("/", "_", gtf_gr$transcript_name)
    } else {
        gtf_gr$transcript_name <- gtf_gr$transcript_id
    }
    if(!("gene_biotype" %in% names(S4Vectors::mcols(gtf_gr)))) {
        gtf_gr$gene_biotype = "protein_coding"
    }
    if(!("transcript_biotype" %in% names(S4Vectors::mcols(gtf_gr)))) {
        gtf_gr$transcript_biotype = "protein_coding"
    }
    if(!("transcript_support_level" %in% names(S4Vectors::mcols(gtf_gr)))) {
        gtf_gr$transcript_support_level = 1
    }
    return(gtf_gr)
}

################################################################################
# Sub

.process_gtf <- function(gtf_gr, reference_path) {
    message("Processing gtf file...")    
    message("...genes")
    Genes_group = .process_gtf_genes(gtf_gr, reference_path)
    message("...transcripts")
    .process_gtf_transcripts(gtf_gr, reference_path)
    message("...CDS")
    .process_gtf_misc(gtf_gr, reference_path)
    message("...exons")
    .process_gtf_exons(gtf_gr, reference_path, Genes_group)
    message("...done\n")
}

.process_gtf_genes <- function(gtf_gr, reference_path) {
    Genes <- gtf_gr[gtf_gr$type == "gene"]
    if(length(Genes) == 0) {
        Genes = as.data.table(gtf_gr)
        Genes = Genes[,
        c("start", "end", "width") := list(
            min(get("start")),
            max(get("end")),
            max(get("end")) - min(get("start")) + 1
        ), by = c("seqnames", "strand", "gene_id", "gene_name", "gene_biotype")]
        Genes = unique(Genes, 
            by = c("seqnames", "strand", 
            "gene_id", "gene_name", "gene_biotype"))
        Genes$type = "gene"
        Genes = .grDT(Genes, keep.extra.columns = TRUE)
        if(length(Genes) == 0) {
            .log("No genes detected in reference!")
        }
    }
    Genes <- GenomeInfoDb::sortSeqlevels(Genes)
    Genes <- sort(Genes)
    Genes$gene_display_name <- paste0(Genes$gene_name, " (", Genes$gene_id, ")")

    # Annotate gene_groups_stranded
    Genes_group.stranded <- as.data.table(reduce(Genes))
    setorder(Genes_group.stranded, seqnames, start, strand)

    Genes_group.stranded[, c("gene_group") := .I]
    OL <- findOverlaps(
        Genes, .grDT(Genes_group.stranded)
    )
    Genes$gene_group_stranded[from(OL)] <-
        Genes_group.stranded$gene_group[to(OL)]

    # Annotate gene_groups_unstranded
    Genes_group.unstranded <- as.data.table(reduce(Genes, ignore.strand = TRUE))
    setorder(Genes_group.unstranded, seqnames, start)
    Genes_group.unstranded[, c("gene_group") := .I]
    OL <- findOverlaps(
        Genes,
        .grDT(Genes_group.unstranded, ignore.strand = TRUE)
    )
    Genes$gene_group_unstranded[from(OL)] <-
        Genes_group.unstranded$gene_group[to(OL)]

    write.fst(
        as.data.frame(Genes),
        file.path(reference_path, "fst", "Genes.fst")
    )
    final = list(
        stranded = Genes_group.stranded,
        unstranded = Genes_group.unstranded
    )
    return(final)
}

.process_gtf_transcripts <- function(gtf_gr, reference_path) {
    Transcripts <- gtf_gr[gtf_gr$type == "transcript"]
    if(length(Transcripts) == 0) {
        Transcripts = as.data.table(gtf_gr)
        Transcripts = Transcripts[,
        c("start", "end", "width") := list(
            min(get("start")),
            max(get("end")),
            max(get("end")) - min(get("start")) + 1
        ), by = c("seqnames", "strand", 
            "gene_id", "gene_name", "gene_biotype",
            "transcript_id", "transcript_name", "transcript_biotype")
        ]
        Transcripts = unique(Transcripts, by = c("seqnames", "strand", 
            "gene_id", "gene_name", "gene_biotype",
            "transcript_id", "transcript_name", "transcript_biotype"))
        Transcripts$type = "transcript"
        Transcripts = .grDT(Transcripts, keep.extra.columns = TRUE)
        if(length(Transcripts) == 0) {
            .log("No transcripts detected in reference!")
        }
    }
    Transcripts <- GenomeInfoDb::sortSeqlevels(Transcripts)
    Transcripts <- sort(Transcripts)
    if ("gene_biotype" %in% names(mcols(Transcripts))) {
        # do nothing
    } else if ("gene_type" %in% names(mcols(Transcripts))) {
        colnames(mcols(Transcripts))[which(colnames(mcols(Transcripts)) ==
            "gene_type")] <- "gene_biotype"
    } else {
        mcols(Transcripts)$gene_biotype <- "protein_coding"
    }
    if ("transcript_biotype" %in% names(mcols(Transcripts))) {
        # do nothing
    } else if ("transcript_type" %in% names(mcols(Transcripts))) {
        colnames(mcols(Transcripts))[which(colnames(mcols(Transcripts)) ==
            "transcript_type")] <- "transcript_biotype"
    } else {
        mcols(Transcripts)$transcript_biotype <- "protein_coding"
    }
    if ("transcript_support_level" %in% names(mcols(Transcripts))) {
        Transcripts$transcript_support_level <-
            tstrsplit(Transcripts$transcript_support_level, split = " ")[[1]]
        Transcripts$transcript_support_level[
            is.na(Transcripts$transcript_support_level)
        ] <- "NA"
    }
    write.fst(
        as.data.frame(Transcripts),
        file.path(reference_path, "fst", "Transcripts.fst")
    )
}

.process_gtf_misc <- function(gtf_gr, reference_path) {
    # Proteins
    Proteins <- gtf_gr[gtf_gr$type == "CDS"]
    if(length(Proteins) == 0) {
        .log("No CDS (proteins) detected in reference!")
    }
    Proteins <- GenomeInfoDb::sortSeqlevels(Proteins)
    Proteins <- sort(Proteins)
    write.fst(
        as.data.frame(Proteins),
        file.path(reference_path, "fst", "Proteins.fst")
    )
    # Misc
    gtf.misc <- gtf_gr[!gtf_gr$type %in% c("gene", "transcript", "exon", "CDS")]
    if(length(gtf.misc) == 0) {
        .log("No start / stop codons detected in reference!")
    }
    gtf.misc <- GenomeInfoDb::sortSeqlevels(gtf.misc)
    gtf.misc <- sort(gtf.misc)
    write.fst(
        as.data.frame(gtf.misc),
        file.path(reference_path, "fst", "Misc.fst")
    )
}

.process_gtf_exons <- function(gtf_gr, reference_path, Genes_group) {
    Exons <- gtf_gr[gtf_gr$type == "exon"]
    if(length(Exons) == 0) {
        .log("No exons detected in reference!")
    }
    Exons <- GenomeInfoDb::sortSeqlevels(Exons)
    Exons <- sort(Exons)

    # transcript_biotype is very important field.
    #   If Gencode, this is transcript_type.
    #   In rare case we do not have this field
    #   This next bit ensures transcript_biotype exists.
    if ("transcript_biotype" %in% names(mcols(Exons))) {
    } else if ("transcript_type" %in% names(mcols(Exons))) {
        colnames(mcols(Exons))[
            which(colnames(mcols(Exons)) == "transcript_type")
        ] <-
            "transcript_biotype"
    } else {
        mcols(Exons)$transcript_biotype <- "protein_coding"
    }
    
    # Assign gene groups then bake exon-groups into Exons
    tmp.Exons_group.stranded <- .process_exon_groups(Exons, Genes_group, 
        stranded = TRUE)
    tmp.Exons_group.unstranded <- .process_exon_groups(Exons, Genes_group, 
        stranded = FALSE)

    # Now annotate all exons in Exons with the gene and exon groups
    OL <- findOverlaps(
        Exons, .grDT(tmp.Exons_group.stranded)
    )
    Exons$gene_group_stranded[from(OL)] <-
        tmp.Exons_group.stranded$gene_group[to(OL)]
    Exons$exon_group_stranded[from(OL)] <-
        tmp.Exons_group.stranded$exon_group[to(OL)]

    OL <- findOverlaps(
        Exons, .grDT(tmp.Exons_group.unstranded,
            ignore.strand = TRUE
        )
    )
    Exons$gene_group_unstranded[from(OL)] <-
        tmp.Exons_group.unstranded$gene_group[to(OL)]
    Exons$exon_group_unstranded[from(OL)] <-
        tmp.Exons_group.unstranded$exon_group[to(OL)]

    write.fst(as.data.frame(Exons), 
        file.path(reference_path, "fst", "Exons.fst"))
    write.fst(
        rbind(tmp.Exons_group.stranded, tmp.Exons_group.unstranded),
        file.path(reference_path, "fst", "Exons.Group.fst")
    )
}

.process_exon_groups <- function(Exons, Genes_group, stranded = TRUE) {
    tmp.exons.exclude <- Exons[!grepl("intron", Exons$transcript_biotype)]
    tmp.Exons_group <- as.data.table(reduce(tmp.exons.exclude,
        ignore.strand = !stranded
    ))
    if(stranded) {
        GG = Genes_group[["stranded"]]
    } else {
        GG = Genes_group[["unstranded"]]
    }
    OL <- findOverlaps(
        .grDT(tmp.Exons_group),
        .grDT(GG, ignore.strand = !stranded)
    )
    tmp.Exons_group$gene_group[from(OL)] <-
        GG$gene_group[to(OL)]
    # Some retained_intron transcripts have terminal exons lying outside of
    #   main transcripts. Include these also
    tmp.exons.exclude.span <- split(
        .grDT(tmp.Exons_group),
        tmp.Exons_group$gene_group
    )
    tmp.exons.exclude.span <-
        unlist(range(tmp.exons.exclude.span), use.names = TRUE)
    tmp.exons.RI <- Exons[grepl("intron", Exons$transcript_biotype)]
    if (length(tmp.exons.RI) > 0) {
        OL <- findOverlaps(
            tmp.exons.RI,
            tmp.exons.exclude.span,
            ignore.strand = !stranded
        )
        tmp.exons.RI <- tmp.exons.RI[-from(OL)]
        tmp.exons.exclude <- c(tmp.exons.exclude, tmp.exons.RI)

        tmp.Exons_group <- as.data.table(reduce(tmp.exons.exclude,
            ignore.strand = !stranded
        ))
        OL <- findOverlaps(
            .grDT(tmp.Exons_group),
            .grDT(GG, ignore.strand = !stranded)
        )
        tmp.Exons_group$gene_group[from(OL)] <-
            GG$gene_group[to(OL)]
    }
    setorder(tmp.Exons_group, seqnames, start, strand)
    tmp.Exons_group[, c("exon_group") :=
        data.table::rowid(get("gene_group"))]
    if(stranded) {
        tmp.Exons_group[get("strand") == "-",
            c("exon_group") := max(get("exon_group")) + 1 - get("exon_group"),
            by = "gene_group"
        ]
    }
    return(tmp.Exons_group)
}

################################################################################
# Sub

.process_introns <- function(reference_path, genome, 
        UseExtendedTranscripts = TRUE) {
    message("Processing introns...")

    message("...data")
    data = .process_introns_data(reference_path, genome, UseExtendedTranscripts)
    gc()
    data[["candidate.introns"]] = .process_introns_annotate(
        data[["candidate.introns"]], data[["Transcripts"]], genome,
        data[["Proteins"]], data[["Exons"]]
    )
    message("...defining flanking exon islands")  
    data[["candidate.introns"]] = .process_introns_group(
        data[["candidate.introns"]], data[["Exons_group.stranded"]],
        data[["Exons_group.unstranded"]]
    )
    gc()
    write.fst(
        data[["candidate.introns"]],
        file.path(reference_path, "fst", "junctions.fst")
    )
    message("done\n")
}

.process_introns_data <- function(reference_path, genome, 
        UseExtendedTranscripts = TRUE) {
    Exons <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Exons.fst")),
    )
    Transcripts <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Transcripts.fst")),
    )
    Proteins <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Proteins.fst")),
    )
    Exons_group <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Exons.Group.fst")),
    )
    Exons_group.stranded <- Exons_group[get("strand") != "*"]
    Exons_group.unstranded <- Exons_group[get("strand") == "*"]

    if (UseExtendedTranscripts == FALSE) {
        candidate.transcripts <- Exons[get("transcript_biotype") %in%
            c("processed_transcript", "protein_coding")]
    } else {
        candidate.transcripts <- Exons
    }
    candidate.introns <- as.data.table(
        .grlGaps(
            split(.grDT(candidate.transcripts), 
            candidate.transcripts$transcript_id)
        )
    )
    candidate.introns[, c("group") := NULL]
    setnames(candidate.introns, "group_name", "transcript_id")
    setorderv(candidate.introns, c("seqnames", "start", "end", "strand"))
    
    final = list(
        Exons = Exons, Transcripts = Transcripts, Proteins = Proteins,
        Exons_group.stranded = Exons_group.stranded,
        Exons_group.unstranded = Exons_group.unstranded,
        candidate.introns = candidate.introns
    )
    return(final)
}

#############################################################################

.process_introns_annotate <- function(candidate.introns, Transcripts, genome,
        Proteins, Exons) {
    # Annotating Introns:
    message("...basic annotations")  
    candidate.introns <- .process_introns_annotate_basics(
        candidate.introns, Transcripts)

    # Grab splice motifs at this point; filter by valid splice motifs
    message("...splice motifs")  
    candidate.introns <- .process_introns_annotate_splice_motifs(
        candidate.introns, genome)

    # Do other annotations here:
    message("...other motifs")  
    candidate.introns <- .process_introns_annotate_others(
        candidate.introns, Transcripts, Proteins, Exons
    )
    gc()
    return(candidate.introns)
}

.process_introns_annotate_basics <- function(
        candidate.introns, Transcripts) {
    candidate.introns[, c("intron_number") :=
        data.table::rowid(get("transcript_id"))]
    candidate.introns[get("strand") == "-",
        c("intron_number") :=
            max(get("intron_number")) + 1 - get("intron_number"),
        by = "transcript_id"
    ]
    
    candidate.introns[, c("intron_id") :=
        paste0(get("transcript_id"), "_Intron", get("intron_number"))]
    
    candidate.introns[Transcripts,
        on = "transcript_id",
        c("gene_name", "gene_id", "transcript_name", "transcript_biotype") :=
            list(
                get("i.gene_name"), get("i.gene_id"),
                get("i.transcript_name"), get("i.transcript_biotype")
            )
    ]
    return(candidate.introns)
}

.process_introns_annotate_splice_motifs <- function(
        candidate.introns, genome) {

    donor.introns <- data.frame(
        seqnames = candidate.introns$seqnames,
        start = ifelse(candidate.introns$strand == "+",
            candidate.introns$start, candidate.introns$end - 1
        ),
        stop = ifelse(candidate.introns$strand == "+",
            candidate.introns$start + 1, candidate.introns$end
        ),
        strand = candidate.introns$strand
    )
    donor_seq <- getSeq(genome, makeGRangesFromDataFrame(donor.introns))
    acceptor.introns <- data.frame(
        seqnames = candidate.introns$seqnames,
        start = ifelse(candidate.introns$strand == "+",
            candidate.introns$end - 1, candidate.introns$start
        ),
        stop = ifelse(candidate.introns$strand == "+",
            candidate.introns$end, candidate.introns$start + 1
        ),
        strand = candidate.introns$strand
    )
    acceptor_seq <- getSeq(genome, makeGRangesFromDataFrame(acceptor.introns))
    candidate.introns$splice_motif <- paste0(donor_seq, acceptor_seq)
    
    return(candidate.introns)
}

.process_introns_annotate_others <- function(
        candidate.introns, Transcripts, Proteins, Exons) {
    candidate.introns[Transcripts,
        on = "transcript_id",
        c("gene_name", "gene_id", "transcript_name") :=
            list(get("i.gene_name"), get("i.gene_id"), get("i.transcript_name"))
    ]
    if ("transcript_support_level" %in% colnames(Transcripts)) {
        candidate.introns[Transcripts,
            on = "transcript_id",
            c("transcript_support_level") :=
                list(get("i.transcript_support_level"))
        ]
        candidate.introns[, c("transcript_support_level") :=
            tstrsplit(get("transcript_support_level"), split = " ")[[1]]]
        candidate.introns[
            is.na(get("transcript_support_level")),
            c("transcript_support_level") := "NA"
        ]
    }
    if ("protein_id" %in% colnames(Proteins)) {
        Proteins.red = unique(Proteins[, c("transcript_id", "protein_id")])
        candidate.introns[Proteins.red,
            on = "transcript_id",
            c("protein_id") := list(get("i.protein_id"))
        ]
    }
    if ("ccds_id" %in% colnames(Exons)) {
        Exons.red = unique(Exons[, c("transcript_id", "ccds_id")])        
        candidate.introns[Exons.red,
            on = "transcript_id",
            c("ccds_id") := list(get("i.ccds_id"))
        ]
    }
    return(candidate.introns)
}

#############################################################################
.process_introns_group <- function(candidate.introns, 
        Exons_group.stranded, Exons_group.unstranded) {
    # Cannot annotate candidate introns by min and max exon_groups
    # because retained introns will overlap one or more exon groups
    # need to walk start -=1, end += 1, then do the overlap thing
    candidate.introns[, c("intron_start") := get("start")]
    candidate.introns[, c("intron_end") := get("end")]

    candidate.introns <- .process_introns_group_overlap(
        candidate.introns, Exons_group.stranded,
        c("gene_group_stranded", "exon_group_stranded_upstream",
            "gene_group_stranded", "exon_group_stranded_downstream"),
        c("gene_group", "exon_group", "gene_group", "exon_group")
    )
    # Need fix for retained_introns or sense_intronic where junction extends
    #   into the obligate introns
    candidate.introns <- .process_introns_group_fix_RI(
        candidate.introns, Exons_group.stranded, 
        c("gene_group_stranded", "exon_group_stranded_upstream",
            "gene_group_stranded", "exon_group_stranded_downstream"),
        c("gene_group", "intron_number", "gene_group", "intron_number")
    )
    
    # Now repeat the same for unstranded condition
    candidate.introns <- .process_introns_group_overlap(
        candidate.introns, Exons_group.unstranded,
        c("gene_group_unstranded", "exon_group_unstranded_upstream",
            "gene_group_unstranded", "exon_group_unstranded_downstream"),
        c("gene_group", "exon_group", "gene_group", "exon_group")
    )

    # Need fix for retained_introns or sense_intronic where
    #   junction extends into the obligate introns
    candidate.introns <- .process_introns_group_fix_RI(
        candidate.introns, Exons_group.unstranded, 
        c("gene_group_unstranded", "exon_group_unstranded_upstream",
            "gene_group_unstranded", "exon_group_unstranded_downstream"),
        c("gene_group", "intron_number", "gene_group", "intron_number")
    )

    # reset
    candidate.introns[, c("start") := get("intron_start")]
    candidate.introns[, c("end") := get("intron_end")]
    candidate.introns[, c("Event") := paste0(
        get("seqnames"), ":", get("intron_start"), "-",
        get("intron_end"), "/", get("strand")
    )]
    return(candidate.introns)
}

.process_introns_group_overlap <- function(target.DT, groups.DT,
    target.columns, groups.columns) {

    OL <- .overlaps_exon_island(
        target.DT,
        groups.DT,
        upstream = TRUE
    )
    set(target.DT, from(OL), 
        target.columns[1],
        groups.DT[, get(groups.columns[1])][to(OL)]
    )
    set(target.DT, from(OL), 
        target.columns[2],
        groups.DT[, get(groups.columns[2])][to(OL)]
    )
    OL <- .overlaps_exon_island(
        target.DT,
        groups.DT,
        upstream = FALSE
    )
    set(target.DT, from(OL), 
        target.columns[3],
        groups.DT[, get(groups.columns[3])][to(OL)]
    )
    set(target.DT, from(OL), 
        target.columns[4],
        groups.DT[, get(groups.columns[4])][to(OL)]
    )
    return(target.DT)
}
    
.process_introns_group_fix_RI <- function(
        target.DT, groups.DT, 
        target.columns, groups.columns) {
    tmp <- .grDT(groups.DT, keep.extra.columns = TRUE)
    tmp.Introns_group <- .grlGaps(
        split(tmp, tmp$gene_group)
    )
    tmp.Introns_group <- as.data.table(tmp.Introns_group)
    setnames(tmp.Introns_group, "group_name", "gene_group")
    tmp.Introns_group[, c("intron_number") :=
        data.table::rowid(get("gene_group"))]
    tmp.Introns_group[get("strand") == "-",
        c("intron_number") :=
            max(get("intron_number")) + 1 - get("intron_number"),
        by = "gene_group"
    ]

    target.DT.subset <- target.DT[is.na(get(target.columns[2]))]
    target.DT <- target.DT[!is.na(get(target.columns[2]))]
    OL <- .overlaps_exon_island(
        target.DT.subset,
        tmp.Introns_group,
        upstream = TRUE
    )
    set(target.DT.subset, from(OL), 
        target.columns[1],
        as.integer(tmp.Introns_group[, get(groups.columns[1])][to(OL)])
    )
    set(target.DT.subset, from(OL), 
        target.columns[2],
        tmp.Introns_group[, get(groups.columns[2])][to(OL)]
    )
    target.DT <- rbind(target.DT, target.DT.subset)

    target.DT.subset <- target.DT[is.na(get(target.columns[4]))]
    target.DT <- target.DT[!is.na(get(target.columns[4]))]
    OL <- .overlaps_exon_island(
        target.DT.subset,
        tmp.Introns_group,
        upstream = FALSE
    )
    set(target.DT.subset, from(OL), 
        target.columns[3],
        as.integer(tmp.Introns_group[, get(groups.columns[3])][to(OL)])
    )
    set(target.DT.subset, from(OL), 
        target.columns[4],
        as.integer(tmp.Introns_group[, get(groups.columns[4])][to(OL)] + 1)
    )
    return(rbind(target.DT, target.DT.subset))
}

.overlaps_exon_island <- function(intron.DT, groups.DT, upstream = TRUE) {
    if(all(c("intron_start", "intron_end") %in% colnames(intron.DT))) {
        int.DT = intron.DT[, c("seqnames", "start", "end", "strand", 
            "intron_start", "intron_end")]
    } else {
        int.DT = intron.DT[, c("seqnames", "start", "end", "strand")]
        int.DT[, c("intron_start", "intron_end") := 
            list(get("start"), get("end"))]
    }
    if (upstream) {
        int.DT[get("strand") == "+", c("start") := get("intron_start") - 1]
        int.DT[get("strand") == "+", c("end") := get("intron_start")]
        int.DT[get("strand") == "-", c("start") := get("intron_end")]
        int.DT[get("strand") == "-", c("end") := get("intron_end") + 1]
    } else {
        int.DT[get("strand") == "-", c("start") := get("intron_start") - 1]
        int.DT[get("strand") == "-", c("end") := get("intron_start")]
        int.DT[get("strand") == "+", c("start") := get("intron_end")]
        int.DT[get("strand") == "+", c("end") := get("intron_end") + 1]
    }
    OL <- findOverlaps(
        .grDT(int.DT),
        .grDT(groups.DT)
    )
    return(OL)
}

################################################################################
# Sub

.gen_irf <- function(reference_path, extra_files, genome,
        convert_chromosome_names) {
    message("Generating IRFinder reference: ref-cover.bed...")

    # Generating IRFinder-base references
    message("...prepping data")  
    data = .gen_irf_prep_data(reference_path)
    data2 = .gen_irf_prep_introns(
        data[["candidate.introns"]], data[["Exons"]], 
        extra_files, convert_chromosome_names
    )
    data2[["introns.unique"]] <- .gen_irf_prep_introns_unique(
        data2[["introns.unique"]], data2[["exclude.directional"]],
        data[["Genes.rev"]], data[["Genes.Extended"]]
    )
    message("...determining measurable introns (directional)")  
    tmpdir.IntronCover.summa = .gen_irf_export_introncover(
        .gen_irf_exclusion_zones(
            data2[["introns.unique"]], data2[["exclude.omnidirectional"]],
            data2[["exclude.directional"]], stranded = TRUE
        ), stranded = TRUE, reference_path, data2[["introns.unique"]]
    )
    message("...determining measurable introns (non-directional)")  
    tmpnd.IntronCover.summa = .gen_irf_export_introncover(
        tmpnd.IntronCover <- .gen_irf_exclusion_zones(
            data2[["introns.unique"]], data2[["exclude.omnidirectional"]],
            data2[["exclude.directional"]], stranded = FALSE
        ), stranded = FALSE, reference_path, data2[["introns.unique"]]
    )
    message("...writing final BED")  
    # Generate final ref-cover.bed
    ref.cover = .gen_irf_refcover(reference_path)
    message("done\n")
    message("Generating IRFinder reference: ref-ROI.bed ...", appendLF = FALSE)
    ref.ROI <- .gen_irf_ROI(reference_path, extra_files,
        convert_chromosome_names, genome, 
        data[["Genes"]], data[["Transcripts"]])
    message("done\n")
    message("Generating IRFinder reference: ref-read-continues.ref ...", 
        appendLF = FALSE)
    readcons = .gen_irf_readcons(reference_path,
        tmpdir.IntronCover.summa, tmpnd.IntronCover.summa
    )
    message("done\n")
    message("Generating IRFinder reference: ref-sj.ref ...", appendLF = FALSE)
    ref.sj <- .gen_irf_sj(reference_path)
    message("done\n")
    .gen_irf_final(reference_path,ref.cover, readcons, ref.ROI, ref.sj)
}
################################################################################

.gen_irf_prep_data <- function(reference_path) {
    Genes <- .grDT(
        read.fst(file.path(reference_path, "fst", "Genes.fst")),
        keep.extra.columns = TRUE
    )
    Genes.rev <- Genes
    strand(Genes.rev) <- ifelse(strand(Genes.rev) == "+", "-",
        ifelse(strand(Genes.rev) == "-", "+", "*")
    ) # Invert strand
    Genes.Extended <- reduce(c(
        flank(Genes.rev, 5000),
        flank(Genes.rev, 1000, start = FALSE)
    ))

    candidate.introns <- as.data.table(
        read.fst(file.path(reference_path, "fst", "junctions.fst"))
    )
    Exons <- .grDT(
        read.fst(file.path(reference_path, "fst", "Exons.fst")),
        keep.extra.columns = TRUE
    )
    Transcripts <- .grDT(
        read.fst(file.path(reference_path, "fst", "Transcripts.fst")),
        keep.extra.columns = TRUE
    )
    candidate.introns <- candidate.introns[get("transcript_biotype") %in%
        c("protein_coding", "processed_transcript",
            "lincRNA", "antisense", "nonsense_mediated_decay")]
    candidate.introns[, c("transcript_biotype") :=
        factor(get("transcript_biotype"),
            c("protein_coding", "processed_transcript",
                "lincRNA", "antisense", "nonsense_mediated_decay"),
            ordered = TRUE
        )]
    if ("transcript_support_level" %in% colnames(candidate.introns)) {
        # Sort by tsl first, then reverse later
        setorderv(
            candidate.introns,
            c("transcript_biotype", "transcript_support_level")
        )
    } else {
        setorderv(candidate.introns, "transcript_biotype")
    }
    setorderv(candidate.introns, c("seqnames", "start", "end", "strand"))
    final = list(
        Genes = Genes, Genes.rev = Genes.rev,
        Genes.Extended = Genes.Extended,
        Exons = Exons, Transcripts = Transcripts,
        candidate.introns = candidate.introns
    )
    return(final)
}

.gen_irf_prep_introns <- function(candidate.introns, Exons,
        extra_files, convert_chromosome_names) {
    tmp.exons.exclude <- Exons[!grepl("intron", Exons$transcript_biotype)]
    introns.unique <- unique(candidate.introns,
        by = c("seqnames", "start", "end", "width", "strand"))
    setorderv(introns.unique, c("seqnames", "start", "end", "strand"))
    introns.unique <- .grDT(introns.unique, keep.extra.columns = TRUE)

    exclude.directional <- as.data.table(tmp.exons.exclude)
    exclude.directional <- unique(exclude.directional,
        by = c("seqnames", "start", "end", "width", "strand"))
    exclude.directional[, c("start") := get("start") - 5]
    exclude.directional[, c("end") := get("end") + 5]

    exclude.omnidirectional <- GRanges(NULL)
    if (extra_files$MappabilityFile != "") {
        exclude.omnidirectional <- c(exclude.omnidirectional,
            rtracklayer::import(extra_files$MappabilityFile, "bed"))
    }
    if (extra_files$BlacklistFile != "") {
        exclude.omnidirectional <- c(exclude.omnidirectional,
            rtracklayer::import(extra_files$BlacklistFile, "bed"))
    }
    if(!is.null(convert_chromosome_names)) {
        exclude.omnidirectional = .gr_convert_seqnames(
            exclude.omnidirectional,convert_chromosome_names)
    }
    # merge with any gaps <= 9
    exclude.omnidirectional <-
        reduce(exclude.omnidirectional, min.gapwidth = 9)
    # clean introns by those lying completely within blacklist regions
    if (length(exclude.omnidirectional) > 0) {
        introns.unique.blacklisted <- findOverlaps(introns.unique,
            exclude.omnidirectional,
            type = "within"
        )
        introns.unique <- introns.unique[-introns.unique.blacklisted@from]
    }
    
    final = list(
        introns.unique = introns.unique,
        exclude.directional = exclude.directional,
        exclude.omnidirectional = exclude.omnidirectional
    )
    return(final)
}

.gen_irf_prep_introns_unique <- function(introns.unique, exclude.directional,
        Genes.rev, Genes.Extended) {
    introns.unique.exon.dir <- findOverlaps(introns.unique,
        .grDT(exclude.directional),
        type = "within"
    )
    introns.unique.exon.nd <- findOverlaps(introns.unique,
        .grDT(exclude.directional),
        type = "within", ignore.strand = TRUE
    )

    introns.unique$known_exon_dir <-
        (seq_len(length(introns.unique)) %in% introns.unique.exon.dir@from)
    introns.unique$known_exon_nd <-
        (seq_len(length(introns.unique)) %in% introns.unique.exon.nd@from)

    introns.unique.antiover <- findOverlaps(introns.unique, Genes.rev)
    introns.unique.antinear <- findOverlaps(introns.unique, Genes.Extended)

    introns.unique$antiover <-
        (seq_len(length(introns.unique)) %in% introns.unique.antiover@from)
    introns.unique$antinear <-
        (seq_len(length(introns.unique)) %in% introns.unique.antinear@from)

    # Now subset introns by punching holes using blacklist regions
    introns.unique$intron_width <- BiocGenerics::width(introns.unique)

    # Remove introns less than 50 bp:
    introns.unique <- introns.unique[BiocGenerics::width(introns.unique) > 50]

    # remove 5 bases from start & end
    BiocGenerics::start(introns.unique) <-
        BiocGenerics::start(introns.unique) + 5
    BiocGenerics::end(introns.unique) <-
        BiocGenerics::end(introns.unique) - 5

    return(introns.unique)
}

.gen_irf_exclusion_zones <- function(introns.unique,
        exclude.omnidirectional, exclude.directional,
        stranded = TRUE) {
    if(!stranded) {
        exclude.directional.reverse <- copy(exclude.directional)
        exclude.directional.reverse[get("strand") == "-", c("strand") := "P"]
        exclude.directional.reverse[get("strand") == "+", c("strand") := "-"]
        exclude.directional.reverse[get("strand") == "P", c("strand") := "+"]
        exclude.directional.gr = c(.grDT(exclude.directional), 
            .grDT(exclude.directional.reverse))
    } else {
        exclude.directional.gr = .grDT(exclude.directional)
    }
    if(length(exclude.omnidirectional) > 0) {
        introns.intersect <- GenomicRanges::intersect(introns.unique,
            c(exclude.omnidirectional, exclude.directional.gr))
    } else if(length(exclude.directional) > 0) {
        introns.intersect <- GenomicRanges::intersect(introns.unique, 
            exclude.directional.gr)
    }
    OL <- findOverlaps(introns.unique, introns.intersect)
    # make a GRanges same size as the number of intersections
    introns.intersect.final <- introns.intersect[to(OL)]
    introns.intersect.final$intron_id <- introns.unique$intron_id[from(OL)]

    introns.unique.ID <- split(introns.unique, introns.unique$intron_id)
    introns.intersect.ID <- split(introns.intersect.final,
            introns.intersect.final$intron_id)
    introns.unique.ID.compare <- introns.unique.ID[
        names(introns.unique.ID) %in% names(introns.intersect.ID)
    ]

    IntronCover <- setdiff(introns.unique.ID.compare, introns.intersect.ID)
    # now add back introns that did not require intersection
    #   (or would have been excluded as known-exons)
    IntronCover <- c(IntronCover,
        introns.unique.ID[!(names(introns.unique.ID) %in% 
            names(introns.intersect.ID))])
    if(stranded) {
        IntronCover = c(IntronCover,
            introns.unique.ID[names(introns.unique.ID) %in%
                introns.unique$intron_id[introns.unique$known_exon_dir == TRUE]]
        )
    } else {
        IntronCover = c(IntronCover,
            introns.unique.ID[names(introns.unique.ID) %in%
                introns.unique$intron_id[introns.unique$known_exon_nd == TRUE]]
        )
    }
    IntronCover = as.data.table(IntronCover)
    IntronCover <- IntronCover[,
        c("seqnames", "start", "end", "strand", "width", "group_name")
    ]
    colnames(IntronCover)[6] <- "intron_id"    
    return(IntronCover)
}

.gen_irf_export_introncover <- function(IntronCover, stranded = TRUE,
        reference_path, introns.unique) {
    IntronCover.summa <- IntronCover
    IntronCover.summa[,
        c("num_blocks", "inclbases") := list(.N, sum(get("width"))),
        by = "intron_id"]
    IntronCover.summa <- unique(
        IntronCover.summa[, c("intron_id", "num_blocks", "inclbases")],
        by = "intron_id")
    IntronCover.summa[as.data.table(introns.unique),
        on = "intron_id",
        c("seqnames", "intron_start", "intron_end",
            "intron_width", "width", "strand", "gene_name", "transcript_id")
        := list(get("i.seqnames"), get("i.intron_start"),
                get("i.intron_end"), get("i.intron_width"),
                get("i.width"), get("i.strand"), get("i.gene_name"),
                get("i.transcript_id"))
    ]
    if(stranded) {
        IntronCover.summa[as.data.table(introns.unique), on = "intron_id", 
            c("known_exon_dir", "GG", "EG_up", "EG_down") := 
            list(get("i.known_exon_dir"), get("i.gene_group_stranded"), 
                get("i.exon_group_stranded_upstream"),
                get("i.exon_group_stranded_downstream"))
        ]
    } else {
        IntronCover.summa[as.data.table(introns.unique), on = "intron_id", 
            c("known_exon_nd", "antiover", "antinear", 
                "GG", "EG_up", "EG_down") := 
            list(get("i.known_exon_nd"), get("i.antiover"), get("i.antinear"),
                get("i.gene_group_unstranded"),
                get("i.exon_group_unstranded_upstream"),
                get("i.exon_group_unstranded_downstream"))
        ]    
    }
    IntronCover.summa[, c("exclbases") :=
        get("intron_width") - get("inclbases")]
    # Exclude exclbases / width > 0.3
    IntronCover.summa <-
        IntronCover.summa[get("exclbases") / get("intron_width") < 0.3]
    IntronCover <- semi_join.DT(IntronCover, IntronCover.summa,
        by = "intron_id")
    IntronCover.summa <- .gen_irf_irfname(IntronCover.summa, 
        stranded = stranded)
    
    IntronCover <- .grDT(IntronCover, 
        keep.extra.columns = TRUE)
    IntronCover <- split(IntronCover, IntronCover$intron_id)
    names(IntronCover) <- IntronCover.summa$IRFname[match(
        names(IntronCover), IntronCover.summa$intron_id)]
    setorderv(IntronCover.summa, 
        c("seqnames", "intron_start", "intron_end", "strand"))
    IntronCover <- IntronCover[IntronCover.summa$IRFname]

    rtracklayer::export(IntronCover, file.path(reference_path,
        ifelse(stranded, "tmpdir.IntronCover.bed", "tmpnd.IntronCover.bed")
    ))
    write.fst(IntronCover.summa, file.path(
        reference_path, "fst",
        ifelse(stranded, "Introns.Dir.fst", "Introns.ND.fst")
    ))
    return(IntronCover.summa)
}

.gen_irf_irfname <- function(IntronCover.summa, stranded = TRUE) {
    if(stranded) {
        IntronCover.summa[, c("IRFname") := paste("dir", 
            get("gene_name"), get("intron_id"), get("strand"), 
            get("num_blocks"), sprintf("%.f", get("intron_start") - 1),
            sprintf("%.f", get("intron_end")), get("inclbases"), 
            get("exclbases"),
            ifelse(get("known_exon_dir"), "known-exon", "clean"), sep = "/"
        )]
    } else {
        IntronCover.summa[, c("IRFname") := paste("nd",
            get("gene_name"), get("intron_id"), get("strand"), 
            get("num_blocks"), sprintf("%.f", get("intron_start") - 1),
            sprintf("%.f", get("intron_end")), get("inclbases"), 
            get("exclbases"), sep = "/"
        )]
        # casewise naming of last condition
        IntronCover.summa[
            get("known_exon_nd") & get("antiover") & get("antinear"),
            c("IRFname") := paste(get("IRFname"),
                "known-exon+anti-over+anti-near", sep = "/")]
        IntronCover.summa[
            get("known_exon_nd") & get("antiover") & !get("antinear"),
            c("IRFname") := paste(get("IRFname"), 
                "known-exon+anti-over", sep = "/")]
        IntronCover.summa[
            get("known_exon_nd") & !get("antiover") & get("antinear"),
            c("IRFname") := paste(get("IRFname"),
                "known-exon+anti-near", sep = "/")]
        IntronCover.summa[
            !get("known_exon_nd") & get("antiover") & get("antinear"),
            c("IRFname") := paste(get("IRFname"), 
                "anti-over+anti-near",sep = "/")]
        IntronCover.summa[
            !get("known_exon_nd") & !get("antiover") & get("antinear"),
            c("IRFname") := paste(get("IRFname"), "anti-near", sep = "/")]
        IntronCover.summa[
            !get("known_exon_nd") & get("antiover") & !get("antinear"),
            c("IRFname") := paste(get("IRFname"), "anti-over", sep = "/")]
        IntronCover.summa[
            get("known_exon_nd") & !get("antiover") & !get("antinear"),
            c("IRFname") := paste(get("IRFname"), "known-exon", sep = "/")]
        IntronCover.summa[
            !get("known_exon_nd") & !get("antiover") & !get("antinear"),
            c("IRFname") := paste(get("IRFname"),"clean", sep = "/")]
    }
    return(IntronCover.summa)
}

.gen_irf_refcover <- function(reference_path) {
    tmpdir.IntronCover <- fread(file.path(
        reference_path, "tmpdir.IntronCover.bed"
    ), sep = "\t")
    tmpdir.IntronCover[, c("cat") := "dir"]
    tmpnd.IntronCover <- fread(file.path(
        reference_path, "tmpnd.IntronCover.bed"
    ), sep = "\t")
    tmpnd.IntronCover[, c("cat") := "nd"]

    ref.cover <- rbind(tmpdir.IntronCover, tmpnd.IntronCover)
    setorderv(ref.cover, c("V1", "V2", "V3", "V6", "cat"))
    ref.cover$cat <- NULL
    ref.cover[, c("V9") := as.character(get("V9"))]
    ref.cover[, c("V9") := "255,0,0"]

    return(ref.cover)
}

.gen_irf_ROI <- function(reference_path, extra_files, 
        convert_chromosome_names, genome,
        Genes, Transcripts) {
    rRNA <- as.data.table(Transcripts[grepl("rRNA", Transcripts$gene_biotype)])
    if(nrow(rRNA) > 0) {
        rRNA[, c("start") := get("start") - 1]
        rRNA[, c("name") := paste("rRNA", get("seqnames"), get("start"), 
            get("end"), get("strand"), get("transcript_id"), 
            get("gene_biotype"), get("gene_id"), get("gene_name"), sep = "/")]
        rRNA <- rRNA[, c("seqnames", "start", "end", "name"), with = FALSE]
    } else {
        rRNA <- c()
    }

    if (extra_files$nonPolyAFile != "") {
        nonPolyA <- rtracklayer::import(extra_files$nonPolyAFile, "bed")
        if(!is.null(convert_chromosome_names)) {
            nonPolyA = .gr_convert_seqnames(nonPolyA, convert_chromosome_names)
        }
        nonPolyA = as.data.table(nonPolyA)
        nonPolyA <- nonPolyA[, c("seqnames", "start", "end"), with = FALSE]
        nonPolyA[, c("name") := "NonPolyA"]
    } else {
        nonPolyA <- c()
    }

    AllChr <- makeGRangesListFromDataFrame(data.frame(
        seqnames = names(seqinfo(genome)),
        start = 1, end = seqlengths(seqinfo(genome)),
        names = names(seqinfo(genome))
    ), split.field = "names")
    Genes.chr <- c(
        Genes, flank(Genes, 10000), flank(Genes, 10000, start = FALSE)
    )
    Genes.chr <- reduce(Genes.chr, min.gapwidth = 1000)
    Genes.chr$chr <- seqnames(Genes.chr)
    Genes.chr <- split(Genes.chr, Genes.chr$chr)
    AllChr <- AllChr[names(Genes.chr)]
    AllChr.split <- setdiff(AllChr, Genes.chr, ignore.strand = TRUE)
    Intergenic <- unlist(AllChr.split)
    if (length(Intergenic) > 0) {
        names(Intergenic) <- seq_len(length(Intergenic))
        Intergenic = as.data.table(Intergenic)
        Intergenic <- Intergenic[, c("seqnames", "start", "end"), with = FALSE]
        Intergenic[, c("name") := paste("Intergenic", Intergenic$seqnames, 
            sep = "/")]
    } else {
        Intergenic <- c()
    }
    ref.ROI <- rbind(rRNA, nonPolyA, Intergenic) 
    if(!is.null(ref.ROI) && nrow(ref.ROI) > 0) {
        # ref.ROI = as.data.table(ref.ROI)
        setorderv(ref.ROI, c("seqnames", "start"))
        ref.ROI[, c("start") := get("start") - 1] # convert back to 0-based
    }
    return(ref.ROI)
}

.gen_irf_readcons <- function(reference_path,
        tmpdir.IntronCover.summa, tmpnd.IntronCover.summa) {
    # ref-read-continues.ref
    introns.unique.readcons <- rbind(
        tmpdir.IntronCover.summa[
            ,
            c("seqnames", "intron_start", "intron_end", "strand")
        ],
        tmpnd.IntronCover.summa[
            ,
            c("seqnames", "intron_start", "intron_end", "strand")
        ]
    )
    # 0-based
    introns.unique.readcons[, c("intron_start") := get("intron_start") - 1]

    readcons.left <- introns.unique.readcons[
        ,
        c("seqnames", "intron_start", "strand")
    ]
    readcons.right <- introns.unique.readcons[
        ,
        c("seqnames", "intron_end", "strand")
    ]
    colnames(readcons.left) <- c("V1", "V2", "V3")
    colnames(readcons.right) <- c("V1", "V2", "V3")
    readcons <- rbind(readcons.left, readcons.right)
    setorderv(readcons, c("V1", "V2", "V3"))
    readcons <- unique(readcons)
    gc()
    return(readcons)    
}
.gen_irf_sj <- function(reference_path) {
    # ref-sj.ref
    # Reload candidate introns here, as we've filtered this before
    candidate.introns <- as.data.table(
        read.fst(file.path(reference_path, "fst", "junctions.fst"))
    )

    ref.sj <- candidate.introns[, c("seqnames", "start", "end", "strand")]
    # annotate NMD-unique junctions

    ref.sj$Is_NMD <- ifelse(
        grepl(
            "nonsense_mediated_decay",
            candidate.introns$transcript_biotype
        ),
        "NMD", ""
    )

    ref.sj <- ref.sj[, lapply(.SD, function(x) {
        ifelse(all(x != ""), "NMD", "")
    }), by = c("seqnames", "start", "end", "strand")]

    ref.sj[, c("start") := get("start") - 1]
    setorderv(ref.sj, c("seqnames", "start", "end", "strand"))
    gc()
    return(ref.sj)
}
.gen_irf_final <- function(reference_path,
        ref.cover, readcons, ref.ROI, ref.sj) {
    IRF_file <- file.path(reference_path, "IRFinder.ref")
    # Concatenate all 4 reference files into one file
    fwrite(list("# ref-cover.bed"), IRF_file,
        sep = "\t", eol = "\n", col.names = FALSE, scipen = 50
    )
    fwrite(ref.cover, IRF_file,
        append = TRUE,
        sep = "\t", eol = "\n", col.names = FALSE, scipen = 50
    )
    fwrite(list("# ref-read-continues.ref"), IRF_file,
        append = TRUE,
        sep = "\t", eol = "\n", col.names = FALSE, scipen = 50
    )
    fwrite(readcons, IRF_file,
        append = TRUE,
        sep = "\t", eol = "\n", col.names = FALSE, scipen = 50
    )
    fwrite(list("# ref-ROI.bed"), IRF_file,
        append = TRUE,
        sep = "\t", eol = "\n", col.names = FALSE, scipen = 50
    )
    if(!is.null(ref.ROI) && nrow(ref.ROI) > 0) {
        fwrite(ref.ROI, IRF_file,
            append = TRUE,
            sep = "\t", eol = "\n", col.names = FALSE, scipen = 50
        )
    }
    fwrite(list("# ref-sj.ref"), IRF_file,
        append = TRUE,
        sep = "\t", eol = "\n", col.names = FALSE, scipen = 50
    )
    fwrite(ref.sj, IRF_file,
        append = TRUE,
        sep = "\t", eol = "\n", col.names = FALSE, scipen = 50
    )

    # Add EOF (to avoid undefined behaviour when there is no termination char)
    fwrite(list("# EOF"), IRF_file,
        append = TRUE,
        sep = "\t", eol = "\n", col.names = FALSE, scipen = 50
    )
    
    gzip(filename = IRF_file, destname = paste0(IRF_file, ".gz"),
        overwrite = TRUE)
    if(file.exists(IRF_file) & file.exists(paste0(IRF_file, ".gz"))) {
        file.remove(IRF_file)
    }
    # cleanup
    if (file.exists(file.path(reference_path, "tmpdir.IntronCover.bed"))) {
        file.remove(file.path(reference_path, "tmpdir.IntronCover.bed"))
    }
    if (file.exists(file.path(reference_path, "tmpnd.IntronCover.bed"))) {
        file.remove(file.path(reference_path, "tmpnd.IntronCover.bed"))
    }
}
################################################################################
.gen_nmd <- function(reference_path, genome) {

    Exons.tr <- .gen_nmd_exons_trimmed(reference_path)
    protein.introns <- .gen_nmd_protein_introns(reference_path, Exons.tr)
    # Exclude introns preceding any ORF exons:
    # protein.introns = protein.introns[get("intron_type") != "UTR5"]
    NMD.Table <- .gen_nmd_determine(Exons.tr, protein.introns, genome, 50)
    protein.introns.red = unique(
        protein.introns[, c("intron_id", "intron_type")])
    NMD.Table[protein.introns.red,
        on = "intron_id",
        c("intron_type") := get("i.intron_type")
    ]

    write.fst(NMD.Table, file.path(reference_path, "fst", "IR.NMD.fst"))
    gc()
}

.gen_nmd_exons_trimmed <- function(reference_path) {
    Exons.tr <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Exons.fst"))
    )
    Misc <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Misc.fst"))
    )
    start.DT <- Misc[get("type") == "start_codon"]
    Exons.tr <- Exons.tr[get("transcript_id") %in% 
        start.DT[,get("transcript_id")]]

    Exons.tr[start.DT,
        on = c("transcript_id"),
        c("sc_start", "sc_end") := list(get("i.start"), get("i.end"))
    ]
    Exons.tr[
        get("start") < get("sc_start") & get("strand") == "+",
        c("start") := get("sc_start")
    ]
    Exons.tr[
        get("end") < get("sc_start") & get("strand") == "+",
        c("end") := get("sc_start")
    ]
    Exons.tr[
        get("start") > get("sc_end") & get("strand") == "-",
        c("start") := get("sc_end")
    ]
    Exons.tr[
        get("end") > get("sc_end") & get("strand") == "-",
        c("end") := get("sc_end")
    ]
    Exons.tr <- Exons.tr[get("start") < get("end")]
    return(Exons.tr)
}

.gen_nmd_protein_introns <- function(reference_path, Exons.tr) {
    candidate.introns <- as.data.table(
        read.fst(file.path(reference_path, "fst", "junctions.fst"))
    )
    Misc <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Misc.fst"))
    )
    protein.introns <- candidate.introns[
        get("transcript_id") %in% Exons.tr$transcript_id
    ]
    # determine here whether protein introns are CDS, 5' or 3' UTR introns
    UTR5 <- Misc[get("type") == "five_prime_utr"]
    UTR5.introns <- .grlGaps(
        split(
            makeGRangesFromDataFrame(as.data.frame(UTR5)),
            UTR5$transcript_id
        )
    )
    UTR5.introns <- as.data.table(UTR5.introns)
    UTR3 <- Misc[get("type") == "three_prime_utr"]
    UTR3.introns <- .grlGaps(
        split(
            makeGRangesFromDataFrame(as.data.frame(UTR3)),
            UTR3$transcript_id
        )
    )
    UTR3.introns <- as.data.table(UTR3.introns)

    CDS.introns <- .grlGaps(
        split(
            makeGRangesFromDataFrame(as.data.frame(Exons.tr)),
            Exons.tr$transcript_id
        )
    )
    CDS.introns <- as.data.table(CDS.introns)

    protein.introns[UTR5.introns,
        on = c("seqnames", "start", "end", "strand"),
        c("intron_type") := "UTR5"
    ]
    protein.introns[UTR3.introns,
        on = c("seqnames", "start", "end", "strand"),
        c("intron_type") := "UTR3"
    ]
    protein.introns[CDS.introns,
        on = c("seqnames", "start", "end", "strand"),
        c("intron_type") := "CDS"
    ]
    return(protein.introns)
}

.gen_nmd_determine <- function(exon.DT, intron.DT, genome, threshold = 50) {
    message("Calculating IR-NMD")
    exon.DT <- exon.DT[, 
        c("seqnames", "start", "end", "strand", "transcript_id")]
    exon_gr <- .grDT(exon.DT)
    message("Computing for Exon Sequences...", appendLF = FALSE)
    set(exon.DT,,"seq", as.character(getSeq(genome, exon_gr)))
    final <- .gen_nmd_determine_spliced_exon(exon.DT, intron.DT, 
        threshold = threshold)
    message("done")

    intron.DT.use = intron.DT[get("intron_type") != "UTR5"]
    exon.DT.skinny <- exon.DT[, -("seq")]
    i_partition <- c(seq(1, nrow(intron.DT.use), by = 10000), 
        nrow(intron.DT.use) + 1)
    message("Computing for Intron Sequences...")
    pb <- txtProgressBar(max = length(i_partition) - 1, style = 3)
    l_seq = 1000
    for (i in seq_len(length(i_partition) - 1)) {
        setTxtProgressBar(pb, i)
        dash_progress("Determining NMD Transcripts: Calculating...", 
            length(i_partition) - 1)
        intron.part <- intron.DT.use[
            seq(i_partition[i], i_partition[i + 1] - 1),
            c("transcript_id", "intron_id",
                "seqnames", "start", "end", "strand")
        ]
        set(intron.part,,"type", "intron")
        
        intron.part.upstream <- .gen_nmd_determine_build_introns_upstream(
            intron.part, exon.DT.skinny, use_short = TRUE)
        intron.part.upstream <- .gen_nmd_determine_retrieve_short_seq(
            exon.DT, intron.part.upstream, genome, l_seq = l_seq)
        final <- .gen_nmd_determine_translate(
            final, intron.part.upstream, use_short = TRUE, 
            threshold = threshold)
        intron_id_exclude = unique(intron.part.upstream$intron_id)
        intron_id_exclude = intron_id_exclude[(
            intron_id_exclude %in% 
            final[get("IRT_is_NMD") == TRUE, get("intron_id")]
        )]
        intron.part <- intron.part[!(get("intron_id") %in% intron_id_exclude)]

        intron.part.upstream <- .gen_nmd_determine_build_introns_upstream(
            intron.part, exon.DT.skinny, use_short = FALSE)
        intron.part.upstream <- .gen_nmd_determine_retrieve_full_seq(
            exon.DT, intron.part.upstream, genome)
        final <- .gen_nmd_determine_translate(
            final, intron.part.upstream, use_short = FALSE, 
            threshold = threshold)
    }
    setTxtProgressBar(pb, i)
    close(pb)
    message("done\n")
    return(final)
}

.gen_nmd_determine_spliced_exon <- function(
        exon.DT, intron.DT, threshold = 50) {
    exon.MLE.DT <- copy(exon.DT)
    setorderv(exon.MLE.DT, "start")
    exon.MLE.DT[, c("elem_number") := data.table::rowid(get("transcript_id"))]
    exon.MLE.DT[get("strand") == "-", 
        c("elem_number") := max(get("elem_number")) + 1 - get("elem_number"),
        by = "transcript_id"]
    exon.MLE.DT[, by = "transcript_id",
        c("is_last_elem") := (get("elem_number") == max(get("elem_number")))]
    exon.MLE.DT <- exon.MLE.DT[get("is_last_elem") == FALSE]

    # sort by order
    setorderv(exon.MLE.DT, c("transcript_id", "elem_number"))
    exon.MLE.DT <- exon.MLE.DT[, c("transcript_id", "seq")]
    splice <- exon.MLE.DT[, lapply(.SD, paste0, collapse = ""), 
        by = "transcript_id"]
    splice[is.na(rowSums(stringr::str_locate(get("seq"), "N"))),
        c("AA") := as.character(
            suppressWarnings(
                Biostrings::translate(as(get("seq"), "DNAStringSet"))
            )
        )
    ]
    # Find nucleotide position of first stop codon
    splice[, c("stop_pos") := 
        stringr::str_locate(get("AA"), "\\*")[, 1] * 3 - 2]
    splice[, c("splice_len") := nchar(get("seq"))]
    splice[!is.na(get("AA")), 
        c("stop_to_EJ") := get("splice_len") - get("stop_pos")]
    intron.DT <- intron.DT[, 
        c("seqnames", "start", "end", "strand", "transcript_id", "intron_id")]
    final <- intron.DT[, c("intron_id", "transcript_id")]
    final[splice[, c("transcript_id", "stop_pos", "splice_len", "stop_to_EJ")],
        on = "transcript_id",
        c("splice_stop_pos", "splice_start_to_last_EJ", 
                "splice_stop_to_last_EJ") :=
            list(get("i.stop_pos"), get("i.splice_len"), get("i.stop_to_EJ"))
    ]
    final[, c("splice_is_NMD") :=
        ifelse(get("splice_start_to_last_EJ") - get("splice_stop_to_last_EJ")
        >= threshold, TRUE, FALSE)]
    final[is.na(get("splice_stop_to_last_EJ")), c("splice_is_NMD") := FALSE]
    final[is.na(get("splice_start_to_last_EJ")), c("splice_is_NMD") := NA]
    return(final)
}

.gen_nmd_determine_build_introns_upstream <- function(
        intron.part, exon.DT.skinny, use_short = FALSE
) {
    intron.part.skinny = intron.part[, c("transcript_id", "intron_id")]
    # join exons with introns to determine phase of intron
    exon.DT.skinny.copy = copy(exon.DT.skinny)
    exon.DT.skinny.copy = exon.DT.skinny.copy[
        intron.part.skinny, on = "transcript_id",
        allow.cartesian=TRUE]
    exon.DT.skinny.copy = exon.DT.skinny.copy[,
        c("transcript_id", "intron_id", "seqnames", "start", "end", "strand")]
    set(exon.DT.skinny.copy, , "type", "exon")
    intron.part.upstream <- rbindlist(list(intron.part, exon.DT.skinny.copy))

    setorderv(intron.part.upstream, c("seqnames", "start"))
    intron.part.upstream[, 
        c("elem_number") := data.table::rowid(get("intron_id"))]
    intron.part.upstream[get("strand") == "-",
        c("elem_number") := 
            max(get("elem_number")) + 1 - get("elem_number"),
        by = "intron_id"
    ]

    # trim exons downstream of intron
    intron.part.upstream.intron = intron.part.upstream[get("type") == "intron"]
    intron.part.upstream[intron.part.upstream.intron,
        on = "intron_id", c("intron_pos") := get("i.elem_number")
    ]
    intron.part.upstream <- intron.part.upstream[!is.na(get("intron_pos"))]
    if(use_short) {
        intron.part.upstream <-
            intron.part.upstream[get("elem_number") < get("intron_pos") |
                get("type") == "intron"]    
    } else {
        # remove last exon: then the terminus is the last exon junction
        intron.part.upstream[,
            by = "transcript_id",
            c("is_last_elem") := (get("elem_number") == max(get("elem_number")))
        ]        
        intron.part.upstream <-
            intron.part.upstream[!get("is_last_elem") |
                get("type") == "intron"]

    }
    return(intron.part.upstream)
}

.gen_nmd_determine_retrieve_full_seq <- function(
        exon.DT, intron.part.upstream, genome
) {
    intron.part.short <- intron.part.upstream[get("type") == "intron"]
    intron.short_gr <- .grDT(intron.part.short)
    intron.part.short[, 
        c("seq") := as.character(getSeq(genome, intron.short_gr))]
    intron.part.upstream[exon.DT,
        on = c("transcript_id", "seqnames", "start", "end", "strand"),
        c("seq") := get("i.seq")
    ]
    intron.part.upstream[intron.part.short,
        on = c("intron_id", "type"),
        c("seq") := get("i.seq")
    ]
    return(intron.part.upstream)
}

.gen_nmd_determine_retrieve_short_seq <- function(
        exon.DT, intron.part.upstream, genome, l_seq = 1000, threshold = 50
) {
    intron.part.short <- intron.part.upstream[get("type") == "intron"]
    # Truncate intron by threshold as its terminus is taken as (EJC - threshold)
    intron.part.short[get("strand") == "+" & 
        get("end") - get("start") > threshold,
        c("end") := get("end") - threshold]
    intron.part.short[get("strand") == "-" & 
        get("end") - get("start") > threshold,
        c("start") := get("start") + threshold]
    intron.part.short[
        get("strand") == "+" & get("end") - get("start") > l_seq,
        c("end") := get("start") + l_seq
    ]
    intron.part.short[
        get("strand") == "-" & get("end") - get("start") > l_seq,
        c("start") := get("end") - l_seq
    ]

    intron.short_gr <- .grDT(intron.part.short)
    intron.part.short[, 
        c("seq") := as.character(getSeq(genome, intron.short_gr))]
    intron.part.upstream[exon.DT,
        on = c("transcript_id", "seqnames", "start", "end", "strand"),
        c("seq") := get("i.seq")
    ]
    intron.part.upstream[intron.part.short,
        on = c("intron_id", "type"),
        c("seq") := get("i.seq")
    ]
    return(intron.part.upstream)
}

.gen_nmd_determine_translate <- function(
        splice_table, elems, use_short = FALSE, threshold = 50
) {
    setorderv(elems, c("transcript_id", "elem_number"))
    elems <- elems[, c("intron_id", "seq")]

    IRT <- elems[, 
        lapply(.SD, paste0, collapse = ""), by = "intron_id"]
    # trim
    IRT[, c("seq") := substr(get("seq"), 1, 
        nchar(get("seq")) - (nchar(get("seq")) %% 3))]
    IRT[
        is.na(rowSums(stringr::str_locate(get("seq"), "N"))),
        c("AA") := as.character(
            # suppressWarnings(
                Biostrings::translate(as(get("seq"), "DNAStringSet"))
            # )
        )
    ]

    # Find nucleotide position of first stop codon
    IRT[, c("stop_pos") := 
        stringr::str_locate(get("AA"), "\\*")[, 1] * 3 - 2]
    IRT[, c("IRT_len") := nchar(get("seq"))]
    IRT[!is.na(get("AA")), c("stop_to_EJ") := 
        get("IRT_len") - get("stop_pos")]
    IRT[, c("use_short") := use_short]

    IRT[, c("IRT_is_NMD") := ifelse(
        get("stop_to_EJ") >= get("threshold"), TRUE, FALSE)]
    IRT[is.na(get("stop_pos")), c("IRT_is_NMD") := FALSE]
    IRT[is.na(get("IRT_len")), c("IRT_is_NMD") := NA]

    # Annotate into splice_table
    splice_table[IRT,
        on = "intron_id",
        c(
            "IRT_stop_pos", "IRT_start_to_last_EJ", "IRT_stop_to_last_EJ",
            "IRT_use_short", "IRT_is_NMD"
        ) :=
            list(
                get("i.stop_pos"), get("i.IRT_len"), get("i.stop_to_EJ"),
                get("i.use_short"), get("i.IRT_is_NMD")
            )
    ]
    return(splice_table)
}

################################################################################
# Sub

.gen_splice <- function(reference_path, genome) {
    message("Annotating Splice Events\n")
    candidate.introns <- as.data.table(
        read.fst(file.path(reference_path, "fst", "junctions.fst"))
    )
    introns.skipcoord <- .gen_splice_skipcoord(
        reference_path, candidate.introns)

    # chrOrder <- names(seqinfo(genome))
    message("Annotating Mutually-Exclusive-Exon Splice Events...",
        appendLF = FALSE
    )
    introns_found_MXE <- .gen_splice_MXE(introns.skipcoord)
    message("done")

    # annotate skipped junctions with two included junctions
    message("Annotating Skipped-Exon Splice Events...", appendLF = FALSE)
    introns_found_SE <- .gen_splice_SE(introns.skipcoord, candidate.introns)
    message("done")

    message("Annotating Alternate 5' / 3' Splice Site Splice Events...",
        appendLF = FALSE)

    introns_found_A5SS = .gen_splice_A5SS(candidate.introns)
    introns_found_A3SS = .gen_splice_A3SS(candidate.introns)
    message("done")

    message("Annotating Alternate First / Last Exon Splice Events...",
        appendLF = FALSE)
    # AFE/ALE

    introns_found_AFE = .gen_splice_AFE(candidate.introns, introns_found_A5SS)
    introns_found_ALE = .gen_splice_ALE(candidate.introns, introns_found_A3SS)
    message("done")
    gc()

################################################################################
    #   Filter for valid splicing

    is_valid_splice_type <- function(x) !is.null(x) && nrow(x) > 0
    tmp_AS <- list(
        introns_found_MXE, introns_found_SE,
        introns_found_AFE, introns_found_ALE,
        introns_found_A5SS, introns_found_A3SS
    )
    tmp_AS <- base::Filter(is_valid_splice_type, tmp_AS)
    AS_Table <- rbindlist(tmp_AS)

    if (nrow(AS_Table) > 0) {
        .gen_splice_save(AS_Table, candidate.introns, reference_path)
        message("Splice Annotations Filtered\n")
    } else {
        message("No splice events found\n")
    }
}

################################################################################

.gen_splice_skipcoord <- function(reference_path, candidate.introns) {
    GeneOrder <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Genes.fst"))
    )
    setorder(GeneOrder, seqnames, start, end, strand)    
    
    introns.skipcoord <- copy(candidate.introns)
    introns.skipcoord[, c("gene_id") := 
        factor(get("gene_id"), GeneOrder[, get("gene_id")], ordered = TRUE)]
        
    setorderv(introns.skipcoord, 
        c("gene_id", "transcript_name", "intron_number"))

    introns.skipcoord[, c("skip_coord") := ""]    
    introns.skipcoord[get("strand") == "+",
        c("skip_coord") := paste0(
                get("seqnames"), ":", get("intron_start"),
                "-", data.table::shift(get("intron_end"), 1, NA, "lead"),
                "/", get("strand")
        ),
        by = "transcript_id"
    ]
    introns.skipcoord[get("strand") == "-",
        c("skip_coord") := paste0(
                get("seqnames"), ":",
                data.table::shift(get("intron_start"), 1, NA, "lead"),
                "-", get("intron_end"), "/", get("strand")
        ),
        by = "transcript_id"
    ]
    introns.skipcoord[grepl("NA", get("skip_coord")), c("skip_coord") := NA]    

    introns.skipcoord[, c("skip_coord_2") := 
        data.table::shift(get("skip_coord"), 1, NA, "lag")    ]
    return(introns.skipcoord)
}

.gen_splice_MXE <- function(introns.skipcoord) {
    introns_search_MXE <- introns.skipcoord[introns.skipcoord[,
        .I[get("intron_number") < max(get("intron_number"))],
        by = "transcript_id"]$V1]
    introns_search_MXE <- introns_search_MXE[
        introns_search_MXE[, .N, by = c("gene_id", "skip_coord")],
        on = c("gene_id", "skip_coord"), 
        c("N") := get("i.N")]
    introns_search_MXE <- introns_search_MXE[get("N") > 1]
    introns_search_MXE_pos <- introns_search_MXE[get("strand") == "+"]
    setorderv(introns_search_MXE_pos, 
        c("seqnames", "intron_start", "intron_end"))
    introns_search_MXE_neg <- introns_search_MXE[get("strand") == "-"]
    setorderv(introns_search_MXE_neg, 
        c("seqnames", "intron_end", "intron_start"),
        order = c(1, 1, -1))
    introns_search_MXE <- rbindlist(
        list(introns_search_MXE_pos, introns_search_MXE_neg))
    introns_search_MXE <- introns_search_MXE[,
        c("skip_coord", "gene_id", "Event", "transcript_id",
            "transcript_name", "intron_number")]
    setnames(introns_search_MXE, old = "Event", new = "Event1")

    introns_search_MXE2 <- introns.skipcoord[,
        c("skip_coord_2", "gene_id", "Event", 
        "transcript_id", "transcript_name")]
    setnames(introns_search_MXE2, old = c("skip_coord_2", "Event"),
        new = c("skip_coord", "Event2"))

    introns_search_MXE[introns_search_MXE2,
        on = c("gene_id", "transcript_id", "transcript_name", "skip_coord"),
        c("Event2") := get("i.Event2")]

    introns_search_MXE <- unique(introns_search_MXE,
        by = c("gene_id", "skip_coord", "Event1"))
    introns_search_MXE <- unique(introns_search_MXE,
        by = c("gene_id", "skip_coord", "Event2"))
    introns_search_MXE <- introns_search_MXE[, if (.N > 1) .SD,
        by = c("gene_id", "skip_coord")]

    if (nrow(introns_search_MXE) > 0) {
        introns_found_MXE <- introns_search_MXE[,
            {
                edge1 <- rep(seq_len(.N), (.N:1) - 1L)
                i <- 2L:(.N * (.N - 1L) / 2L + 1L)
                o <- cumsum(c(0, (.N - 2L):1))
                edge2 <- i - o[edge1]
                .(
                    gene_id = get("gene_id")[edge1],
                    gene_id_b = get("gene_id")[edge2],
                    Event1a = get("Event1")[edge1],
                    Event1b = get("Event1")[edge2],
                    Event2a = get("Event2")[edge1],
                    Event2b = get("Event2")[edge2],
                    transcript_id_a = get("transcript_id")[edge1],
                    transcript_id_b = get("transcript_id")[edge2],
                    transcript_name_a = get("transcript_name")[edge1],
                    transcript_name_b = get("transcript_name")[edge2],
                    intron_number_a = get("intron_number")[edge1],
                    intron_number_b = get("intron_number")[edge2]
                )
            },by = "skip_coord"
        ]
        setorderv(introns_found_MXE, c("gene_id", "transcript_name_a"))
        introns_found_MXE[, c("EventName") := paste0(
                "MXE:", get("transcript_name_a"), "-exon",
                (1 + get("intron_number_a")), ";",
                get("transcript_name_b"), "-exon", 
                (1 + get("intron_number_b")))
        ]
        introns_found_MXE[, c("EventID") := paste0("MXE#", seq_len(.N))]
        setnames(introns_found_MXE,
            old = "skip_coord", new = "EventRegion")
        introns_found_MXE[, c("EventType") := "MXE"]
        introns_found_MXE <- introns_found_MXE[,
            c("EventType", "EventID", "EventName", "Event1a", "Event1b",
                "Event2a", "Event2b", "gene_id", "gene_id_b", "EventRegion",
                "transcript_id_a", "transcript_name_a", "intron_number_a",
                "transcript_id_b", "transcript_name_b", "intron_number_b")]
        introns_found_MXE <- unique(introns_found_MXE,
            by = c("Event1a", "Event1b", "Event2a", "Event2b"))
    } else {
        introns_found_MXE <- c()
    }
    return(introns_found_MXE)
}

.gen_splice_SE <- function(introns.skipcoord, candidate.introns) {
    introns.skippedJn <- introns.skipcoord[
        get("skip_coord") %in% get("Event"),
        c("gene_id", "gene_name", "skip_coord")]
    introns.skippedJn <- unique(introns.skippedJn)
    
    introns_found_SE <- introns.skippedJn[, "skip_coord"]
    introns_search_SE <- candidate.introns[,
        c("gene_id", "Event", "transcript_id",
            "transcript_name", "intron_number")]
    setnames(introns_search_SE,
        old = c("Event", "transcript_id", "transcript_name", "intron_number"),
        new = c("skip_coord", "skip_transcript_id", "skip_transcript_name",
            "skip_intron_number"))
    introns_found_SE[introns_search_SE, on = "skip_coord",
        c("gene_id_b", "skip_transcript_id",
            "skip_transcript_name", "skip_intron_number") :=
            list(get("i.gene_id"), get("i.skip_transcript_id"),
                get("i.skip_transcript_name"), get("i.skip_intron_number"))
    ]
    introns_found_SE <- unique(introns_found_SE,
        by = c("gene_id_b", "skip_coord"))

    introns_search_SE2 <- introns.skipcoord[,
        c("skip_coord", "gene_id", "Event", "transcript_id", 
            "transcript_name", "intron_number")]
    introns_found_SE[introns_search_SE2,
        on = "skip_coord",
        c("gene_id", "inc_coord_upst", "inc_transcript_id",
            "inc_transcript_name", "inc_intron_number") :=
            list(get("i.gene_id"), get("i.Event"), get("i.transcript_id"),
                get("i.transcript_name"), get("i.intron_number"))]

    introns_search_SE3 <- introns.skipcoord[,
        c("skip_coord_2", "gene_id", "Event",
            "transcript_id", "transcript_name")]
    setnames(introns_search_SE3,
        old = c("skip_coord_2", "transcript_id"),
        new = c("skip_coord", "inc_transcript_id"))
    introns_found_SE[introns_search_SE3,
        on = c("skip_coord", "inc_transcript_id"),
        c("inc_coord_downst") := get("i.Event")]
    introns_found_SE <- unique(introns_found_SE,
        by = c("gene_id", "skip_coord",
            "inc_transcript_id", "inc_transcript_name"))
    setorderv(introns_found_SE, c("gene_id", "inc_transcript_name"))
    introns_found_SE[, c("EventName") := paste0(
        "SE:", get("inc_transcript_name"), "-exon", 
        (1 + get("inc_intron_number")), ";",
        get("skip_transcript_name"), "-int", get("skip_intron_number"))]
    introns_found_SE[, c("EventID") := paste0("SE#", seq_len(.N))]
    introns_found_SE[, c("EventType") := "SE"]
    introns_found_SE[, c("Event2b") := NA]
    introns_found_SE[, c("EventRegion") := get("skip_coord")]
    introns_found_SE <- introns_found_SE[,
        c("EventType", "EventID", "EventName", "inc_coord_upst",
            "skip_coord", "inc_coord_downst", "Event2b",
            "gene_id", "gene_id_b", "EventRegion",
            "inc_transcript_id", "inc_transcript_name", "inc_intron_number",
            "skip_transcript_id", "skip_transcript_name", "skip_intron_number")]
    setnames(introns_found_SE,
        old = c("inc_coord_upst", "inc_coord_downst", "skip_coord"),
        new = c("Event1a", "Event2a", "Event1b"))
    setnames(introns_found_SE,
        old = c("inc_transcript_id", "inc_transcript_name",
            "skip_transcript_id", "skip_transcript_name"),
        new = c("transcript_id_a", "transcript_name_a",
            "transcript_id_b", "transcript_name_b"))
    setnames(introns_found_SE,
        new = c("intron_number_a", "intron_number_b"),
        old = c("inc_intron_number", "skip_intron_number"))
    return(introns_found_SE)
}

.gen_splice_AFE <- function(candidate.introns, introns_found_A5SS) {
    introns_search_AFE <- candidate.introns[get("intron_number") == 1]
    introns_search_AFE_pos <- introns_search_AFE[get("strand") == "+"]
    setorderv(introns_search_AFE_pos,
        c("seqnames", "intron_end", "intron_start"),
        order = c(1, 1, -1))
    introns_search_AFE_pos <- introns_search_AFE_pos[,
        c("seqnames", "intron_end", "strand", "Event", "gene_id",
            "transcript_id", "transcript_name", "intron_number")]
    setnames(introns_search_AFE_pos, old = "intron_end", new = "intron_coord")

    introns_search_AFE_neg <- introns_search_AFE[get("strand") == "-"]
    setorderv(introns_search_AFE_neg, 
        c("seqnames", "intron_start", "intron_end"))
    introns_search_AFE_neg <- introns_search_AFE_neg[,
        c("seqnames", "intron_start", "strand", "Event", "gene_id",
            "transcript_id", "transcript_name", "intron_number")]
    setnames(introns_search_AFE_neg, old = "intron_start", new = "intron_coord")

    introns_search_AFE <- rbindlist(
        list(introns_search_AFE_pos, introns_search_AFE_neg))
    introns_search_AFE <- unique(introns_search_AFE, by = "Event")

    introns_found_AFE <- introns_search_AFE[,
        {
            edge1 <- rep(seq_len(.N), (.N:1) - 1L)
            i <- 2L:(.N * (.N - 1L) / 2L + 1L)
            o <- cumsum(c(0, (.N - 2L):1))
            edge2 <- i - o[edge1]
            .(
                gene_id = get("gene_id")[edge1],
                gene_id_b = get("gene_id")[edge2],
                Event1a = get("Event")[edge1],
                Event1b = get("Event")[edge2],
                Event2a = NA, Event2b = NA, EventRegion = get("Event")[edge2],
                transcript_id_a = get("transcript_id")[edge1],
                transcript_id_b = get("transcript_id")[edge2],
                transcript_name_a = get("transcript_name")[edge1],
                transcript_name_b = get("transcript_name")[edge2],
                intron_number_a = get("intron_number")[edge1],
                intron_number_b = get("intron_number")[edge2]
            )
        }, by = c("seqnames", "intron_coord")
    ]
    introns_found_AFE <- introns_found_AFE[!is.na(get("gene_id"))]
    introns_found_AFE <- unique(introns_found_AFE, by = c("Event1a", "Event1b"))

    setorderv(introns_found_AFE, c("gene_id", "transcript_name_a"))
    introns_found_AFE <- introns_found_AFE[,
        c("EventName") := paste0(
            "AFE:", get("transcript_name_a"), "-exon",
            get("intron_number_a"), ";",
            get("transcript_name_b"), "-exon", get("intron_number_b")
        )]
    introns_found_AFE <- introns_found_AFE[, c("EventType") := "AFE"]
    introns_found_AFE <- introns_found_AFE[, c("EventRegion") := get("Event1b")]

    introns_found_AFE <- introns_found_AFE[!introns_found_A5SS,
        on = c("Event1a", "Event1b")]
        
    introns_found_AFE[, c("EventID") := paste0("AFE#", seq_len(.N))]
    introns_found_AFE <- introns_found_AFE[,
        c("EventType", "EventID", "EventName",
            "Event1a", "Event1b", "Event2a", "Event2b",
            "gene_id", "gene_id_b", "EventRegion",
            "transcript_id_a", "transcript_name_a", "intron_number_a",
            "transcript_id_b", "transcript_name_b", "intron_number_b")]
    return(introns_found_AFE)
}
.gen_splice_ALE <- function(candidate.introns, introns_found_A3SS) {
    introns_search_ALE <- candidate.introns[candidate.introns[,
        .I[get("intron_number") == max(get("intron_number"))],
        by = "transcript_id"]$V1]
    introns_search_ALE_pos <- introns_search_ALE[get("strand") == "+"]
    setorderv(introns_search_ALE_pos, 
        c("seqnames", "intron_start", "intron_end"))
    introns_search_ALE_pos <- introns_search_ALE_pos[,
        c("seqnames", "intron_start", "strand", "Event", "gene_id",
            "transcript_id", "transcript_name", "intron_number")]
    setnames(introns_search_ALE_pos,
        old = "intron_start", new = "intron_coord")

    introns_search_ALE_neg <- introns_search_ALE[get("strand") == "-"]
    setorderv(introns_search_ALE_neg, 
        c("seqnames", "intron_end", "intron_start"),
        order = c(1, 1, -1))
    introns_search_ALE_neg <- introns_search_ALE_neg[,
        c("seqnames", "intron_end", "strand", "Event", "gene_id",
            "transcript_id", "transcript_name", "intron_number")]
    setnames(introns_search_ALE_neg, old = "intron_end", new = "intron_coord")

    introns_search_ALE <- rbindlist(
        list(introns_search_ALE_pos, introns_search_ALE_neg))
    introns_search_ALE <- unique(introns_search_ALE, by = "Event")

    introns_found_ALE <- introns_search_ALE[,
        {
            edge1 <- rep(seq_len(.N), (.N:1) - 1L)
            i <- 2L:(.N * (.N - 1L) / 2L + 1L)
            o <- cumsum(c(0, (.N - 2L):1))
            edge2 <- i - o[edge1]
            .(
                gene_id = get("gene_id")[edge1],
                gene_id_b = get("gene_id")[edge2],
                Event1a = get("Event")[edge1],
                Event1b = get("Event")[edge2],
                Event2a = NA, Event2b = NA, EventRegion = get("Event")[edge2],
                transcript_id_a = get("transcript_id")[edge1],
                transcript_id_b = get("transcript_id")[edge2],
                transcript_name_a = get("transcript_name")[edge1],
                transcript_name_b = get("transcript_name")[edge2],
                intron_number_a = get("intron_number")[edge1],
                intron_number_b = get("intron_number")[edge2]
            )
        }, by = c("seqnames", "intron_coord")
    ]
    introns_found_ALE <- introns_found_ALE[!is.na(get("gene_id"))]
    introns_found_ALE <- unique(introns_found_ALE, by = c("Event1a", "Event1b"))

    setorderv(introns_found_ALE, c("gene_id", "transcript_name_a"))
    introns_found_ALE <- introns_found_ALE[,
        c("EventName") := paste0("ALE:", get("transcript_name_a"), "-exon",
            get("intron_number_a"), ";",
            get("transcript_name_b"), "-exon", get("intron_number_b"))]
    introns_found_ALE <- introns_found_ALE[, c("EventType") := "ALE"]
    introns_found_ALE <- introns_found_ALE[, c("EventRegion") := get("Event1b")]

    introns_found_ALE <- introns_found_ALE[!introns_found_A3SS,
        on = c("Event1a", "Event1b")]

    introns_found_ALE[, c("EventID") := paste0("ALE#", seq_len(.N))]
    introns_found_ALE <- introns_found_ALE[,
        c("EventType", "EventID", "EventName",
            "Event1a", "Event1b", "Event2a", "Event2b",
            "gene_id", "gene_id_b", "EventRegion",
            "transcript_id_a", "transcript_name_a", "intron_number_a",
            "transcript_id_b", "transcript_name_b", "intron_number_b")]
    return(introns_found_ALE)
}
.gen_splice_ASS_common <- function(candidate.introns) {
    candidate.introns.ASS <- candidate.introns[
        !is.na(get("exon_group_stranded_upstream")) &
            !is.na(get("exon_group_stranded_downstream"))]
    setnames(candidate.introns.ASS,
        old = c("exon_group_stranded_upstream",
            "exon_group_stranded_downstream"),
        new = c("exon_groups_start", "exon_groups_end"))
    return(candidate.introns.ASS)
}

.gen_splice_A5SS <- function(candidate.introns) {
    introns_search_A5SS <- copy(.gen_splice_ASS_common(candidate.introns))
    introns_search_A5SS_pos <- introns_search_A5SS[get("strand") == "+"]
    
    setorderv(introns_search_A5SS_pos, 
        c("seqnames", "intron_end", "intron_start"),
        order = c(1, 1, -1))
    introns_search_A5SS_pos <- introns_search_A5SS_pos[,
        c("seqnames", "intron_end", "strand", "Event", "gene_id",
            "transcript_id", "transcript_name", "intron_number",
            "exon_groups_start", "exon_groups_end")]
    setnames(introns_search_A5SS_pos,
        old = "intron_end", new = "intron_coord")

    introns_search_A5SS_neg <- introns_search_A5SS[get("strand") == "-"]
    setorderv(introns_search_A5SS_neg, 
        c("seqnames", "intron_start", "intron_end"))
    introns_search_A5SS_neg <- introns_search_A5SS_neg[,
        c("seqnames", "intron_start", "strand", "Event", "gene_id",
            "transcript_id", "transcript_name", "intron_number",
            "exon_groups_start", "exon_groups_end")]
    setnames(introns_search_A5SS_neg,
        old = "intron_start", new = "intron_coord")

    introns_search_A5SS <- rbindlist(
        list(introns_search_A5SS_pos, introns_search_A5SS_neg))
    introns_search_A5SS <- unique(introns_search_A5SS, by = "Event")
    introns_found_A5SS <- introns_search_A5SS[,
        {
            edge1 <- rep(seq_len(.N), (.N:1) - 1L)
            i <- 2L:(.N * (.N - 1L) / 2L + 1L)
            o <- cumsum(c(0, (.N - 2L):1))
            edge2 <- i - o[edge1]
            .(
                gene_id = get("gene_id")[edge1],
                gene_id_b = get("gene_id")[edge2],
                Event1a = get("Event")[edge1],
                Event1b = get("Event")[edge2],
                Event2a = NA, Event2b = NA, EventRegion = get("Event")[edge2],
                transcript_id_a = get("transcript_id")[edge1],
                transcript_id_b = get("transcript_id")[edge2],
                transcript_name_a = get("transcript_name")[edge1],
                transcript_name_b = get("transcript_name")[edge2],
                intron_number_a = get("intron_number")[edge1],
                intron_number_b = get("intron_number")[edge2],
                exon_groups_start_a = get("exon_groups_start")[edge1],
                exon_groups_start_b = get("exon_groups_start")[edge2],
                exon_groups_end_a = get("exon_groups_end")[edge1],
                exon_groups_end_b = get("exon_groups_end")[edge2]
            )
        }, by = "intron_coord"
    ]
    introns_found_A5SS <- introns_found_A5SS[!is.na(get("gene_id"))]
    introns_found_A5SS <- unique(introns_found_A5SS,
        by = c("Event1a", "Event1b"))
    # filter by same exon group starts and ends:
    introns_found_A5SS <- introns_found_A5SS[
        get("exon_groups_start_a") == get("exon_groups_start_b")]
    introns_found_A5SS <- introns_found_A5SS[
        get("exon_groups_end_a") == get("exon_groups_end_b")]
    setorderv(introns_found_A5SS, c("gene_id", "transcript_name_a"))
    introns_found_A5SS <- introns_found_A5SS[,
        c("EventName") := paste0(
            "A5SS:", get("transcript_name_a"), "-exon",
            get("intron_number_a"), ";",
            get("transcript_name_b"), "-exon", get("intron_number_b"))]
    introns_found_A5SS <- introns_found_A5SS[, c("EventType") := "A5SS"]
    introns_found_A5SS <- introns_found_A5SS[, 
        c("EventRegion") := get("Event1b")]
    # introns_found_A5SS <- introns_found_A5SS[!introns_found_AFE,
        # on = c("Event1a", "Event1b")]
    introns_found_A5SS[, c("EventID") := paste0("A5SS#", seq_len(.N))]
    introns_found_A5SS <- introns_found_A5SS[,
        c("EventType", "EventID", "EventName",
            "Event1a", "Event1b", "Event2a", "Event2b",
            "gene_id", "gene_id_b", "EventRegion",
            "transcript_id_a", "transcript_name_a", "intron_number_a",
            "transcript_id_b", "transcript_name_b", "intron_number_b")]
    return(introns_found_A5SS)
}

.gen_splice_A3SS <- function(candidate.introns) {
    introns_search_A3SS <- copy(
        .gen_splice_ASS_common(candidate.introns))
    introns_search_A3SS_pos <- introns_search_A3SS[get("strand") == "+"]
    setorderv(introns_search_A3SS_pos, 
        c("seqnames", "intron_start", "intron_end"))
    introns_search_A3SS_pos <- introns_search_A3SS_pos[,
        c("seqnames", "intron_start", "strand", "Event", "gene_id",
            "transcript_id", "transcript_name", "intron_number",
            "exon_groups_start", "exon_groups_end")]
    setnames(introns_search_A3SS_pos,
        old = "intron_start", new = "intron_coord")

    introns_search_A3SS_neg <- introns_search_A3SS[get("strand") == "-"]
    setorderv(introns_search_A3SS_neg, 
        c("seqnames", "intron_end", "intron_start"),
        order = c(1, 1, -1))
    introns_search_A3SS_neg <- introns_search_A3SS_neg[,
        c("seqnames", "intron_end", "strand", "Event", "gene_id",
            "transcript_id", "transcript_name", "intron_number",
            "exon_groups_start", "exon_groups_end")]
    setnames(introns_search_A3SS_neg,
        old = "intron_end", new = "intron_coord")

    introns_search_A3SS <- rbindlist(
        list(introns_search_A3SS_pos, introns_search_A3SS_neg))
    introns_search_A3SS <- unique(introns_search_A3SS, by = "Event")

    introns_found_A3SS <- introns_search_A3SS[,
        {
            edge1 <- rep(seq_len(.N), (.N:1) - 1L)
            i <- 2L:(.N * (.N - 1L) / 2L + 1L)
            o <- cumsum(c(0, (.N - 2L):1))
            edge2 <- i - o[edge1]
            .(
                gene_id = get("gene_id")[edge1],
                gene_id_b = get("gene_id")[edge2],
                Event1a = get("Event")[edge1],
                Event1b = get("Event")[edge2],
                Event2a = NA, Event2b = NA, EventRegion = get("Event")[edge2],
                transcript_id_a = get("transcript_id")[edge1],
                transcript_id_b = get("transcript_id")[edge2],
                transcript_name_a = get("transcript_name")[edge1],
                transcript_name_b = get("transcript_name")[edge2],
                intron_number_a = get("intron_number")[edge1],
                intron_number_b = get("intron_number")[edge2],
                exon_groups_start_a = get("exon_groups_start")[edge1],
                exon_groups_start_b = get("exon_groups_start")[edge2],
                exon_groups_end_a = get("exon_groups_end")[edge1],
                exon_groups_end_b = get("exon_groups_end")[edge2]
            )
        }, by = "intron_coord"
    ]
    introns_found_A3SS <- introns_found_A3SS[!is.na(get("gene_id"))]
    introns_found_A3SS <- unique(introns_found_A3SS,
        by = c("Event1a", "Event1b"))
    # filter by same exon group starts and ends:
    introns_found_A3SS <- introns_found_A3SS[
        get("exon_groups_start_a") == get("exon_groups_start_b")]
    introns_found_A3SS <- introns_found_A3SS[
        get("exon_groups_end_a") == get("exon_groups_end_b")]

    setorderv(introns_found_A3SS, c("gene_id", "transcript_name_a"))
    introns_found_A3SS <- introns_found_A3SS[,
        c("EventName") := paste0(
            "A3SS:", get("transcript_name_a"), "-exon",
            get("intron_number_a"), ";",
            get("transcript_name_b"), "-exon", get("intron_number_b"))
    ]
    introns_found_A3SS <- introns_found_A3SS[, c("EventType") := "A3SS"]
    introns_found_A3SS <- introns_found_A3SS[, 
        c("EventRegion") := get("Event1b")]
    # introns_found_A3SS <- introns_found_A3SS[!introns_found_ALE,
        # on = c("Event1a", "Event1b")]
    introns_found_A3SS[, c("EventID") := paste0("A3SS#", seq_len(.N))]
    introns_found_A3SS <- introns_found_A3SS[,
        c("EventType", "EventID", "EventName",
            "Event1a", "Event1b", "Event2a", "Event2b",
            "gene_id", "gene_id_b", "EventRegion",
            "transcript_id_a", "transcript_name_a", "intron_number_a",
            "transcript_id_b", "transcript_name_b", "intron_number_b")]
    return(introns_found_A3SS)
}

.gen_splice_save <- function(AS_Table, candidate.introns, reference_path) {
    candidate.introns.order <- copy(candidate.introns)
    if (!("transcript_support_level" %in% colnames(candidate.introns))) {
        candidate.introns.order[, c("transcript_support_level") := "NA"]
    }
    candidate.introns.order[,
        c("is_protein_coding") := !is.na(get("protein_id"))]
    candidate.introns.order[, by = "transcript_id",
        c("is_last_intron") :=
            (get("intron_number") == max(get("intron_number")))]

    AS_Table <- .gen_splice_prep_events(AS_Table, candidate.introns.order,
        reference_path)
    AS_Table <- .gen_splice_name_events(AS_Table, reference_path)
}

.gen_splice_prep_events <- function(AS_Table, candidate.introns.order,
        reference_path) {
    AS_Table_search.a <- AS_Table[,
        c("EventType", "EventID", "Event1a", "Event2a")]
    AS_Table_search.a[, c("Event") := get("Event1a")]
    AS_Table_search.a <- candidate.introns.order[AS_Table_search.a,
        on = "Event",
        c("EventType", "EventID", "Event1a", "Event2a", "transcript_id",
            "transcript_support_level", "is_protein_coding",
            "is_last_intron", "intron_number")]
    setnames(AS_Table_search.a, "intron_number", "in_1a")
    AS_Table_search.a <-
        AS_Table_search.a[get("EventType") != "AFE" | get("in_1a") == 1]
    AS_Table_search.a <-
        AS_Table_search.a[get("EventType") != "ALE" | get("is_last_intron")]
    AS_Table_search.a[, c("Event") := get("Event2a")]
    AS_Table_search.a[is.na(get("Event")), c("Event") := get("Event1a")]
    AS_Table_search.a <- candidate.introns.order[AS_Table_search.a,
        on = c("Event", "transcript_id", "transcript_support_level"),
        c("EventType", "EventID", "Event1a", "Event2a", "transcript_id",
            "transcript_support_level", "is_protein_coding",
            "is_last_intron", "in_1a", "intron_number")]
    AS_Table_search.a <- AS_Table_search.a[!is.na(get("intron_number"))]
    setnames(AS_Table_search.a, "intron_number", "in_2a")

    AS_Table_search.b <- AS_Table[,
        c("EventType", "EventID", "Event1b", "Event2b")]
    AS_Table_search.b[, c("Event") := get("Event1b")]
    AS_Table_search.b <- candidate.introns.order[AS_Table_search.b,
        on = "Event",
        c("EventType", "EventID", "Event1b", "Event2b", "transcript_id",
            "transcript_support_level", "is_protein_coding",
            "is_last_intron", "intron_number")]
    setnames(AS_Table_search.b, "intron_number", "in_1b")
    AS_Table_search.b <-
        AS_Table_search.b[get("EventType") != "AFE" | get("in_1b") == 1]
    AS_Table_search.b <-
        AS_Table_search.b[get("EventType") != "ALE" | get("is_last_intron")]
    AS_Table_search.b[, c("Event") := get("Event2b")]
    AS_Table_search.b[is.na(get("Event")), c("Event") := get("Event1b")]
    AS_Table_search.b <- candidate.introns.order[AS_Table_search.b,
        on = c("Event", "transcript_id", "transcript_support_level"),
        c("EventType", "EventID", "Event1b", "Event2b", "transcript_id",
            "transcript_support_level", "is_protein_coding",
            "is_last_intron", "in_1b", "intron_number")]
    AS_Table_search.b <- AS_Table_search.b[!is.na(get("intron_number"))]
    setnames(AS_Table_search.b, "intron_number", "in_2b")

    AS_Table_search.a[candidate.introns.order, on = "transcript_id",
        c("transcript_name") := get("i.transcript_name")]
    AS_Table_search.b[candidate.introns.order, on = "transcript_id",
        c("transcript_name") := get("i.transcript_name")]
    setorderv(AS_Table_search.a,
        c("transcript_support_level", "is_protein_coding", 
            "transcript_name"),
        order = c(1, -1, 1))
    setorderv(AS_Table_search.b,
        c("transcript_support_level", "is_protein_coding", 
            "transcript_name"),
        order = c(1, -1, 1))
    AS_Table.find.a <- unique(AS_Table_search.a, by = "EventID")
    AS_Table.find.a <- AS_Table.find.a[AS_Table[, "EventID"], 
        on = "EventID"]
    AS_Table.find.b <- unique(AS_Table_search.b, by = "EventID")
    AS_Table.find.b <- AS_Table.find.b[AS_Table[, "EventID"], 
        on = "EventID"]

    AS_Table$transcript_id_a <- AS_Table.find.a$transcript_id
    AS_Table$transcript_name_a <- AS_Table.find.a$transcript_name
    AS_Table$intron_number_a <- AS_Table.find.a$in_1a
    AS_Table$transcript_id_b <- AS_Table.find.b$transcript_id
    AS_Table$transcript_name_b <- AS_Table.find.b$transcript_name
    AS_Table$intron_number_b <- AS_Table.find.b$in_1b

    setnames(AS_Table_search.a,
        old = c("Event1a", "Event2a", "in_1a", "in_2a"),
        new = c("Event1", "Event2", "in_1", "in_2"))
    AS_Table_search.a[, c("isoform") := "A"]
    setnames(AS_Table_search.b,
        old = c("Event1b", "Event2b", "in_1b", "in_2b"),
        new = c("Event1", "Event2", "in_1", "in_2"))
    AS_Table_search.b[, c("isoform") := "B"]

    write.fst(as.data.frame(rbind(AS_Table_search.a, AS_Table_search.b)),
        file.path(reference_path, "fst", "Splice.options.fst"))
        
    return(AS_Table)
}

.gen_splice_name_events <- function(AS_Table, reference_path) {
    AS_Table[get("EventType") == "MXE",
        c("EventName") := paste0(
            "MXE:", get("transcript_name_a"), "-exon",
            as.character(as.numeric(get("intron_number_a")) + 1), ";",
            get("transcript_name_b"), "-exon",
            as.character(as.numeric(get("intron_number_b")) + 1))]
    AS_Table[
        get("EventType") == "SE",
        c("EventName") := paste0(
            "SE:", get("transcript_name_a"), "-exon",
            as.character(as.numeric(get("intron_number_a")) + 1), ";",
            get("transcript_name_b"), "-int",
            as.character(as.numeric(get("intron_number_b"))))]
    AS_Table[
        get("EventType") == "AFE",
        c("EventName") := paste0(
            "AFE:", get("transcript_name_a"), "-exon1;",
            get("transcript_name_b"), "-exon1")]
    AS_Table[
        get("EventType") == "ALE",
        c("EventName") := paste0(
            "ALE:", get("transcript_name_a"), "-exon",
            as.character(as.numeric(get("intron_number_a")) + 1), ";",
            get("transcript_name_b"), "-exon",
            as.character(as.numeric(get("intron_number_b")) + 1))]
    AS_Table[
        get("EventType") == "A5SS",
        c("EventName") := paste0(
            "A5SS:", get("transcript_name_a"), "-exon",
            as.character(as.numeric(get("intron_number_a"))), ";",
            get("transcript_name_b"), "-exon",
            as.character(as.numeric(get("intron_number_b"))))]
    AS_Table[
        get("EventType") == "A3SS",
        c("EventName") := paste0(
            "A3SS:", get("transcript_name_a"), "-exon",
            as.character(as.numeric(get("intron_number_a") + 1)), ";",
            get("transcript_name_b"), "-exon",
            as.character(as.numeric(get("intron_number_b") + 1)))]

    write.fst(as.data.frame(AS_Table), 
        file.path(reference_path, "fst", "Splice.fst"))
        
    return(AS_Table)
}

################################################################################
# Sub
.gen_splice_proteins <- function(reference_path, genome) {
    message("Translating Alternate Splice Peptides...", appendLF = FALSE)

    AS_Table <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Splice.fst"))
    )
    Proteins_Splice <- as.data.table(
        read.fst(file.path(reference_path, "fst", "Proteins.fst"))
    )
    AS_Table.Extended <- copy(AS_Table)

    Proteins_Splice$exon_number <- as.numeric(Proteins_Splice$exon_number)
    # make phase easier for me to understand
    Proteins_Splice[, c("phase") := -get("phase") %% 3]

    AS_Table.Extended = .gen_splice_proteins_upstream(AS_Table, 
        AS_Table.Extended, Proteins_Splice, genome, isoform = "A")
    AS_Table.Extended = .gen_splice_proteins_upstream(AS_Table, 
        AS_Table.Extended, Proteins_Splice, genome, isoform = "B")
    AS_Table.Extended = .gen_splice_proteins_downstream(AS_Table, 
        AS_Table.Extended, Proteins_Splice, genome, isoform = "A")
    AS_Table.Extended = .gen_splice_proteins_downstream(AS_Table, 
        AS_Table.Extended, Proteins_Splice, genome, isoform = "B")    
    AS_Table.Extended = .gen_splice_proteins_casette(AS_Table, 
        AS_Table.Extended, Proteins_Splice, genome, isoform = "A")
    AS_Table.Extended = .gen_splice_proteins_casette(AS_Table, 
        AS_Table.Extended, Proteins_Splice, genome, isoform = "B")

    AS_Table.Extended[, c("AA_full_A") := ""]
    AS_Table.Extended[!is.na(get("AA_upstr_A")),
        c("AA_full_A") := paste0(get("AA_full_A"), get("AA_upstr_A"))]
    AS_Table.Extended[!is.na(get("AA_casette_A")),
        c("AA_full_A") := paste0(get("AA_full_A"), get("AA_casette_A"))]
    AS_Table.Extended[!is.na(get("AA_downstr_A")),
        c("AA_full_A") := paste0(get("AA_full_A"), get("AA_downstr_A"))]
    AS_Table.Extended[, c("AA_full_B") := ""]
    AS_Table.Extended[!is.na(get("AA_upstr_B")),
        c("AA_full_B") := paste0(get("AA_full_B"), get("AA_upstr_B"))]
    AS_Table.Extended[!is.na(get("AA_casette_B")),
        c("AA_full_B") := paste0(get("AA_full_B"), get("AA_casette_B"))]
    AS_Table.Extended[!is.na(get("AA_downstr_B")),
        c("AA_full_B") := paste0(get("AA_full_B"), get("AA_downstr_B"))]
    write.fst(as.data.frame(AS_Table.Extended),
        file.path(reference_path, "fst", "Splice.Extended.fst"))
    message("done\n")
}

.gen_splice_proteins_upstream <- function(AS_Table, AS_Table.Extended,
        Proteins_Splice, genome, isoform = c("A", "B")) {
    # Upstream applicable for MXE, SE, ALE, A3SS
    cols = c("EventType", "EventID", 
        paste0(c("transcript_id_", "intron_number_"), tolower(isoform)))   
    Upstream <- AS_Table[get("EventType") %in% c("MXE", "SE", "ALE", "A3SS"),
        cols, with = FALSE]
    Upstream[, c("transcript_id", "exon_number") :=
        list(
            get(paste0("transcript_id_", tolower(isoform))), 
            get(paste0("intron_number_", tolower(isoform)))
        )]
    # left_join with Exons
    Upstream <- Proteins_Splice[Upstream,
        on = c("transcript_id", "exon_number"),
        c("EventID", "seqnames", "start", "end", "width", "strand", "phase")
    ]
    Upstream_gr <- .grDT(na.omit(Upstream), keep.extra.columns = TRUE)
    Upstream_seq <- getSeq(genome, Upstream_gr)
    Upstream[!is.na(get("seqnames")), c("seq") := as.character(Upstream_seq)]
    # Trim sequence by phase
    seq <- substr(
        Upstream$seq[!is.na(Upstream$seqnames)],
        1 + (3 - Upstream$phase[!is.na(Upstream$seqnames)]) %% 3,
        nchar(Upstream$seq[!is.na(Upstream$seqnames)])
    )
    # trim last n bases
    seq <- substr(seq, 1, nchar(seq) - (nchar(seq) %% 3))
    # translate
    prot <- Biostrings::translate(as(seq, "DNAStringSet"))
    Upstream[!is.na(get("seqnames")), 
        c("DNA_seq", "AA_seq") := list(
        seq, as.character(prot))
    ]
    cols = c("EventID", "DNA_seq", "AA_seq")
    AS_Table.Extended[Upstream[, cols, with = FALSE], on = "EventID",
        c(paste0(c("DNA_upstr_", "AA_upstr_"), toupper(isoform))) := 
        list(get("DNA_seq"), get("AA_seq"))
    ]
    return(AS_Table.Extended)
}

.gen_splice_proteins_downstream <- function(AS_Table, AS_Table.Extended,
        Proteins_Splice, genome, isoform = c("A", "B")) {
    # Add EventType as exon_number is conditional on this
    cols = c("EventType", "EventID", 
        paste0(c("transcript_id_", "intron_number_"), tolower(isoform)))
    Downstream <- AS_Table[get("EventType") %in% c("MXE", "SE", "AFE", "A5SS"),
        cols, with = FALSE]
    Downstream[, c("transcript_id", "exon_number") :=
        list(get(paste0("transcript_id_", tolower(isoform))), 
            get(paste0("intron_number_", tolower(isoform))))]
    # Modify downstream exon number
    if(toupper(isoform) == "A") {
        Downstream[get("EventType") %in% c("MXE", "SE"),
            c("exon_number") := get("exon_number") + 2]
        Downstream[get("EventType") %in% c("AFE", "A5SS"),
            c("exon_number") := get("exon_number") + 1]    
    } else {
        Downstream[get("EventType") %in% c("MXE"), 
            c("exon_number") := get("exon_number") + 2]
        Downstream[get("EventType") %in% c("SE", "AFE", "A5SS"),
            c("exon_number") := get("exon_number") + 1]    
    }
    Downstream <- Proteins_Splice[Downstream,
        on = c("transcript_id", "exon_number"),
        c("EventID", "seqnames", "start", "end", "width", "strand", "phase")
    ] # left_join with Exons
    Downstream_gr <- .grDT(na.omit(Downstream), keep.extra.columns = TRUE)
    Downstream_seq <- getSeq(genome, Downstream_gr)
    Downstream[!is.na(get("seqnames")), 
        c("seq") := as.character(Downstream_seq)]
    seq <- substr(
        Downstream$seq[!is.na(Downstream$seqnames)],
        1 + (3 - Downstream$phase[!is.na(Downstream$seqnames)]) %% 3,
        nchar(Downstream$seq[!is.na(Downstream$seqnames)])
    ) # Trim sequence by phase
    seq <- substr(seq, 1, nchar(seq) - (nchar(seq) %% 3)) # trim last n bases
    prot <- Biostrings::translate(as(seq, "DNAStringSet"))
    Downstream[!is.na(get("seqnames")), 
        c("DNA_seq", "AA_seq") := list(
        seq, as.character(prot))]
    cols = c("EventID", "DNA_seq", "AA_seq")
    AS_Table.Extended[Downstream[, cols, with = FALSE], on = "EventID",
        c(paste0(c("DNA_downstr_", "AA_downstr_"), toupper(isoform))) := 
        list(get("DNA_seq"), get("AA_seq"))] # translate
    return(AS_Table.Extended)
}

.gen_splice_proteins_casette <- function(AS_Table, AS_Table.Extended,
        Proteins_Splice, genome, isoform = c("A", "B")) {
    cols = c("EventType", "EventID", 
        paste0(c("transcript_id_", "intron_number_"), tolower(isoform)))
    Casette <- AS_Table[, cols, with = FALSE]
    Casette[, c("transcript_id", "exon_number") :=
        list(get(paste0("transcript_id_", tolower(isoform))), 
            get(paste0("intron_number_", tolower(isoform))))]
    if(toupper(isoform) == "A") {
        Casette[get("EventType") %in% c("MXE", "SE", "ALE", "A3SS"),
            c("exon_number") := get("exon_number") + 1]
    } else {
        Casette = Casette[get("EventType") != "SE"]
        Casette[get("EventType") %in% c("MXE", "ALE", "A3SS"),
            c("exon_number") := get("exon_number") + 1]    
    }
    Casette <- Proteins_Splice[Casette,
        on = c("transcript_id", "exon_number"),
        c("EventID", "seqnames", "start", "end", "width", "strand", "phase")]
    Casette_gr <- .grDT(na.omit(Casette), keep.extra.columns = TRUE)
    Casette_seq <- getSeq(genome, Casette_gr)
    Casette[!is.na(get("seqnames")),
        c("casette_seq") := as.character(Casette_seq)]
    setnames(Casette, "phase", "phase_casette")
    # Add nucleotides from upstream and downstream
    upstream = paste0("DNA_upstr_", toupper(isoform))
    downstream = paste0("DNA_downstr_", toupper(isoform))
    cols = c("EventID", upstream, downstream)
    AS_Table.seq = AS_Table.Extended[, cols, with = FALSE]
    cols = c("EventID", "phase_casette", "casette_seq", 
        upstream, downstream)
    Casette <- AS_Table.seq[Casette,
        on = "EventID",
        cols, with = FALSE]

    # Construct extended casette sequence:
    Casette[, c("casette_seq_extended") := get("casette_seq")]
    # Trim casette_seq_extended if upstream sequence does not exists
    Casette[!is.na(get("phase_casette")) & is.na(get(upstream)),
        c("casette_seq_extended") := substr(
            get("casette_seq_extended"),
            get("phase_casette") + 1,
            nchar(get("casette_seq_extended")))]
    Casette[!is.na(get("phase_casette")) & get("phase_casette") > 0 &
            !is.na(get(upstream)),
        c("casette_seq_extended") := paste0(
            substr(get(upstream),
                nchar(get(upstream)) + 1 - get("phase_casette"),
                nchar(get(upstream))
            ),
            get("casette_seq_extended")
        )
    ]
    Casette[nchar(get("casette_seq_extended")) %% 3 > 0 & 
            !is.na(get(downstream)),
        c("casette_seq_extended") := paste0(
            get("casette_seq_extended"),
            substr(get(downstream), 1, 
            3 - (nchar(get("casette_seq_extended")) %% 3))
        )
    ]
    # Translate:
    seq <- Casette$casette_seq_extended[
        !is.na(Casette$casette_seq_extended)]
    # trim out-of-phase to be tidy:
    seq <- substr(seq, 1, nchar(seq) - (nchar(seq) %% 3))
    prot <- Biostrings::translate(as(seq, "DNAStringSet"))
    Casette[!is.na(get("casette_seq_extended")),
        c("DNA_seq", "AA_seq") := list(
        seq, as.character(prot))]
    cols = c("EventID", "DNA_seq", "AA_seq")
    AS_Table.Extended[Casette[, cols, with = FALSE], on = "EventID",
        c(paste0(c("DNA_casette_", "AA_casette_"), toupper(isoform))) := 
        list(get("DNA_seq"), get("AA_seq"))]
    return(AS_Table.Extended)
}

#' @describeIn BuildReference One-step function that fetches resources,
#'   creates a STAR reference (including mappability calculations), then
#'   creates the NxtIRF reference
#' @export
BuildReference_Full <- function(
        reference_path,
        fasta, gtf,
        convert_chromosome_names = NULL,
        overwrite_resource = FALSE,
        genome_type = genome_type,
        nonPolyARef = GetNonPolyARef(genome_type), 
        BlacklistRef = "", 
        UseExtendedTranscripts = TRUE,
        n_threads = 4
) {
    GetReferenceResource(reference_path = reference_path,
        fasta = fasta, gtf = gtf,
        generate_mappability_reads = TRUE,
        convert_chromosome_names = convert_chromosome_names,
        overwrite_resource = overwrite_resource)
    
    STAR_buildRef(reference_path = reference_path, 
        also_generate_mappability = TRUE, 
        n_threads = n_threads)

    BuildReference(reference_path = reference_path,
        genome_type = genome_type,
        nonPolyARef = nonPolyARef, 
        # MappabilityRef = MappabilityRef, 
        BlacklistRef = BlacklistRef, 
        UseExtendedTranscripts = UseExtendedTranscripts)
}

#' Creates a custom Mappability Exclusion BED file
#'
#' @description
#' These function create a custom Mappability Exclusion BED file from the 
#' supplied genome FASTA and gene annotation GTF file. 
#' 
#' @details Creating a Mappability Exclusion BED file is a three-step process.
#' \cr\cr
#' First, using `Mappability_GenReads()`,
#' reads are systematically generated by the given genome FASTA file 
#' (or the `reference_path` containing the genome resource).\cr\cr
#' Second, an aligner
#' such as STAR (preferably the same aligner used for the subsequent RNA-seq 
#' experiment) is required to align these reads to the source genome. Poorly
#' mapped regions of the genome will be reflected by regions of low coverage 
#' depth.\cr\cr
#' Finally, the BAM file containing the aligned reads is analysed 
#' using `Mappability_CalculateExclusions()`, to identify
#' low-mappability regions to compile the Mappability Exclusion BED file.
#' 
#' @param fasta_file The path to the user-supplied genome fasta file
#' @param reference_path The directory to store the reference files
#' @param read_len The nucleotide length of the generated reads
#' @param read_stride The nucleotide distance between sequentially generated
#'   reads
#' @param error_pos The position of the generated error within the reads
#' @param verbose Whether additional status messages are shown
#' @param aligned_bam The BAM file of reads (generated by 
#'   `Mappability_GenReads()`), aligned by the user using their
#'   aligner of choice.
#' @param threshold Regions with this read depth (or below) are defined as low
#'   mappability regions.
#' @param n_threads The number of threads used to calculate mappability
#'   exclusion regions from aligned bam file of synthetic reads.
#' @return None. `Mappability_GenReads()` writes `MappabilityReads.fa`
#'   to the given `reference_path`. `Mappability_CalculateExclusions()` writes
#'   its output BED file as named by `output_file` inside `reference_path`,
#'   appended by the ".txt" suffix. This BED file can be used as input
#'   as `MappabilityRef` in `BuildReference()`
#' @examples
#' 
#' # (1) Creates resource files and systematically generate reads based on
#' # the NxtIRF mock genome:
#' 
#' GetReferenceResource(
#'     reference_path = file.path(tempdir(), "Reference"),
#'     fasta = mock_genome(), gtf = mock_gtf()
#' )
#' Mappability_GenReads(
#'     reference_path = file.path(tempdir(), "Reference"),
#' )
#' 
#' \dontrun{
#' # Setting generate_mappability_reads = TRUE to GetReferenceResource() will
#' # run Mappability_GenReads after resource retrieval. The following code
#' is equivalent to that of above:
#' 
#' GetReferenceResource(
#'     reference_path = file.path(tempdir(), "Reference"),
#'     fasta = mock_genome(), gtf = mock_gtf(),
#'     generate_mappability_reads = TRUE
#' )
#' 
#' # (2) Align the generated reads using Rsubread:
#' 
#' setwd(file.path(tempdir(), "Reference"))
#' Rsubread::buildindex(basename = "./reference_index", 
#'     reference = mock_genome())
#' Rsubread::subjunc(
#'     index = "./reference_index", 
#'     readfile1 = "MappabilityReads.fa", 
#'     output_file = "MappabilityReads.bam", 
#'     useAnnotation = TRUE, 
#'     annot.ext = mock_gtf(), 
#'     isGTF = TRUE
#' )
#' 
#' # (3) Analyse the aligned reads for low-mappability regions:
#' 
#' Mappability_CalculateExclusions(aligned_bam = "MappabilityReads.bam",
#'     output_file = file.path(tempdir(), "Reference", "Mappability.bed")
#' )
#' }
#' @name Mappability-methods
#' @aliases 
#' Mappability_GenReads
#' Mappability_CalculateExclusions
#' @seealso [BuildReference]
#' <https://github.com/williamritchie/IRFinder/blob/master/bin/util/generateReadsError.pl>
#' @md
NULL

#' @describeIn Mappability-methods Generates Mappability reads from a 
#' genome FASTA file. This function replicates the functionality of 
#' generateReadsError.pl in vanilla IRFinder.
#' @export
Mappability_GenReads <- function(reference_path, fasta_file,
        read_len = 70, read_stride = 10, error_pos = 35,
        verbose = TRUE) {
    .gmr_check_params(read_len, read_stride, error_pos)
    if(missing(fasta_file)) {
        fasta_file = .STAR_get_FASTA(reference_path)
        if(!file.exists(fasta_file)) {
            .log(paste("In Mappability_GenReads,",
                "failed to generate genome fasta file from given reference"))
        }
    } else if(!file.exists(fasta_file)) {
        .log(paste("In Mappability_GenReads,",
            "given fasta file", fasta_file, "not found"))
    }
    .validate_path(file.path(normalizePath(reference_path), "Mappability"))
    # Run map read generator:
    run_IRFinder_GenerateMapReads(
        normalizePath(fasta_file),
        file.path(normalizePath(reference_path), "Mappability", "Reads.fa"),
        read_len, read_stride, error_pos
    )
}

#' @describeIn Mappability-methods Generate a BED file defining 
#' low mappability regions, using reads generated by 
#' \code{Mappability_GenReads()}, aligned to the genome.
#' @export
Mappability_CalculateExclusions <- function(reference_path, 
        aligned_bam = file.path(reference_path, "Mappability", 
            "Aligned.out.bam"), 
        threshold = 4, n_threads = 1) {
    if(!file.exists(aligned_bam)) {
        .log(paste("In Mappability_CalculateExclusions(),",
            aligned_bam, "BAM file does not exist"))
    }

    .validate_path(file.path(normalizePath(reference_path), "Mappability"))
    output_file = file.path(normalizePath(reference_path), "Mappability",
        "MappabilityExclusion.bed")
        
    .log(paste("Calculating Mappability Exclusion regions from:",
        aligned_bam), type = "message")
    run_IRFinder_MapExclusionRegions(
        bamfile = normalizePath(aligned_bam),
        output_file = output_file,
        threshold = threshold,
        n_threads = n_threads
    )
}

.gmr_check_params <- function(read_len, read_stride, error_pos) {
    if(!is.numeric(read_len) || read_len < 30) {
        .log(paste("In Mappability_GenReads,",
            "read_len must be numerical and at least 30"))
    }
    if(!is.numeric(read_stride) || read_stride > read_len) {
        .log(paste("In Mappability_GenReads,",
            "read_stride must be numerical and less than read_len"))
    }
    if(!is.numeric(error_pos) || error_pos > read_len) {
        .log(paste("In Mappability_GenReads,",
            "error_pos must be numerical and less than read_len"))
    }
}

.genmapreads_validate <- function(fasta, reference_path) {
    ah_genome = ""
    if (!dir.exists(dirname(reference_path))) {
        .log(paste("In Mappability_GenReads",
            dirname(reference_path), "does not exist"))
    }
    if(file.exists(fasta)) {
        if({
            tryCatch({
                genome = rtracklayer::import(fasta, "fasta")
                FALSE
            }, error = function(e) TRUE)
        }) {
            .log(paste("In Mappability_GenReads",
                fasta, "does not appear to be a valid FASTA file"))
        }
    } else {
        if (!file.exists(file.path(reference_path, "settings.Rds"))) {
            .log(paste("In Mappability_GenReads",
                "settings.Rds is not found.",
                "Run GetReferenceResource() first or supply a FASTA file"))
        }
        settings.Rds = readRDS(file.path(reference_path, "settings.Rds"))
        ah_genome = settings.Rds[["ah_genome"]]
        if (ah_genome != "" && substr(ah_genome, 1, 2) != "AH") {
            .log(paste("In Mappability_GenReads",
                ah_genome, "is not a valid AnnotationHub record name"))
        } else if(ah_genome == "" && !file.exists(
                file.path(reference_path, "resource", "genome.2bit"))) {
            .log(paste("In Mappability_GenReads",
                "resource/genome.2bit", "was not found in", reference_path,
                ". Perhaps run GetReferenceResource() again or supply",
                "a FASTA file"))
        }
    }
    if (!dir.exists(file.path(reference_path))) {
        dir.create(file.path(reference_path))
    }
    if (!dir.exists(file.path(reference_path, "resource"))) {
        dir.create(file.path(reference_path, "resource"))
    }
    return(ah_genome)
}




