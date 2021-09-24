#' NxtIRF Example BAMs and NxtSE Experiment Object
#'
#' 6 example bam files based on samples from the 
#' Leucegene dataset (GSE67039). Bam files are constructed
#' based on the complete bam files of 6 samples from Leucegene,
#' subsetted by regions containing the 7 above genes. Then, the reads of these 
#' subsetted BAMs were realigned to the NxtIRF example reference using 
#' STAR.\cr\cr
#' A NxtSE SummarizedExperiment object was constructed by running `IRFinder()`,
#' `CollateData()` and `MakeSE()`
#' on the example BAM files, using the NxtIRF example reference. A full pipeline
#' of how this NxtSE was generated can be found in
#' `system.file("scripts", "make_data.R", package = "NxtIRFcore")`
#' @return See `NxtIRF_example_bams()` returns a 2-column data frame containing
#'   sample names and BAM paths of the example dataset. `NxtIRF_example_NxtSE()`
#'   returns a NxtSE object.
#' @examples
#'
#' # returns a data frame with the first column as sample names, and the 
#' # second column as BAM paths
#'
#' NxtIRF_example_bams() 
#'
#' # Returns a NxtSE object created by the example bams aligned to the
#' # mock NxtSE reference
#'
#' se = NxtIRF_example_NxtSE() 
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
#' NxtIRF_example_bams NxtIRF_example_NxtSE
#' @keywords package
#' @seealso [BuildReference()], [IRFinder()], [CollateData()], [MakeSE()]
#' @md
NULL

#' @describeIn example-NxtIRF-data Returns a 2-column data frame, containing 
#'   sample names and sample paths (in tempdir()) of example BAM files
#' @export
NxtIRF_example_bams <- function() {
    example_bams()
    return(Find_Bams(tempdir()))
}

#' @describeIn example-NxtIRF-data Returns a (in-memory / realized) NxtSE object 
#'   generated using the NxtIRF mock reference and example BAM files
#' @export
NxtIRF_example_NxtSE <- function() {
    se = readRDS(system.file("extdata", 
        "example_NxtSE.Rds", package = "NxtIRFcore"))
    covs = Find_Samples(system.file("extdata", package = "NxtIRFcore"), ".cov")
    covfile(se) <- covs$path
    se
}
