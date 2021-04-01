#' NxtIRF Examples
#'
#' This package contains files that provides a workable example for the 
#' NxtIRF package.\cr\cr
#' A mock reference, with genome sequence (FASTA) and gene annotation (GTF)
#' files are provided, based on the genes SRSF1, SRSF2, SRSF3, TRA2A, TRA2B, 
#' TP53 and NSUN5, of which sequences are used to construct an artificial 
#' chromosome Z. This was generated based on release-94 of Ensembl GRCh38 (hg38)
#' reference.\cr\cr
#' NxtIRFdata contains 6 example bam files based on samples from the 
#' Leucegene dataset (GSE67039). Bam files are constructed
#' based on the complete bam files of 6 samples from Leucegene,
#' subsetted by regions containing the 7 above genes. Then, the reads of these 
#' subsetted BAMs were realigned to the mock reference using STAR.\cr\cr
#' Additionally, NxtIRFdata contains Mappability exclusion regions generated
#' using NxtIRF, suitable for use in generating references based on hg38,
#' hg19, mm10 and mm9 genomes.
#' @param genome_type Either one of `hg38`, `hg19`, `mm10` or `mm9`
#' @return See Examples section below.
#' @examples
#' mock_genome() # returns the location of the genome.fa file of the mock reference
#'
#' mock_gtf() # returns the location of the transcripts.gtf file of the mock reference
#'
#' example_bams() # returns the locations of the 6 example bam files
#'
#' get_mappability_exclusion("hg38") # returns the location of the Mappability exclusion BED for hg38
#' @references
#' Generation of the mappability files was performed using NxtIRF using
#' a method analogous to that described in:
#' 
#' Middleton R, Gao D, Thomas A, Singh B, Au A, Wong JJ, Bomane A, Cosson B, Eyras E, Rasko JE, Ritchie W.
#' IRFinder: assessing the impact of intron retention on mammalian gene expression.
#' Genome Biol. 2017 Mar 15;18(1):51.
#' \url{https://doi.org/10.1186/s13059-017-1184-4}
#' @name NxtIRF-example-data
#' @aliases 
#' mock_genome
#' mock_gtf
#' example_bams
#' get_mappability_exclusion
#' @keywords package
#' @md
NULL

#' @export
mock_genome <- function()
{
    system.file("extdata", "genome.fa",
        package="NxtIRFdata", mustWork=TRUE)
}

#' @export
mock_gtf <- function()
{
    system.file("extdata", "transcripts.gtf",
        package="NxtIRFdata", mustWork=TRUE)
}

#' @export
example_bams <- function()
{
    c(
        system.file("extdata", "02H003_chrZ.bam",
            package="NxtIRFdata", mustWork=TRUE),
        system.file("extdata", "02H025_chrZ.bam",
            package="NxtIRFdata", mustWork=TRUE),    
        system.file("extdata", "02H026_chrZ.bam",
            package="NxtIRFdata", mustWork=TRUE),    
        system.file("extdata", "02H033_chrZ.bam",
            package="NxtIRFdata", mustWork=TRUE),    
        system.file("extdata", "02H043_chrZ.bam",
            package="NxtIRFdata", mustWork=TRUE),    
        system.file("extdata", "02H046_chrZ.bam",
            package="NxtIRFdata", mustWork=TRUE)  
    )
}

#' @export
get_mappability_exclusion <- function(
        genome_type = c("hg38", "hg19", "mm10", "mm9")) {
    genome_type = match.arg(genome_type)
    if(genome_type == "hg38") {
        system.file("extdata", "Mappability_Regions_hg38_v94.txt.gz",
            package="NxtIRFdata", mustWork=TRUE)    
    } else if(genome_type == "hg19") {
        system.file("extdata", "Mappability_Regions_hg19_v75.txt.gz",
            package="NxtIRFdata", mustWork=TRUE)        
    } else if(genome_type == "mm10") {
        system.file("extdata", "Mappability_Regions_mm10_v94.txt.gz",
            package="NxtIRFdata", mustWork=TRUE)        
    } else if(genome_type == "mm9") {
        system.file("extdata", "Mappability_Regions_mm9_v67.txt.gz",
            package="NxtIRFdata", mustWork=TRUE)        
    } else {
        stop(paste("In get_mappability_exclusion():",
            "genome_type = ", genome_type, "is not recogised"
        ), call. = FALSE)
    }
}

#' @importFrom ExperimentHub ExperimentHub
get_eh_data <- function(id, type,
        localHub = FALSE, eh = ExperimentHub(localHub = localHub)) {

}