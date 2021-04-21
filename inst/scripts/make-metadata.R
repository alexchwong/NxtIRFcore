bam_samples <- c("02H003", "02H025", "02H026", "02H033", "02H043", "02H046")

source_repo <- "https://raw.github.com/alexchwong/NxtIRFdata/main/inst/"

genes <- "SRSF1, SRSF2, SRSF3, TRA2A, TRA2B, TP53 and NSUN5"

df_refs <- data.frame(
    Title = sprintf(paste("NxtIRF mock %s containing human genes:", genes),
        c("genome", "gene annotations")), 
    Description = paste(
        "Sequences from genes (", genes, ") were pasted together to create",
        "an artificial 'chrZ' chromosome. Chromosomal coordinates for gene",
        "annotations were modified accordingly."
    ),
    BiocVersion="3.13", 
    Genome="NxtIRF_chrZ", 
    SourceType=c("FASTA", "GTF"),
    SourceUrl=
        sprintf(paste0(source_repo, "NxtIRF/%s"),
            c("genome.fa", "transcripts.gtf")
        ),
    SourceVersion="0.99.0",
    Species="Homo sapiens",
    TaxonomyId="9606",
    Coordinate_1_based=1,
    DataProvider="NxtIRF",
    Maintainer="Alex Wong <a.wong@centenary.org.au>",
    RDataClass=c("FaFile", "GRanges"),
    DispatchClass=c("FaFile", "GFFFile"),
    Location_Prefix = source_repo,
    RDataPath = sprintf("NxtIRF/%s", c("genome.fa", "transcripts.gtf")),
    Tags = "Annotation:MockGenome:NxtIRF",
    stringsAsFactors = FALSE
)

df_bams <- data.frame(
    Title = sprintf(paste("RNA-seq reads aligned to the NxtIRF mock genome",
        "from %s of the Leucegene dataset (GSE67039)"), 
        c("02H003", "02H025", "02H026", "02H033", "02H043", "02H046")), 
    Description = paste(
        c("02H003", "02H025", "02H026", "02H033", "02H043", "02H046"),
        "- aligned reads (from this sample in GSE67039) were filtered by",
        "those mapping to genes in the NxtIRF mock genome.",
        "These were re-aligned to the mock genome using STAR."
    ),
    BiocVersion="3.13", 
    Genome="NxtIRF_chrZ", 
    SourceType="BAM",
    SourceUrl=
        sprintf(
            "%s.bam",
            paste0(source_repo, "NxtIRF/", bam_samples)
        ),
    SourceVersion="0.99.0",
    Species="Homo sapiens",
    TaxonomyId="9606",
    Coordinate_1_based=1,
    DataProvider="NxtIRF",
    Maintainer="Alex Wong <a.wong@centenary.org.au>",
    RDataClass="BamFile",
    DispatchClass="BamFile",
    Location_Prefix = source_repo,
    RDataPath = sprintf(
            paste("%s.bam", "%s.bam.bai", sep = ", "),
            paste0("NxtIRF/", bam_samples),
            paste0("NxtIRF/", bam_samples)
        ),
    Tags = "ExperimentData:MockGenome:Leucegene:NxtIRF",
    stringsAsFactors = FALSE
)

df_mappa <- data.frame(
    Title = sprintf("Mappability Exclusion Regions for %s", 
        c(
            "Ensembl GRCh38 (hg38) release-94", 
            "Ensembl GRCh37 (hg19) release-75",
            "Ensembl GRCm38 (mm10) release-94",
            "Ensembl NCBIM37 (mm9) release-67"
        )
    ),
    Description = paste(
        c(
            "Ensembl GRCh38 (hg38) release-94", 
            "Ensembl GRCh37 (hg19) release-75",
            "Ensembl GRCm38 (mm10) release-94",
            "Ensembl NCBIM37 (mm9) release-67"
        ),
        "mappability exclusion regions were generated using a modified script",
        "of the method as described in",
        paste0("https://github.com/williamritchie/IRFinder/",
            "blob/master/bin/util/Mapability")
    ),
    BiocVersion="3.13", 
    Genome=c("hg38", "hg19", "mm10", "mm9"),
    SourceType="BED",
    SourceUrl=paste0(
        source_repo, "NxtIRF/",
        c(
            "Mappability_Regions_hg38_v94.txt.gz",
            "Mappability_Regions_hg19_v75.txt.gz",
            "Mappability_Regions_mm10_v94.txt.gz",
            "Mappability_Regions_mm9_v67.txt.gz"
        )
    ),
    SourceVersion="0.99.0",
    Species=c(rep("Homo sapiens", 2),rep("Mus musculus", 2)),
    TaxonomyId=c(rep("9606", 2),rep("10090", 2)),
    Coordinate_1_based=0,
    DataProvider="NxtIRF",
    Maintainer="Alex Wong <a.wong@centenary.org.au>",
    RDataClass="GRanges",
    DispatchClass="BEDFile",
    Location_Prefix = source_repo,
    RDataPath = c(
        "NxtIRF/Mappability_Regions_hg38_v94.txt.gz",
        "NxtIRF/Mappability_Regions_hg19_v75.txt.gz",
        "NxtIRF/Mappability_Regions_mm10_v94.txt.gz",
        "NxtIRF/Mappability_Regions_mm9_v67.txt.gz"
    ),
    Tags = "Annotation:MappabilityExclusion:Mappability:NxtIRF", 
    stringsAsFactors = FALSE
)

df = rbind(df_refs, df_bams, df_mappa)

write.csv(file="../extdata/metadata.csv", df, row.names=FALSE)
