df <- data.frame(
    Title = sprintf("Mappability Exclusion Regions for %s", 
        c("hg38 Ens v94", "hg19 Ens v75", "mm10 Ens v94", "mm9 Ens v67")), 
    Description = paste(c("hg38 Ens v94", "hg19 Ens v75", "mm10 Ens v94", "mm9 Ens v67"),
        "mappability exclusion regions were generated using a modified script of the method ",
        "as described in https://github.com/williamritchie/IRFinder/blob/master/bin/util/Mapability"
    ),
    RDataPath = file.path("NxtIRF", "mappability", 
        c("hg38_v94.rds", "hg19_v75.rds", "mm10_v94.rds", "mm9_v67.rds")),
    BiocVersion="3.12", 
    Genome=c("hg38", "hg19", "mm10", "mm9"), 
    SourceType="BED", 
    SourceUrl=paste(
        "https://raw.github.com/alexw-gsct/NxtIRF_resources/main/data",
        c(
            "Mappability_Regions_hg38_v94.txt.gz",
            "Mappability_Regions_hg19_v75.txt.gz",
            "Mappability_Regions_mm10_v94.txt.gz",
            "Mappability_Regions_mm9_v67.txt.gz"
        )
    ),
    SourceVersion="0.98.0",
    Species=c(rep("Homo sapiens", 2),rep("Mus musculus", 2)),
    TaxonomyId=c(rep("9606", 2),rep("10090", 2)),
    Coordinate_1_based=FALSE,
    DataProvider="Alex Wong",
    Maintainer="Alex Wong <a.wong@centenary.org.au>",
    RDataClass="GRanges",
    DispatchClass="BEDFile",
    Location_Prefix = "https://raw.github.com/alexw-gsct/NxtIRF_resources/main/data",
    RDataPath = c(
        "Mappability_Regions_hg38_v94.txt.gz",
        "Mappability_Regions_hg19_v75.txt.gz",
        "Mappability_Regions_mm10_v94.txt.gz",
        "Mappability_Regions_mm9_v67.txt.gz"
    ),
    Tags = "MappabilityExclusion:Mappability:NxtIRF", 
    stringsAsFactors = FALSE
)

write.csv(file="../extdata/metadata.csv", df, row.names=FALSE)
