#' @useDynLib NxtIRF, .registration = TRUE 
#' @import NxtIRFdata
#' @import data.table
#' @importFrom fst read.fst write.fst
#' @importFrom DBI dbConnect dbDisconnect
#' @importFrom RSQLite SQLite dbWriteTable dbReadTable dbListTables
#' @import shiny
#' @import shinydashboard
#' @import shinyFiles
#' @import shinyWidgets
#' @import rhandsontable
#' @import ggplot2
#' @importFrom boot logit inv.logit
#' @importFrom methods as is coerce callNextMethod
#' @importFrom graphics text
#' @importFrom stats as.formula model.matrix qt runif na.omit prcomp
#' @importFrom utils download.file packageVersion getFromNamespace
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Rcpp evalCpp
#' @importFrom AnnotationHub AnnotationHub cache
#' @importFrom BiocFileCache BiocFileCache bfcrpath
#' @importFrom BiocGenerics start end width strand
#' @importFrom Biostrings getSeq readDNAStringSet DNAStringSet translate replaceAmbiguities type
#' @importFrom BiocParallel bpparam bplapply SnowParam MulticoreParam SerialParam
#' @importFrom DT datatable selectRows 
#' @importFrom genefilter rowttests
#' @importFrom grDevices colorRampPalette
#' @importFrom GenomeInfoDb sortSeqlevels seqinfo seqlengths seqlevels<- seqlevels
#' @importFrom GenomicRanges GRanges reduce findOverlaps 
#' @importFrom GenomicRanges makeGRangesFromDataFrame 
#' @importFrom GenomicRanges makeGRangesListFromDataFrame mcols split strand 
#' @importFrom GenomicRanges flank setdiff seqnames psetdiff disjoin mcols<- 
#' @importFrom GenomicRanges strand<- seqnames<-
#' @importFrom heatmaply heatmaply
#' @importFrom IRanges IRanges Views RleList
#' @importFrom matrixStats rowSds
#' @importFrom openssl md5
#' @importFrom parallel detectCores
#' @importFrom plotly config layout plotlyOutput event_data ggplotly plotlyProxy plotlyProxyInvoke renderPlotly subplot highlight
#' @importFrom R.utils gzip
#' @importFrom rappdirs user_cache_dir
#' @importFrom rtracklayer import export TwoBitFile track
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom stringr str_locate
#' @importFrom SummarizedExperiment SummarizedExperiment 
#' @importFrom SummarizedExperiment rowData colData rowData<- colData<-
#' @importFrom SummarizedExperiment assay assays assay<- assays<-
#' @importFrom SummarizedExperiment assayNames assayNames<-
#' @importFrom S4Vectors coolcat metadata Rle metadata<- SimpleList 
#' @importFrom S4Vectors endoapply from to setValidity2
#' @importFrom XML getHTMLLinks
NULL

S4_disableValidity <- getFromNamespace("disableValidity", "S4Vectors")
BG_replaceSlots <- getFromNamespace("replaceSlots", "BiocGenerics")
