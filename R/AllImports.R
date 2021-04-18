#' @useDynLib NxtIRF, .registration = TRUE 
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
#' @importFrom methods as is coerce callNextMethod new
#' @importFrom graphics text
#' @importFrom stats as.formula model.matrix qt runif na.omit prcomp
#' @importFrom utils download.file packageVersion getFromNamespace
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Rcpp evalCpp
#' @importFrom AnnotationHub AnnotationHub cache
#' @importFrom BiocFileCache BiocFileCache bfcrpath bfcquery
#' @importFrom BiocGenerics start end width strand
#' @importFrom BiocGenerics nrow ncol rbind cbind
#' @importFrom Biostrings getSeq readDNAStringSet DNAStringSet translate 
#' @importFrom Biostrings replaceAmbiguities type
#' @importFrom BiocParallel bpparam bplapply SnowParam MulticoreParam 
#' @importFrom BiocParallel SerialParam
#' @importFrom DelayedArray qlogis plogis rowMeans DelayedArray
#' @importFrom DelayedMatrixStats rowSds colVars
#' @importFrom DT datatable selectRows 
#' @importFrom genefilter rowttests
#' @importFrom grDevices colorRampPalette
#' @importFrom GenomeInfoDb sortSeqlevels seqinfo seqlengths seqlevels<- 
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges GRanges reduce findOverlaps 
#' @importFrom GenomicRanges makeGRangesFromDataFrame 
#' @importFrom GenomicRanges makeGRangesListFromDataFrame mcols split strand 
#' @importFrom GenomicRanges flank setdiff seqnames psetdiff disjoin mcols<- 
#' @importFrom GenomicRanges strand<- seqnames<-
#' @importFrom HDF5Array writeHDF5Array loadHDF5SummarizedExperiment
#' @importFrom HDF5Array HDF5Array saveHDF5SummarizedExperiment
#' @importFrom heatmaply heatmaply
#' @importFrom httr HEAD GET status_code
#' @importFrom IRanges IRanges Views RleList
#' @importFrom openssl md5
#' @importFrom parallel detectCores
#' @importFrom plotly config layout plotlyOutput event_data ggplotly 
#' @importFrom plotly plotlyProxy plotlyProxyInvoke renderPlotly subplot 
#' @importFrom plotly highlight
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
#' @importFrom S4Vectors endoapply from to setValidity2 DataFrame
#' @importFrom S4Vectors bindCOLS bindROWS
#' @importFrom XML getHTMLLinks
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom S4Vectors DataFrame
NULL

S4_disableValidity <- getFromNamespace("disableValidity", "S4Vectors")
BG_replaceSlots <- getFromNamespace("replaceSlots", "BiocGenerics")
