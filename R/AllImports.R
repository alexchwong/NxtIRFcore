#' @useDynLib NxtIRF, .registration = TRUE 
#' @importFrom magrittr %>%
#' @import data.table
#' @import shiny
#' @import shinydashboard
#' @import ggplot2
#' @importFrom shinyFiles getVolumes parseDirPath parseFilePaths parseSavePath
#' @importFrom shinyFiles shinyDirButton shinyDirChoose shinyFileChoose
#' @importFrom shinyFiles shinyFilesButton shinyFileSave shinySaveButton
#' @importFrom shinyWidgets sliderTextInput updateSliderTextInput
#' @importFrom shinyWidgets radioGroupButtons updateRadioGroupButtons
#' @importFrom shinyWidgets switchInput actionBttn
#' @importFrom shinyWidgets sendSweetAlert ask_confirmation
#' @importFrom rhandsontable rhandsontable hot_to_r
#' @importFrom rhandsontable renderRHandsontable rHandsontableOutput
#' @importFrom tools R_user_dir
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
#' @importFrom BiocParallel SnowParam MulticoreParam SerialParam
#' @importFrom BiocParallel bpparam bplapply
#' @importFrom DelayedArray qlogis plogis rowMeans DelayedArray
#' @importFrom DelayedMatrixStats rowSds colVars
#' @importFrom DT datatable selectRows 
#' @importFrom ExperimentHub ExperimentHub
#' @importFrom fst read.fst write.fst
#' @importFrom genefilter rowttests
#' @importFrom grDevices colorRampPalette
#' @importFrom GenomeInfoDb sortSeqlevels seqinfo seqlengths seqlevels<- 
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges GRanges reduce findOverlaps 
#' @importFrom GenomicRanges makeGRangesFromDataFrame 
#' @importFrom GenomicRanges makeGRangesListFromDataFrame mcols split strand 
#' @importFrom GenomicRanges flank setdiff seqnames psetdiff disjoin mcols<- 
#' @importFrom GenomicRanges strand<- seqnames<-
#' @importFrom HDF5Array HDF5Array writeHDF5Array h5writeDimnames 
#' @importFrom heatmaply heatmaply
#' @importFrom httr HEAD GET status_code
#' @importFrom IRanges IRanges Views RleList
#' @importFrom openssl md5
#' @importFrom parallel detectCores
#' @importFrom plotly config layout plotlyOutput event_data ggplotly 
#' @importFrom plotly plotlyProxy plotlyProxyInvoke renderPlotly subplot 
#' @importFrom plotly highlight
#' @importFrom R.utils gzip
#' @importFrom rtracklayer import export TwoBitFile track
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom rhdf5 h5createFile h5createDataset h5delete h5write h5createGroup
#' @importFrom stringr str_locate
#' @importFrom SummarizedExperiment SummarizedExperiment 
#' @importFrom SummarizedExperiment rowData colData rowData<- colData<-
#' @importFrom SummarizedExperiment assay assays assay<- assays<-
#' @importFrom SummarizedExperiment assayNames assayNames<-
#' @importFrom S4Vectors coolcat metadata Rle metadata<- SimpleList 
#' @importFrom S4Vectors endoapply from to setValidity2 DataFrame
#' @importFrom S4Vectors bindCOLS bindROWS getListElement setListElement
#' @importFrom XML getHTMLLinks
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom S4Vectors DataFrame RleList
NULL

S4_disableValidity <- getFromNamespace("disableValidity", "S4Vectors")
BG_replaceSlots <- getFromNamespace("replaceSlots", "BiocGenerics")
SE_charbound <- 
    function(idx, txt, fmt)
{
    orig <- idx
    idx <- match(idx, txt)
    if (any(bad <- is.na(idx))) {
        msg <- paste(S4_selectSome(orig[bad]), collapse=" ")
        stop(sprintf(fmt, msg))
    }
    idx
}

S4_selectSome <- getFromNamespace("selectSome", "S4Vectors")
