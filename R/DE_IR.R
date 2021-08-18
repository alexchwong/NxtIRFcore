.ASE_check_args <- function(colData, test_factor, 
        test_nom, test_denom, batch1, batch2, n_threads) {
    if(!is_valid(test_factor) | !is_valid(test_nom) | !is_valid(test_denom)) {
        .log("test_factor, test_nom, test_denom must be defined")
    }
    if(!(test_factor %in% colnames(colData))) {
        .log("test_factor is not a condition in colData")
    }
    if(!any(colData[, test_factor] == test_nom)) {
        .log("test_nom is not found in any samples")
    }
    if(!any(colData[, test_factor] == test_denom)) {
        .log("test_denom is not found in any samples")
    }
    if(batch1 != "") {
        if(!(batch1 %in% colnames(colData))) {
            .log("batch1 is not a condition in colData")
        }
        if(test_factor == batch1) {
            .log("batch1 and test_factor are the same")
        }
    }
    if(batch2 != "") {
        if(!(batch2 %in% colnames(colData))) {
            .log("batch2 is not a condition in colData")
        }
        if(test_factor == batch2) {
            .log("batch2 and test_factor are the same")
        }
    }
    if(batch1 != "" & batch2 != "") {
        if(batch1 == batch2) {
            .log("batch1 and batch2 are the same")
        }
    }
    return(TRUE)
}



#' Use Limma, DESeq2 or DoubleExpSeq to test for differential ASEs 
#'   (Alternative Splice Events)
#'
#' @param se The SummarizedExperiment object created by `MakeSE()`. To reduce
#'   runtime and false negatives due to multiple testing issues, please filter
#'   the object using `apply_filter()`
#' @param test_factor A string for the column name which contains the 
#'   contrasting variable
#' @param test_nom The condition in which to test for differential ASE. Usually
#'   the "treatment" condition
#' @param test_denom The condition in which to test against for differential 
#'   ASE. Usually the "control" condition
#' @param batch1,batch2 (Limma and DESeq2 only) One or two columns containing 
#'   batch information to normalise against (can be omitted).
#' @param filter_antiover Whether to filter out IR events that overlap 
#'   antisense genes (for unstranded RNAseq protocols)
#' @param filter_antinear Whether to filter out IR events near but not 
#'   overlapping antisense genes (for unstranded RNAseq protocols)
#' @param filter_annotated_IR Whether to filter out IR events that are already 
#'   annotated exons (after doing so, all IR events will be unannotated - 
#'   i.e. constitutionally spliced introns))
#' @param n_threads (DESeq_ASE only) How many threads to use for DESeq2
#'   based analysis.
#' @return A data table containing the following:
#'   EventName: The name of the ASE event\cr\cr
#'   EventType: The type of event. IR = intron retention, MXE = mutually 
#'     exclusive event, SE = skipped exon, AFE = alternate first exon, ALE = 
#'     alternate last exon, A5SS / A3SS = alternate 5' / 3' splice site\cr\cr
#'   EventRegion: The genomic coordinates the event occupies.\cr\cr  
#'   NMD_direction: Indicates whether one isoform is a NMD substrate. +1 means 
#'     included isoform is NMD, -1 means the excluded isoform is NMD, and 0 
#'     means neither (or both) are NMD substrates\cr\cr
#'   AvgPSI_nom, Avg_PSI_denom: the average percent spliced in / percent intron
#'     retention levels for the two conditions being contrasted\cr\cr
#'   (LIMMA SPECIFIC OUTPUT)\cr\cr
#'   logFC, AveExpr, t, P.Value, adj.P.Val, B: limma topTable columns of 
#'     limma results. See `?limma::topTable`\cr\cr
#'   inc/exc_(logFC, AveExpr, t, P.Value, adj.P.Val, B): limma results 
#'     for differential testing for raw included / excluded counts only\cr\cr
#'   (DESEQ2 SPECIFIC OUTPUT)\cr\cr
#'   baseMean, log2FoldChange, lfcSE, stat, pvalue, padj: 
#'     DESeq2 results columns See `?DESeq2::results`\cr\cr
#'   inc/exc_(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj): 
#'     DESeq2 results for differential testing for
#'     raw included / excluded counts only\cr\cr
#'   (DoubleExp SPECIFIC OUTPUT)\cr\cr
#'     MLE: Expectation values for the two groups\cr\cr
#'     MLE_LFC: Log2-fold change of the MLE\cr\cr
#'     P.Value, adj.P.Val: Nominal and BH-adjusted P values\cr\cr
#'     n_eff: Number of effective samples (i.e. non-zero or non-unity PSI)\cr\cr
#'     mDepth: Mean Depth of splice coverage\cr\cr
#'     Dispersion_Reduced, Dispersion_Full: Dispersion values for reduced and 
#'     full models
#' @examples
#' # see ?MakeSE on example code of generating this NxtSE object
#' se = NxtIRF_example_NxtSE()
#' 
#' colData(se)$treatment = rep(c("A", "B"), each = 3)
#' 
#' if("limma" %in% rownames(installed.packages())) {
#'     res_limma = limma_ASE(se, "treatment", "A", "B")
#' }
#'
#' if("DESeq2" %in% rownames(installed.packages())) {
#'     res_DESeq = DESeq_ASE(se, "treatment", "A", "B")
#' }
#' 
#' if("DoubleExpSeq" %in% rownames(installed.packages())) {
#'     res_DES = DoubleExpSeq_ASE(se, "treatment", "A", "B")
#' }
#' @name ASE-methods
#' @md
NULL

#' @describeIn ASE-methods Use limma to perform differential ASE analysis of
#'   a filtered NxtSE object
#' @export
limma_ASE <- function(se, test_factor, test_nom, test_denom, 
        batch1 = "", batch2 = "",
        filter_antiover = TRUE, filter_antinear = FALSE, 
        filter_annotated_IR = FALSE) {
    
    NxtIRF.CheckPackageInstalled("limma", "3.44.0")
    .ASE_check_args(colData(se), test_factor, 
        test_nom, test_denom, batch1, batch2)
    se_use <- .ASE_filter(
        se, filter_antiover, filter_antinear, filter_annotated_IR)
    # se_use <- realize_NxtSE(se_use)
    
    res.limma2 <- .ASE_limma_contrast(se_use, 
        test_factor, test_nom, test_denom,
        batch1, batch2)
    res.inc = res.limma2[grepl(".Included", get("EventName"))]
    res.inc[, 
        c("EventName") := sub(".Included","",get("EventName"), fixed=TRUE)]
    res.inc = res.inc[get("AveExpr") > 1]
    res.exc = res.limma2[grepl(".Excluded", get("EventName"))]
    res.exc[, 
        c("EventName") := sub(".Excluded","",get("EventName"), fixed=TRUE)]
    res.exc = res.exc[get("AveExpr") > 1]

    rowData = as.data.frame(rowData(se_use))
    se_use = se_use[rowData$EventName %in% res.inc$EventName &
        rowData$EventName %in% res.exc$EventName,]
    res.ASE <- .ASE_limma_contrast_ASE(se_use, 
        test_factor, test_nom, test_denom,
        batch1, batch2)
    res.ASE[res.inc, on = "EventName",
        paste("Inc", colnames(res.inc)[seq_len(6)], sep=".") := 
        list(get("i.logFC"), get("i.AveExpr"), get("i.t"), 
            get("i.P.Value"), get("i.adj.P.Val"), get("i.B"))]
    res.ASE[res.exc, on = "EventName",
        paste("Exc", colnames(res.exc)[seq_len(6)], sep=".") := 
        list(get("i.logFC"), get("i.AveExpr"), get("i.t"), 
            get("i.P.Value"), get("i.adj.P.Val"), get("i.B"))]
    setorderv(res.ASE, "B", order = -1)
    res.ASE <- .ASE_add_diag(res.ASE, se_use, test_factor, test_nom, test_denom)
    return(res.ASE)
}

#' @describeIn ASE-methods Use DESeq2 to perform differential ASE analysis of
#'   a filtered NxtSE object
#' @export
DESeq_ASE <- function(se, test_factor, test_nom, test_denom, 
        batch1 = "", batch2 = "",
        n_threads = 1,
        filter_antiover = TRUE, filter_antinear = FALSE, 
        filter_annotated_IR = FALSE) { 
    NxtIRF.CheckPackageInstalled("DESeq2", "1.30.0")
    .ASE_check_args(colData(se), 
        test_factor, test_nom, test_denom, batch1, batch2)
    BPPARAM_mod <- .validate_threads(n_threads)
    se_use <- .ASE_filter(
        se, filter_antiover, filter_antinear, filter_annotated_IR)
    # se_use <- realize_NxtSE(se_use)

    # Inc / Exc mode
    res.IncExc <- .ASE_DESeq2_contrast(se_use, 
        test_factor, test_nom, test_denom,
        batch1, batch2, BPPARAM_mod)
    res.inc = res.IncExc[grepl(".Included", get("EventName"))]
    res.inc[, c("EventName") := 
        sub(".Included","",get("EventName"), fixed=TRUE)]
    res.exc = res.IncExc[grepl(".Excluded", get("EventName"))]
    res.exc[, c("EventName") := 
        sub(".Excluded","",get("EventName"), fixed=TRUE)]

    # ASE mode
    rowData = as.data.frame(SummarizedExperiment::rowData(se_use))
    se_use = se_use[rowData$EventName %in% res.inc$EventName &
        rowData$EventName %in% res.exc$EventName,]

    res.ASE <- .ASE_DESeq2_contrast_ASE(se_use, 
        test_factor, test_nom, test_denom,
        batch1, batch2, BPPARAM_mod)
    res.ASE[res.inc, on = "EventName",
        paste("Inc", colnames(res.inc)[seq_len(6)], sep=".") := 
        list(get("i.baseMean"), get("i.log2FoldChange"), get("i.lfcSE"), 
            get("i.stat"), get("i.pvalue"), get("i.padj"))]
    res.ASE[res.exc, on = "EventName",
        paste("Exc", colnames(res.exc)[seq_len(6)], sep=".") := 
        list(get("i.baseMean"), get("i.log2FoldChange"), get("i.lfcSE"), 
            get("i.stat"), get("i.pvalue"), get("i.padj"))]
    res.ASE = res.ASE[!is.na(get("pvalue"))]
    setorder(res.ASE, "pvalue")
    res.ASE <- .ASE_add_diag(res.ASE, se_use, test_factor, test_nom, test_denom)
    return(res.ASE)
}

#' @describeIn ASE-methods Use DoubleExpSeq to perform differential ASE analysis 
#'   of a filtered NxtSE object (uses double exponential beta-binomial model)
#'   to estimate group dispersions, followed by LRT
#' @export
DoubleExpSeq_ASE <- function(se, test_factor, test_nom, test_denom, 
        # batch1 = "", batch2 = "",
        filter_antiover = TRUE, filter_antinear = FALSE, 
        filter_annotated_IR = FALSE) {
    
    NxtIRF.CheckPackageInstalled("DoubleExpSeq", "1.1")
    .ASE_check_args(colData(se), test_factor, 
        test_nom, test_denom, "", "")
    se_use <- .ASE_filter(
        se, filter_antiover, filter_antinear, filter_annotated_IR)

    res.ASE <- .ASE_DoubleExpSeq_contrast_ASE(se_use, 
        test_factor, test_nom, test_denom)
    
    res.cols = c(
        paste("MLE", test_nom, sep="_"), paste("MLE", test_denom, sep="_"),
        "P.Value", "adj.P.Val", "n_eff",
        paste("mDepth", test_nom, sep="_"), paste("mDepth", test_denom, sep="_"),
        "Dispersion_Reduced", "Dispersion_Full"
    )
    colnames(res.ASE)[-1] = res.cols
    
    res.ASE[, c("MLE_LFC") := (
        qlogis(res.ASE[,get(paste("MLE", test_nom, sep="_"))]) - 
        qlogis(res.ASE[,get(paste("MLE", test_denom, sep="_"))])
    ) / log(2)]
    
    res.ASE = res.ASE[, c("EventName", res.cols[1:2], "MLE_LFC",
        res.cols[3:9]), with = FALSE]

    res.ASE = res.ASE[!is.na(get("P.Value"))]
    setorderv(res.ASE, "pVal")
    res.ASE <- .ASE_add_diag(res.ASE, se_use, test_factor, test_nom, test_denom)
    return(res.ASE)
}



################################################################################
# helper functions:

.ASE_filter <- function(se, filter_antiover, filter_antinear,
        filter_annotated_IR) {
    se_use = se
    if(filter_antiover) {
        se_use = se_use[
            !grepl("anti-over", rowData(se_use)$EventName),]
    }
    if(filter_antinear) {
        se_use = se_use[
            !grepl("anti-near", rowData(se_use)$EventName),]
    }
    if(filter_annotated_IR) {
        se_use = se_use[
            !grepl("known-exon", rowData(se_use)$EventName),]
    }
    return(se_use)
}

.ASE_limma_contrast <- function(se, test_factor, test_nom, test_denom,
        batch1, batch2) {
    countData = rbind(assay(se, "Included"), 
        assay(se, "Excluded"))
    rowData = as.data.frame(rowData(se))
    colData = colData(se)
    rownames(colData) = colnames(se)
    colnames(countData) = rownames(colData)
    rownames(countData) = c(
        paste(rowData$EventName, "Included", sep="."),
        paste(rowData$EventName, "Excluded", sep=".")
    )
    
    condition_factor = factor(colData[, test_factor])
    if(batch2 != "") {    
        batch2_factor = colData[, batch2]
        batch1_factor = colData[, batch1]
        design1 = model.matrix(~0 + batch1_factor + batch2_factor + 
            condition_factor)
    } else if(batch1 != "") {
        batch1_factor = colData[, batch1]    
        design1 = model.matrix(~0 + batch1_factor + condition_factor)
    } else {
        design1 = model.matrix(~0 + condition_factor)    
    }
    contrast = rep(0, ncol(design1))
    contrast_a = paste0("condition_factor", test_nom)
    contrast_b = paste0("condition_factor", test_denom)
    contrast[which(colnames(design1) == contrast_b)] = -1
    contrast[which(colnames(design1) == contrast_a)] = 1

    countData_use = limma::voom(countData, design1)
    fit = limma::lmFit(countData_use$E, design = design1)
    fit = limma::contrasts.fit(fit, contrast)
    fit = limma::eBayes(fit)

    res = limma::topTable(fit, number = nrow(countData_use$E))
    res$EventName = rownames(res)

    res$AveExpr = res$AveExpr - min(res$AveExpr)
    res = as.data.table(res)
    return(res)
}

.ASE_limma_contrast_ASE <- function(se, test_factor, test_nom, test_denom,
        batch1, batch2) {
    countData = cbind(assay(se, "Included"), 
        assay(se, "Excluded"))

    rowData = as.data.frame(rowData(se))
    colData = as.data.frame(colData(se))
    colData = rbind(colData, colData)
    rownames(colData) = c(
        paste(colnames(se), "Included", sep="."),
        paste(colnames(se), "Excluded", sep=".")
    )
    colData$ASE = rep(c("Included", "Excluded"), each = ncol(se))
    colnames(countData) = rownames(colData)
    rownames(countData) = rowData$EventName
    
    condition_factor = factor(colData[, test_factor])
    ASE = colData[, "ASE"]
    if(batch2 != "") {    
        batch2_factor = colData[, batch2]
        batch1_factor = colData[, batch1]
        design1 = model.matrix(~0 + batch1_factor + batch2_factor + 
            condition_factor + condition_factor:ASE)
    } else if(batch1 != "") {
        batch1_factor = colData[, batch1]    
        design1 = model.matrix(~0 + batch1_factor + condition_factor + 
            condition_factor:ASE)
    } else {
        design1 = model.matrix(~0 + condition_factor + condition_factor:ASE)    
    }
    colnames(design1) = sub(":",".",colnames(design1))
    contrast = rep(0, ncol(design1))
    contrast_a = paste0("condition_factor", test_nom, ".ASEIncluded")
    contrast_b = paste0("condition_factor", test_denom, ".ASEIncluded")
    contrast[which(colnames(design1) == contrast_b)] = -1
    contrast[which(colnames(design1) == contrast_a)] = 1

    countData_use = limma::voom(countData, design1, lib.size = 1)

    fit = limma::lmFit(countData_use$E, design = design1)

    fit = limma::contrasts.fit(fit, contrast)
    fit = limma::eBayes(fit)

    res = limma::topTable(fit, number = nrow(countData_use$E))
    res$EventName = rownames(res)
    res$AveExpr = res$AveExpr - min(res$AveExpr)
    res = as.data.table(res)
    return(res)
}

.ASE_DESeq2_contrast <- function(se, test_factor, test_nom, test_denom,
        batch1, batch2, BPPARAM) {
    countData = rbind(assay(se, "Included"), 
        assay(se, "Excluded"))
    rowData = as.data.frame(rowData(se))
    colData = colData(se)
    rownames(colData) = colnames(se)
    colnames(countData) = rownames(colData)
    rownames(countData) = c(
        paste(rowData$EventName, "Included", sep="."),
        paste(rowData$EventName, "Excluded", sep=".")
    )
    if(batch2 != "") {
        dds_formula = paste0("~", paste(
            batch1, batch2, test_factor,
            sep="+"))

    } else if(batch1 != "") {
        dds_formula = paste0("~", paste(
            batch1, test_factor,
            sep="+"))
    } else {
        dds_formula = paste0("~", test_factor)
    }
    
    countData = as.matrix(countData)
    mode(countData) = "integer"
    dds = DESeq2::DESeqDataSetFromMatrix(
        countData = round(countData),
        colData = .DESeq_colData_sanitise(colData),
        design = as.formula(dds_formula)
    )
    message("DESeq_ASE: Profiling expression of Included and Excluded counts")
    
    dds = DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM)
    res = as.data.frame(DESeq2::results(dds,
        contrast = c(test_factor, test_nom, test_denom), 
        parallel = TRUE, BPPARAM = BPPARAM)
    )
    res$EventName = rownames(res)
    return(as.data.table(res))   
}

.ASE_DESeq2_contrast_ASE <- function(se, test_factor, test_nom, test_denom,
        batch1, batch2, BPPARAM) {
    countData = cbind(assay(se, "Included"), 
        assay(se, "Excluded"))
    rowData = as.data.frame(rowData(se))
    colData = as.data.frame(colData(se))
    colData = rbind(colData, colData)
    rownames(colData) = c(
        paste(colnames(se), "Included", sep="."),
        paste(colnames(se), "Excluded", sep=".")
    )
    colData$ASE = rep(c("Included", "Excluded"), each = ncol(se))
    colnames(countData) = rownames(colData)
    rownames(countData) = rowData$EventName
    
    if(batch2 != "") {
        dds_formula = paste0("~", paste(
            batch1, batch2, test_factor,
                paste0(test_factor, ":ASE"),
            sep="+"))
    } else if(batch1 != "") {
        dds_formula = paste0("~", paste(
            batch1, test_factor, 
            paste0(test_factor, ":ASE"),
            sep="+"))
    } else {
        dds_formula = paste0("~", paste(
            test_factor, 
            paste0(test_factor, ":ASE"),
            sep="+"))
    }
    countData = as.matrix(countData)
    mode(countData) = "integer"
    dds = DESeq2::DESeqDataSetFromMatrix(
        countData = countData,
        colData = .DESeq_colData_sanitise(colData),
        design = as.formula(dds_formula)
    )
    message("DESeq_ASE: Profiling differential ASE")
    
    dds = DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM)
    res = as.data.frame(DESeq2::results(dds,
        list(
            paste0(test_factor, test_nom, ".ASEIncluded"),
            paste0(test_factor, test_denom, ".ASEIncluded")
        ),
        parallel = TRUE, BPPARAM = BPPARAM)
    )
    res$EventName = rownames(res)
    return(as.data.table(res))    
}

.ASE_DoubleExpSeq_contrast_ASE <- function(se, test_factor, 
    test_nom, test_denom) {
    
    y = assay(se, "Included")
    m = assay(se, "Included") + assay(se, "Excluded")
    
    colData = as.data.frame(colData(se))
    groups = factor(colData[, test_factor])
    shrink.method="WEB"
    
    contrast.first = which(levels(groups) == test_nom)
    contrast.second = which(levels(groups) == test_denom)
    
    res = DoubleExpSeq::DBGLM1(
        as.maatrix(y), as.matrix(m), groups, shrink.method,
        contrast=c(contrast.first,contrast.second), 
        fdr.level=0.05, use.all.groups=TRUE)
    
    return(cbind(data.table(EventName = rownames(res$All)),
        as.data.table(res$All)))
}


.ASE_add_diag <- function(res, se, test_factor, test_nom, test_denom) {
    rowData = as.data.frame(rowData(se))
    rowData.DT = as.data.table(rowData[,
        c("EventName","EventType","EventRegion", "NMD_direction")])
    diag = make_diagonal(se, res$EventName, 
        test_factor, test_nom, test_denom)
    colnames(diag)[2:3] = c(paste0("AvgPSI_", test_nom), 
        paste0("AvgPSI_", test_denom))
    res = cbind(
        res[,c("EventName")], 
        as.data.table(diag[,2:3]), 
        res[,-c("EventName")]
    )
    res = rowData.DT[res, on = "EventName"]
    return(res)
}

.DESeq_colData_sanitise <- function(colData) {
    for(i in seq_len(ncol(colData))) {
        if(is(colData[,i], "character")) {
            colData[, i] <- factor(colData[, i])
        }
    }
    colData
}