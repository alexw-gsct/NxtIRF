DE_assert <- function(colData, test_factor, test_nom, test_denom, batch1, batch2) {
  assert_that(is_valid(test_factor) & is_valid(test_nom) & is_valid(test_denom),
    msg = "test_factor, test_nom, test_denom must be defined")
  assert_that(test_factor %in% colnames(colData),
    msg = "test_factor is not a condition in colData")
  assert_that(any(colData[, test_factor] == test_nom),
    msg = "test_nom is not found in any samples")
  assert_that(any(colData[, test_factor] == test_denom),
    msg = "test_denom is not found in any samples")
  if(batch1 != "") {
    assert_that(batch1 %in% colnames(colData),
      msg = "batch1 is not a condition in colData")
    assert_that(test_factor != batch1, msg = "batch1 and test_factor are the same")      
  }
  if(batch2 != "") {
    assert_that(batch2 %in% colnames(colData),
      msg = "batch1 is not a condition in colData")
    assert_that(test_factor != batch2, msg = "batch2 and test_factor are the same")      
  }
  if(batch1 != "" & batch2 != "") {
    assert_that(batch2 != batch1, msg = "batch1 and batch2 are the same")  
  }
  return(TRUE)
}

#' Use Limma to test for differential ASE (Alternative Splice Event)
#'
#' @param se The SummarizedExperiment object created by `MakeSE()`. To reduce runtime and false negatives
#'   due to multiple testing issues, please filter the object using `apply_filter()`
#' @param test_factor A string for the column name which contains the contrasting variable
#' @param test_nom The condition in which to test for differential ASE. Usually the "treatment" condition
#' @param test_denom The condition in which to test against for differential ASE. Usually the "control" condition
#' @param batch1,batch2 One or two columns containing batch information to normalise against (can be omitted)
#' @param filter_antiover Whether to filter out IR events that overlap antisense genes (for unstranded RNAseq protocols)
#' @param filter_antinear Whether to filter out IR events near but not overlapping antisense genes 
#'   (for unstranded RNAseq protocols)
#' @param filter_annotated_IR Whether to filter out IR events that are already annotated exons
#'   (after doing so, all IR events will be unannotated - i.e. constitutionally spliced introns))
#' @return A data table containing the following:
#'   EventName: The name of the ASE event\cr\cr
#'   EventType: The type of event. IR = intron retention, MXE = mutually exclusive event, SE = skipped exon,
#'     AFE = alternate first exon, ALE = alternate last exon, A5SS / A3SS = alternate 5' / 3' splice site\cr\cr
#'   EventRegion: The genomic coordinates the event occupies.\cr\cr  
#'   NMD_direction: Indicates whether one isoform is a NMD substrate. +1 means included isoform is NMD, 
#'     -1 means the excluded isoform is NMD, and 0 means neither (or both) are NMD substrates\cr\cr
#'   AvgPSI_nom, Avg_PSI_denom: the average percent spliced in / percent intron retention levels for the
#'   two conditions being contrasted\cr\cr
#'   logFC, AveExpr, t, P.Value, adj.P.Val, B: limma topTable columns of limma results. See `?limma::topTable`\cr\cr
#'   inc/exc_(logFC, AveExpr, t, P.Value, adj.P.Val, B): limma results for differential testing for
#'     raw included / excluded counts only\cr\cr
#' @export
limma_ASE <- function(se, test_factor, test_nom, test_denom, batch1 = "", batch2 = "",
    filter_antiover = TRUE, filter_antinear = FALSE, filter_annotated_IR = FALSE) {
    
    NxtIRF.CheckPackageInstalled("limma", "3.44.0")
    test_assert = FALSE
    test_assert = DE_assert(SummarizedExperiment::colData(se), test_factor, test_nom, test_denom, batch1, batch2)

    if(!test_assert) {
        return(NULL)
    }

    se_use = se
    if(filter_antiover) {
        se_use = se_use[!grepl("anti-over", SummarizedExperiment::rowData(se_use)$EventName),]
    }
    if(filter_antinear) {
        se_use = se_use[!grepl("anti-near", SummarizedExperiment::rowData(se_use)$EventName),]
    }
    if(filter_annotated_IR) {
        se_use = se_use[!grepl("known-exon", SummarizedExperiment::rowData(se_use)$EventName),]
    }
    
    # Inc / Exc mode
    countData = rbind(SummarizedExperiment::assay(se_use, "Included"), 
        SummarizedExperiment::assay(se_use, "Excluded"))
    rowData = as.data.frame(SummarizedExperiment::rowData(se_use))
    
    colData = SummarizedExperiment::colData(se_use)
    rownames(colData) = colnames(se_use)
    colnames(countData) = rownames(colData)
    rownames(countData) = c(
        paste(rowData$EventName, "Included", sep="."),
        paste(rowData$EventName, "Excluded", sep=".")
    )
    
    condition_factor = factor(colData[, test_factor])
    if(batch2 != "") {    
      batch2_factor = colData[, batch2]
      batch1_factor = colData[, batch1]
      design1 = model.matrix(~0 + batch1_factor + batch2_factor + condition_factor)
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

    res.limma2 = limma::topTable(fit, number = nrow(countData_use$E))
    res.limma2$EventName = rownames(res.limma2)

    res.limma2$AveExpr = res.limma2$AveExpr - min(res.limma2$AveExpr)
    res.limma2 = as.data.table(res.limma2)

    res.inc = res.limma2[grepl(".Included", EventName)]
    res.inc[, EventName := sub(".Included","",EventName, fixed=TRUE)]
    res.inc = res.inc[AveExpr > 1]
    res.exc = res.limma2[grepl(".Excluded", EventName)]
    res.exc[, EventName := sub(".Excluded","",EventName, fixed=TRUE)]
    res.exc = res.exc[AveExpr > 1]

    # ASE mode
    rowData = as.data.frame(SummarizedExperiment::rowData(se_use))
    se_use = se_use[rowData$EventName %in% res.inc$EventName &
        rowData$EventName %in% res.exc$EventName,]
    rowData = as.data.frame(SummarizedExperiment::rowData(se_use))
    countData = cbind(SummarizedExperiment::assay(se_use, "Included"), 
        SummarizedExperiment::assay(se_use, "Excluded"))

    colData = as.data.frame(SummarizedExperiment::colData(se_use))
    colData = rbind(colData, colData)
    rownames(colData) = c(
        paste(colnames(se_use), "Included", sep="."),
        paste(colnames(se_use), "Excluded", sep=".")
    )
    colData$ASE = rep(c("Included", "Excluded"), each = ncol(se_use))
    colnames(countData) = rownames(colData)
    rownames(countData) = rowData$EventName
    
    condition_factor = factor(colData[, test_factor])
    ASE = colData[, "ASE"]
    if(batch2 != "") {    
      batch2_factor = colData[, batch2]
      batch1_factor = colData[, batch1]
      design1 = model.matrix(~0 + batch1_factor + batch2_factor + condition_factor + condition_factor:ASE)
    } else if(batch1 != "") {
      batch1_factor = colData[, batch1]    
      design1 = model.matrix(~0 + batch1_factor + condition_factor + condition_factor:ASE)
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

    res.limma = limma::topTable(fit, number = nrow(countData_use$E))
    res.limma$EventName = rownames(res.limma)

    res.limma$AveExpr = res.limma$AveExpr - min(res.limma$AveExpr)
    res.limma = as.data.table(res.limma)

    # Merge tables together:
    res.ASE = res.limma
    
    res.ASE[res.inc, on = "EventName",
      paste("Inc", colnames(res.inc)[seq_len(6)], sep=".") := 
        list(i.logFC, i.AveExpr, i.t, i.P.Value, i.adj.P.Val, i.B)]
      
    res.ASE[res.exc, on = "EventName",
      paste("Exc", colnames(res.exc)[seq_len(6)], sep=".") := 
        list(i.logFC, i.AveExpr, i.t, i.P.Value, i.adj.P.Val, i.B)]
      
    setorder(res.ASE, -B)
            
    rowData.DT = as.data.table(rowData[,c("EventName","EventType","EventRegion", "NMD_direction")])

    diag = make_diagonal(se, res.ASE$EventName, test_factor, test_nom, test_denom)
    colnames(diag)[2:3] = c(paste0("AvgPSI_", test_nom), paste0("AvgPSI_", test_denom))
    
    res.ASE = cbind(res.ASE[,c("EventName")], as.data.table(diag[,2:3]), res.ASE[,-c("EventName")])
 
    res.ASE = rowData.DT[res.ASE, on = "EventName"]

    res.ASE
}

#' Use DESeq2 to test for differential ASE (Alternative Splice Event)
#'
#' @param se The SummarizedExperiment object created by `MakeSE()`. To reduce runtime and false negatives
#'   due to multiple testing issues, please filter the object using `apply_filter()`
#' @param test_factor A string for the column name which contains the contrasting variable
#' @param test_nom The condition in which to test for differential ASE. Usually the "treatment" condition
#' @param test_denom The condition in which to test against for differential ASE. Usually the "control" condition
#' @param batch1,batch2 One or two columns containing batch information to normalise against (can be omitted)
#' @param n_threads The number of threads to use for DESeq2
#' @param filter_antiover Whether to filter out IR events that overlap antisense genes (for unstranded RNAseq protocols)
#' @param filter_antinear Whether to filter out IR events near but not overlapping antisense genes 
#'   (for unstranded RNAseq protocols)
#' @param filter_annotated_IR Whether to filter out IR events that are already annotated exons
#'   (after doing so, all IR events will be unannotated - i.e. constitutionally spliced introns))
#' @return A data table containing the following:
#'   EventName: The name of the ASE event\cr\cr
#'   EventType: The type of event. IR = intron retention, MXE = mutually exclusive event, SE = skipped exon,
#'     AFE = alternate first exon, ALE = alternate last exon, A5SS / A3SS = alternate 5' / 3' splice site\cr\cr
#'   EventRegion: The genomic coordinates the event occupies.\cr\cr  
#'   NMD_direction: Indicates whether one isoform is a NMD substrate. +1 means included isoform is NMD, 
#'     -1 means the excluded isoform is NMD, and 0 means neither (or both) are NMD substrates\cr\cr
#'   AvgPSI_nom, Avg_PSI_denom: the average percent spliced in / percent intron retention levels for the
#'   two conditions being contrasted\cr\cr
#'   baseMean, log2FoldChange, lfcSE, stat, pvalue, padj: 
#'   DESeq2 results columns See `?DESeq2::results`\cr\cr
#'   inc/exc_(baseMean, log2FoldChange, lfcSE, stat, pvalue, padj): 
#'   DESeq2 results for differential testing for
#'   raw included / excluded counts only\cr\cr
#' @export
DESeq_ASE <- function(se, test_factor, test_nom, test_denom, batch1 = "", batch2 = "",
    n_threads = 1,
    filter_antiover = TRUE, filter_antinear = FALSE, filter_annotated_IR = FALSE) {
    
    NxtIRF.CheckPackageInstalled("DESeq2", "1.30.0")
    test_assert = FALSE
    test_assert = DE_assert(colData(se), 
        test_factor, test_nom, test_denom, batch1, batch2)

    if(!test_assert) {
        return(NULL)
    }
    
    BPPARAM = BiocParallel::bpparam()

    if(n_threads == 1) {
      BPPARAM_mod = BiocParallel::SerialParam()
      message(paste("DESeq_ASE:", "Using SerialParam - 1 thread"))
    } else if(Sys.info()["sysname"] == "Windows") {
      BPPARAM_mod = BiocParallel::SnowParam(n_threads)
      message(paste("DESeq_ASE:", "Using SnowParam", BPPARAM_mod$workers, "threads"))
    } else {
      BPPARAM_mod = BiocParallel::MulticoreParam(n_threads)
      message(paste("DESeq_ASE:", "Using MulticoreParam", 
        BPPARAM_mod$workers, "threads"))
    }

    se_use = se
    if(filter_antiover) {
        se_use = se_use[!grepl("anti-over", rowData(se_use)$EventName),]
    }
    if(filter_antinear) {
        se_use = se_use[!grepl("anti-near", rowData(se_use)$EventName),]
    }
    if(filter_annotated_IR) {
        se_use = se_use[!grepl("known-exon", rowData(se_use)$EventName),]
    }
    
    # Inc / Exc mode
    countData = rbind(assay(se_use, "Included"), 
        assay(se_use, "Excluded"))
    rowData = as.data.frame(rowData(se_use))
    
    colData = colData(se_use)
    rownames(colData) = colnames(se_use)
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
        colData = colData,
        design = as.formula(dds_formula)
    )
    message("DESeq_ASE: Profiling expression of Included and Excluded counts")
    
    dds = DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM_mod)

    res.IncExc = as.data.frame(DESeq2::results(dds,
        contrast = c(test_factor, test_nom, test_denom), 
        parallel = TRUE, BPPARAM = BPPARAM_mod)
    )
    
    res.IncExc$EventName = rownames(res.IncExc)

    res.IncExc = as.data.table(res.IncExc)

    res.inc = res.IncExc[grepl(".Included", EventName)]
    res.inc[, EventName := sub(".Included","",EventName, fixed=TRUE)]
    # res.inc = res.inc[baseMean > 1]
    res.exc = res.IncExc[grepl(".Excluded", EventName)]
    res.exc[, EventName := sub(".Excluded","",EventName, fixed=TRUE)]
    # res.exc = res.exc[baseMean > 1]

    # ASE mode
    rowData = as.data.frame(SummarizedExperiment::rowData(se_use))
    se_use = se_use[rowData$EventName %in% res.inc$EventName &
        rowData$EventName %in% res.exc$EventName,]
    rowData = as.data.frame(SummarizedExperiment::rowData(se_use))
    countData = cbind(SummarizedExperiment::assay(se_use, "Included"), 
        SummarizedExperiment::assay(se_use, "Excluded"))

    colData = as.data.frame(SummarizedExperiment::colData(se_use))
    colData = rbind(colData, colData)
    rownames(colData) = c(
        paste(colnames(se_use), "Included", sep="."),
        paste(colnames(se_use), "Excluded", sep=".")
    )
    colData$ASE = rep(c("Included", "Excluded"), each = ncol(se_use))
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
        colData = colData,
        design = as.formula(dds_formula)
    )
    message("DESeq_ASE: Profiling differential ASE")
    
    dds = DESeq2::DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM_mod)

    res.ASE = as.data.frame(DESeq2::results(dds,
        list(
            paste0(test_factor, test_nom, ".ASEIncluded"),
            paste0(test_factor, test_denom, ".ASEIncluded")
        ),
        parallel = TRUE, BPPARAM = BPPARAM_mod)
    )

    res.ASE$EventName = rownames(res.ASE)

    res.ASE = as.data.table(res.ASE)
    
    res.ASE[res.inc, on = "EventName",
      paste("Inc", colnames(res.inc)[seq_len(6)], sep=".") := 
        list(i.baseMean, i.log2FoldChange, i.lfcSE, i.stat, i.pvalue, i.padj)]
      
    res.ASE[res.exc, on = "EventName",
      paste("Exc", colnames(res.exc)[seq_len(6)], sep=".") := 
        list(i.baseMean, i.log2FoldChange, i.lfcSE, i.stat, i.pvalue, i.padj)]
    
    res.ASE = res.ASE[!is.na(pvalue)]
    
    setorder(res.ASE, pvalue)
            
    rowData.DT = as.data.table(rowData[,c("EventName","EventType","EventRegion", "NMD_direction")])
 
    # Add average PIR / PSI values
    
    diag = make_diagonal(se, res.ASE$EventName, test_factor, test_nom, test_denom)
    colnames(diag)[2:3] = c(paste0("AvgPSI_", test_nom), paste0("AvgPSI_", test_denom))
    
    res.ASE = cbind(res.ASE[,c("EventName")], as.data.table(diag[,2:3]), res.ASE[,-c("EventName")])

    res.ASE = rowData.DT[res.ASE, on = "EventName"]

    res.ASE
}