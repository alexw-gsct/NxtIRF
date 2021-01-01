limma_assert <- function(colData, test_factor, test_nom, test_denom, batch1, batch2) {
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
}

#' @export
limma_ASE <- function(se, test_factor, test_nom, test_denom, batch1 = "", batch2 = "",
    filter_antiover = TRUE, filter_antinear = FALSE, filter_anootated_IR = FALSE) {
    limma_assert(SummarizedExperiment::colData(se), test_factor, test_nom, test_denom, batch1, batch2)

    se_use = se
    if(filter_antiover) {
        se_use = se_use[grepl("anti-over", SummarizedExperiment::rowData(se_use)$EventName),]
    }
    if(filter_antinear) {
        se_use = se_use[grepl("anti-near", SummarizedExperiment::rowData(se_use)$EventName),]
    }
    if(filter_anootated_IR) {
        se_use = se_use[grepl("known-exon", SummarizedExperiment::rowData(se_use)$EventName),]
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
      paste("Inc", colnames(res.inc)[1:6], sep=".") := 
        list(i.logFC, i.AveExpr, i.t, i.P.Value, i.adj.P.Val, i.B)]
      
    res.ASE[res.exc, on = "EventName",
      paste("Exc", colnames(res.exc)[1:6], sep=".") := 
        list(i.logFC, i.AveExpr, i.t, i.P.Value, i.adj.P.Val, i.B)]
      
    setorder(res.ASE, -B)
            
    rowData.DT = as.data.table(rowData[,c("EventName","EventType","EventRegion", "NMD_direction")])
 
    res.ASE = rowData.DT[res.ASE, on = "EventName"]

    
    res.ASE
}

