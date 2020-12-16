
limma_DT <- function(countData, colData, test_factor, test_nom, test_denom, batch1, batch2, useASE = FALSE) {
  assert_that(is_valid(test_factor) & is_valid(test_nom) & is_valid(test_denom),
    msg = "test_factor, test_nom, test_denom must be defined")
  assert_that(test_factor %in% colnames(colData),
    msg = "test_factor is not a condition in colData")
  assert_that(any(colData[, test_factor] == test_nom),
    msg = "test_nom is not found in any samples")
  assert_that(any(colData[, test_factor] == test_denom),
    msg = "test_denom is not found in any samples")
  if(!missing(batch1)) {
    assert_that(batch1 %in% colnames(colData),
      msg = "batch1 is not a condition in colData")
    assert_that(test_factor != batch1, msg = "batch1 and test_factor are the same")      
  }
  if(!missing(batch2)) {
    assert_that(batch2 %in% colnames(colData),
      msg = "batch1 is not a condition in colData")
    assert_that(test_factor != batch2, msg = "batch2 and test_factor are the same")      
  }
  if(!missing(batch1) & !missing(batch2)) {
    assert_that(batch2 != batch1, msg = "batch1 and batch2 are the same")  
  }
  
  # Contrast construction:
  if(!useASE) {
    condition_factor = colData[, test_factor]
    if(!missing(batch2)) {    
      batch2_factor = colData[, batch2]
      batch1_factor = colData[, batch1]
      design1 = model.matrix(~0 + batch1_factor + batch2_factor + condition_factor)
    } else if(!missing(batch1)) {
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
  } else {
    condition_factor = colData[, test_factor]
    ASE = colData[, "ASE"]
    if(!missing(batch2)) {    
      batch2_factor = colData[, batch2]
      batch1_factor = colData[, batch1]
      design1 = model.matrix(~0 + batch1_factor + batch2_factor + condition_factor + condition_factor:ASE)
    } else if(!missing(batch1)) {
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
  }

  countData_use = limma::voom(countData, design1, lib.size = 1)

  fit = limma::lmFit(countData_use$E, design = design1)

  fit = limma::contrasts.fit(fit, contrast)
  fit = limma::eBayes(fit)

  res.limma = limma::topTable(fit, number = nrow(countData_use$E))
  res.limma$EventName = rownames(res.limma)

  res.limma$AveExpr = res.limma$AveExpr - min(res.limma$AveExpr)
  
  return(as.data.table(res.limma))
}

#' @export
limma_optimised = function(se, DE_Var, nom_DE, denom_DE) {
    colData = as.data.frame(SummarizedExperiment::colData(se))
    assert_that(DE_Var %in% colnames(colData),
        msg = paste(DE_Var, "is not a valid parameter"))
    assert_that(nom_DE %in% unname(unlist(colData[,DE_Var])),
        msg = paste(nom_DE, "does not exist in column", DE_Var))
    assert_that(denom_DE %in% unname(unlist(colData[,DE_Var])),
        msg = paste(denom_DE, "does not exist in column", DE_Var))
    
    countData = cbind(SummarizedExperiment::assay(se, "Included"), 
        SummarizedExperiment::assay(se, "Excluded"))
    rowData = as.data.frame(SummarizedExperiment::rowData(se))
    
    colData_use = rbind(colData, colData)
    rownames(colData_use) = c(
        paste(rownames(colData), "Included", sep="."),
        paste(rownames(colData), "Excluded", sep=".")
    )
    colData_use$ASE = rep(c("Included", "Excluded"), each = nrow(colData))
    colnames(countData) = rownames(colData_use)
    rownames(countData) = rowData$EventName

    res.ASE = limma_DT(countData, colData_use, DE_Var, 
      nom_DE, denom_DE, useASE = TRUE)
    res.inc = limma_DT(countData[, seq_len(nrow(colData))], 
      colData_use[seq_len(nrow(colData)),], DE_Var, nom_DE, denom_DE)
    res.exc = limma_DT(countData[, seq(nrow(colData) + 1, nrow(colData) * 2)], 
      colData_use[seq(nrow(colData) + 1, nrow(colData) * 2),], 
        DE_Var, nom_DE, denom_DE)

    res.ASE[res.inc, on = "EventName",
      paste("Inc", colnames(res.inc)[1:6], sep=".") := 
        list(i.logFC, i.AveExpr, i.t, i.P.Value, i.adj.P.Val, i.B)]
      
    res.ASE[res.exc, on = "EventName",
      paste("Exc", colnames(res.exc)[1:6], sep=".") := 
        list(i.logFC, i.AveExpr, i.t, i.P.Value, i.adj.P.Val, i.B)]
      
    setorder(res.ASE, -B)
            
    rowData.DT = as.data.table(rowData[,c("EventName","EventType","EventRegion", "NMD_direction")])
 
    res.ASE = rowData.DT[res.ASE, on = "EventName"]
    
    # Apply custom filters
    res.ASE = res.ASE[Inc.AveExpr > 1 & Exc.AveExpr > 1]
    
    res.ASE
}