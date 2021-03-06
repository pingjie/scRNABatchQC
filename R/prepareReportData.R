prepareReportData <- function(sampleTable, organism, cachePrefix) {
  cat("Preparing sample statistics data ...\n")
  if (!missing(cachePrefix)) {
    scesFile <- paste0(cachePrefix, "_sces.rdata")
    if (file.exists(scesFile)) {
      cat("Loading from cache file", scesFile, " ...\n")
      load(scesFile)
    } else {
      sces <- prepareSCRNADataSet(sampleTable, organism)
      save(sces, file=scesFile)
    }
  } else {
    sces <- prepareSCRNADataSet(sampleTable, organism)
  }
  
  cat("Preparing PCA and tSNE data ...\n")
  if (!missing(cachePrefix)) {
    scesAllFile <- paste0(cachePrefix, "_scesall.rdata")
    if (file.exists(scesAllFile)) {
      cat("Loading from cache file", scesAllFile, " ...\n")
      load(scesAllFile)
    } else {
      scesall <- preparePCATSNEData(sces)
      save(scesall, file=scesAllFile)
    }
  } else {
    scesall <- preparePCATSNEData(sces)
  }

  cat("Preparing differential expression analysis data ...\n")
  if (!missing(cachePrefix)) {
    diffFCFile <- paste0(cachePrefix, "_diffFC.rdata")
    if (file.exists(diffFCFile)) {
      cat("Loading from cache file", diffFCFile, " ...\n")
      load(diffFCFile)
    } else {
      diffFC <- .getDiffGenes(scesall, organism = organism,  FDR = 0.01, geneNo = 50)
      save(diffFC, file = diffFCFile)
    }
  }else{
    diffFC <- .getDiffGenes(scesall, organism = organism,  FDR = 0.01, geneNo = 50)
  }
  
  hvgPathways <- .getMultiplePathway(sces, metaObjectName = "hvgPathway")
  pc1Pathways <- .getMultiplePathway(sces, metaObjectName = "pc1Pathway")
  hvgBiologicalSimilarity <- .getBiologicalSimilarity(sces, objectName = "hvg", 
                                                      filterName = "FDR", valueName = "bio")
  pc1geneBiologicalSimilarity <- .getBiologicalSimilarity(sces, objectName = "pc1genes",
                                                          filterName = "adj.P.Val", valueName = "logFC")

  pw <- .prepareTableSummary(sces)
  
  plotData <- list(sces = sces, 
                   scesall = scesall, 
                   diffFC = diffFC, 
                   hvgPathways = hvgPathways, 
                   pc1Pathways = pc1Pathways,
                   hvgBiologicalSimilarity = hvgBiologicalSimilarity,
                   pc1geneBiologicalSimilarity = pc1geneBiologicalSimilarity,
                   tableSummary = pw)
  
  cat("Report data prepared.\n")
  return(plotData)
}
