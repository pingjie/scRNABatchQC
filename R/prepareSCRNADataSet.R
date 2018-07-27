##' prepareSCRNADataSet
##'
##' The function read multiple count table and prepare statistics information for each table
##'
##' @param sampleTable the sample table with first column as sample name, second column as file location and an optional third column as transform
##' @return a named list of makeSCRNAdata result
##' @export prepareSCRNADataSet
##' @examples 
##' #sampleTable<-data.frame(Sample=c("S1", "S2", "S3"),
##' #                        File=file.path(("Z:/shengq1/20180214_scRNABatchQC/", c("S1.csv", "S2.csv", "S3.csv")),
##' #                        Transform=c(1,1,1))
##' #sces<-prepareSCRNADataSet(sampleTable)
prepareSCRNADataSet <- function(sampleTable, organism){
  result <- list()
  
  n <- 1
  for (n in 1:nrow(sampleTable)) {
    sampleName<-as.character(sampleTable[n,1])
    countFile<-as.character(sampleTable[n,2])
    cat("Preparing ", sampleName, "\n")
    counts<-as.matrix(read.csv(countFile, row.names=1, header=F))
    result[[n]] <- prepareSCRNAData(counts, organism)
  }
  names(result) <- sampleTable[, 1]
  return(result)
}

##' preparePCATSNEData
##'
##' The function prepare statistics information from multiple scRNA dataset.
##'
##' @param sces a named list of makeSCRNAdata result
##' @return a sce:SingleCellExperiment data with PCA and TSNE
##' @importFrom SingleCellExperiment SingleCellExperiment reducedDim
##' @importFrom Scater calculateQCMetrics isOutlier calcAverage nexprs normalize runPCA .get_palette
##' @importFrom Scran quickCluster computeSumFactors trendVar decomposeVar
##' @importFrom Rtsne Rtsne
##' @importFrom data.table data.table
##' @export preparePCATSNEData
##' @examples 
##' #sces <- prepareSCRNADataSet(sampleTable)
##' #sceall <- preparePCATSNEData(sces)
preparePCATSNEData <- function(sces, ncomponents = 10, perplexity = 20) {
	dt <- data.table(cbind(rownames(counts(sces[[1]]$sce)), counts(sces[[1]]$sce)), key = "V1")

	conditions <- rep(names(sces)[1], dim(sces[[1]]$sce)[2])

	if(length(sces) > 1){
		for (i in 2:length(sces)) {
			dt2 <- data.table(cbind(rownames(counts(sces[[i]]$sce)), counts(sces[[i]]$sce)), key = "V1")
			dt <- merge(dt, dt2, by = "V1", all = T)
			conditions <- c(conditions, rep(names(sces)[i], dim(sces[[i]]$sce)[2]))
		}
	}
	dt[is.na(dt)] <- 0
	ac <- apply(dt[, -1], 2, as.numeric)
	rownames(ac) <- dt$V1

	sceall <- SingleCellExperiment(list(counts = (as.matrix(ac))))
	colData(sceall)$condition <- conditions

	ave.counts <- calcAverage(sceall)
	high.ave <- ave.counts >= 0.1
	clusters <- quickCluster(sceall, subset.row = high.ave, method = "igraph")
	sceall <- computeSumFactors(sceall, cluster = clusters, subset.row = high.ave, min.mean = 0, positive=TRUE)
	sceall <- normalize(sceall)

	sceall <- runPCA(sceall, ncomponents = ncomponents)

	set.seed(100)

	tsne_out <- Rtsne(reducedDim(sceall, "PCA"), initial_dims = ncol(reducedDim(sceall, "PCA")), pca = FALSE, perplexity = perplexity)
	reducedDim(sceall, "TSNE") <- tsne_out$Y

	return(sceall)
}
