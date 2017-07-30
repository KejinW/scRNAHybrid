library(Matrix)
library(dplyr)

useSeurat<-function(data){
  
  library(Seurat)
  
  ?#pbmc.data <- Read10X("~/Downloads/outs/filtered_gene_bc_matrices/hg19/")
  source("./configure.R")
  
  pbmc.data<-data
  pbmc <- new("seurat", raw.data = pbmc.data)
  pbmc <- Setup(pbmc, min.cells = min.cells, min.genes = min.genes, do.logNormalize = F, total.expr = 1e4, project = "10X_PBMC")
  
  #pbmc <- SubsetData(pbmc, subset.name = "nGene", accept.low = 1700)
  # This removes unwanted sources of variation
  pbmc <- RegressOut(pbmc, latent.vars = c("nUMI", "nGene"), do.scale = do.scale, do.center = do.center)
  
  pbmc <- MeanVarPlot(pbmc ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = x.low.cutoff, x.high.cutoff = x.high.cutoff, y.cutoff = y.cutoff, do.contour = F)
  
  pbmc <- PCA(pbmc, pc.genes = pbmc@var.genes, do.print = FALSE, pcs.print = 5, genes.print = 5)
  
  pbmc <- FindClusters(pbmc, pc.use = pc.use, resolution = resolution, print.output = 0, save.SNN = T)
  
  pbmc <- RunTSNE(pbmc, dims.use = dims.use, do.fast = T)
  
  
  #pbmc.data<-data
  #pbmc <- new("seurat", raw.data = pbmc.data)
  #pbmc <- Setup(pbmc, min.cells = 3, min.genes = 200, do.logNormalize = F, total.expr = 1e4, project = "10X_PBMC")
  #pbmc <- SubsetData(pbmc, subset.name = "nGene", accept.low = 1700)
  # This removes unwanted sources of variation
  #pbmc <- RegressOut(pbmc, latent.vars = c("nUMI", "nGene"), do.scale = T, do.center = T)
  #pbmc <- MeanVarPlot(pbmc ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)
  #pbmc <- PCA(pbmc, pc.genes = pbmc@var.genes, do.print = FALSE, pcs.print = 5, genes.print = 5
  #pbmc <- JackStraw(pbmc, num.replicate = 100, do.print = FALSE)
  #pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = 1, print.output = 0, save.SNN = T)
  #pbmc <- RunTSNE(pbmc, dims.use = 1:10, do.fast = T)
  
  detach("package:Seurat", unload=TRUE)
  
  return(pbmc)
}

inputSeurat<-function(data){
  
  d<-as(data,"dgTMatrix")
  
  row=as.character(seq(1,nrow(data)))
  col=as.character(seq(1,ncol(data)))

  d@Dimnames<-list(row,col)
  
  return(d)
}
