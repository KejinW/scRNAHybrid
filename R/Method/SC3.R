library(knitr)
library(scater)


useSC3<-function(Data,ClusterNum){
  
  library("SC3")
  
  Data=input$data
  ClusterNum=4
  
  # cell expression
  tmp <- Data

  rownames(tmp) <- seq(1,nrow(tmp))
  colnames(tmp) <- seq(1,ncol(tmp))
  # SCESEt object
  
  source("./configure.R")
  
  sceset <- newSCESet(fpkmData = tmp, logExprsOffset = 1)
  sceset <- calculateQCMetrics(sceset,nmads = nmads, pct_feature_controls_threshold = pct_feature_controls_threshold)
  sceset <- sc3(sceset, ks = ClusterNum, biology = TRUE,
                pct_dropout_min = pct_dropout_min,pct_dropout_max = pct_dropout_max,
                d_region_min = d_region_min, d_region_max =d_region_max,svm_max = svm_max,
                kmeans_nstart=kmeans_nstart,kmeans_iter_max=kmeans_iter_max)
  
  
  detach("package:SC3", unload=TRUE)
  
return(sceset)

}