library(Matrix)
library(parallel)
library(igraph)
library(grDevices)


useSIMLR<-function(Data,ClusterNum,large_scale){
  
  library(SIMLR)
  
  #Data=input$data
  #ClusterNum=4
  
  source("./configure.R")
  set.seed(11111)
  #ClusterNum should be change
  if (large_scale == TRUE) {
    res_example1 = SIMLR_Large_Scale(X=Data,c=ClusterNum,k = k)
  }
  else {
    res_example1 = SIMLR(X=Data,c=ClusterNum,k = k, cores.ratio = cores.ratio)
  }
  
  detach("package:SIMLR", unload=TRUE)
  
  return(res_example1)
}