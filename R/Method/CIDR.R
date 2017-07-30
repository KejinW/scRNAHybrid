useCIDR<-function(Data){
  
library(cidr)
  
source("./configure.R")
  
tags <- as.matrix(Data)

sData <- scDataConstructor(tags)

sData <- determineDropoutCandidates(sData,min1 =min1, min2 = min2,
                                    N = N, alpha = alpha)

sData <- wThreshold(sData,cutoff = cutoff)

sData <- scDissim(sData)

sData <- scPCA(sData)

sData <- nPC(sData)

nCluster(sData)

sData <- scCluster(sData)

detach("package:cidr", unload=TRUE)

return(sData)

}

