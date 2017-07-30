source("./R/function.R")
library(gplots)
library(lattice)
library(mclust)
library(readr)


Hybrid<-function(data){
  source("./configure.R")
  
  source("./R/Method/SIMLR.R")
  print("clustering with SIMLR")
  SIMLR = useSIMLR(Data=data, ClusterNum=c, large_scale = TRUE)
  
  source("./R/Method/CIDR.R")
  print("clustering with CIDR")
  CIDR<-useCIDR(Data=data)
  
  source("./R/Method/SC3.R")
  print("clustering with SC3")
  SC3 = useSC3(Data=data,  ClusterNum=c)
  
  source("./R/Method/Seurat.R")
  print("clustering with Seurat")
  Seurat <- inputSeurat(data)
  Seurat <- useSeurat(Seurat)
  
  
  hybrid<-list()
  
  hybrid$CIDR<-CIDR
  hybrid$SIMLR<-SIMLR
  hybrid$SC3<-SC3
  hybrid$Seurat<-Seurat
  
  
  #pdf("plots.pdf")
  #plotCIDR(hybrid$CIDR,Lab=lab)
  #plotSIMLR(hybrid$SIMLR,Lab=lab)
  #plotSC3(hybrid$SC3,Lab=lab)
  #plotSeurat(hybrid$Seurat,Lab=lab)
  #dev.off()
  
  #hybrid$CIDR$TFM<-getTrueFalseMatrix(CIDR@clusters,lab)
  #hybrid$SIMLR$TFM<-getTrueFalseMatrix(h$SIMLR$y$cluster,lab)
  #hybrid$SC3$TFM<-getTrueFalseMatrix(SC3[[21]],lab)
  #hybrid$Seurat$TFM<-getTrueFalseMatrix(Seurat@ident,lab)
  
  
  
  
  return(hybrid)
}



ConsistentCluster<-function(hybrid,plot=FALSE,PDF=NULL){
  
  cluster=list(hybrid$CIDR@clusters,hybrid$SIMLR$y$cluster,hybrid$SC3[[21]],hybrid$Seurat@ident)
  
  d <- matrix(unlist(cluster), ncol = length(cluster))
  
  #number of cluster vector
  n=ncol(d)
  
  #length of each cluster
  l=nrow(d)
  
  
  
  ccm=matrix(nrow=l, ncol=l)
  for (i in 1:dim(ccm)[1]){
    for (j in 1:dim(ccm)[2]){
      c1=d[i,]
      c2=d[j,] 
      
      #Similarity
      s=0
      
      for (k in 1:n){
        if (c1[k]==c2[k]){
          s=s+1
        }
      }
      
      s=s/n
      ccm[i,j]=s
      
    }
  }
  
  hybrid$ccm<-ccm
  
  
  if(plot){
    
    if (!is.null(PDF)){
      pdf(PDF)
    }
    
    print(levelplot(hybrid$ccm, scales=list(tck=0, x=list(rot=90)), col.regions=colorpanel(4,"white", "grey10"), 
            at=c(0,0.25,0.5,0.75,1), main="Consistent Cluster Matrix", xlab=NULL, ylab=NULL))
    
    
    if (!is.null(PDF)){
      dev.off()
    }
  
  }
  
  return(hybrid)
  
}


TPFNMatrix<-function(hybrid,Lab,plot=FALSE,PDF=NULL){
  #Lab=get(Lab)
  
  hybrid$TFM$CIDR<-getTrueFalseMatrix(hybrid$CIDR@clusters,Lab)
  hybrid$TFM$SIMLR<-getTrueFalseMatrix(hybrid$SIMLR$y$cluster,Lab)
  hybrid$TFM$SC3<-getTrueFalseMatrix(hybrid$SC3[[21]],Lab)
  hybrid$TFM$Seurat<-getTrueFalseMatrix(hybrid$Seurat@ident,Lab)
  
  if(plot){
    
    if (!is.null(PDF)){
      pdf(PDF)
    }
    
    plotTFM(hybrid$TFM$CIDR,main="CIDR")
    plotTFM(hybrid$TFM$SIMLR,main="SIMLR")
    plotTFM(hybrid$TFM$SC3,main="SC3")
    plotTFM(hybrid$TFM$Seurat,main="Seurat")
    
    
    if (!is.null(PDF)){
      dev.off()
    }
    
  }
  
  return(hybrid)
}



plotHybrid<-function(hybrid,Lab=NULL,PDF=NULL){
  if(!is.null(Lab)){
    #Lab=get(Lab)
  }
  
  if(!is.null(PDF)){
    
  pdf(PDF)
    
  }
  
    
  plotCIDR(hybrid$CIDR,Lab=Lab)
  plotSIMLR(hybrid$SIMLR,Lab=Lab)
  plotSC3(hybrid$SC3,Lab=Lab)
  plotSeurat(hybrid$Seurat,Lab=Lab)
  
  
  if(!is.null(PDF)){
  dev.off()
  
  }
}


#' Adjusted Rand Index
#'

#' @param cluVec predict clustering vector
#' @param Truelab true clustering vector
#' @importFrom both cluVec and Truelab should be a vector
#' @return numeric value -1:1, ARI of a random cluster will near to 0
#' @export

ARI<-function(hybrid,Lab,plot=FALSE,PDF=NULL){
  
  #Lab=get(Lab)
  
  Lab=unlist(Lab)
  
  hybrid$ARI$CIDR=adjustedRandIndex(hybrid$CIDR@clusters,Lab)
  
  hybrid$ARI$SIMLR=adjustedRandIndex(hybrid$SIMLR$y$cluster,Lab)
  
  hybrid$ARI$SC3=adjustedRandIndex(hybrid$SC3[[21]],Lab)
  
  hybrid$ARI$Seurat=adjustedRandIndex(hybrid$Seurat@ident,Lab)
  
  if(plot){
    
    if (!is.null(PDF)){
      pdf(PDF)
    }
    
    ari=c(hybrid$ARI$CIDR,hybrid$ARI$SIMLR, hybrid$ARI$SC3,hybrid$ARI$Seurat)
    barplot(ari,horiz=TRUE,main="Adjusted Rand Index",names.arg=c("CIDR", "SIMLR", "SC3","Seurat"))
    
    
    
    if (!is.null(PDF)){
      dev.off()
    }
    
  }
  
  
  
  return(hybrid)

}


dataTrans<-function(csv,Lab=TRUE){
  data=read_csv(csv)
  
  if(Lab){
    d<- as.matrix(data[-1,-1])
    gene <- as.vector(t(data[-1,1]))
    rownames(d) <-  gene
  
    label <-as.data.frame(t(data[1,-1]))
    
    input=list()
    input$data=d
    input$label=label
    
    return(input)
  
  }
  
  else{
    d<- as.matrix(data[,-1])
    gene <- as.vector(t(data[,1]))
    rownames(d) <-  gene
    
    input=list()
    input$data=d
    input$label=NULL
    
    return(input)
  }

  
  
}




