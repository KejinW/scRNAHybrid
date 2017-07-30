#' get consistent cluster
#'
#'This function aim to find the clustering relation 
#'between two cell in different method
#'if all of the methods assign two cell into same clustering
#'the consistent of those two cell equal to 1
#' @param cluster predict clustering vector
#' @importFrom cluster should be a list include all the clustering vector
#' @return consistent matrix
#' @export


  



#' get True False Positive Negative Matrix
#' two cells were assigned to same of different cluster
#' TP predict cluster:same, true cluster:same
#' TN predict cluster:different, true cluster:different
#' FP predict cluster:different, true cluster:same
#' FN predict cluster:same, true cluster:different
#'
#' @param cluVec predict clustering vector
#' @param Truelab true clustering vector
#' @importFrom both cluVec and Truelab should be a vector
#' @return true false matrix
#' @export


getTrueFalseMatrix<-function(cluVec,Lab){
  cluVec=as.matrix(cluVec)
  Lab=as.matrix(Lab)
  
  l=nrow(cluVec)
  
  tfm=matrix(nrow=l, ncol=l)
  
  for (i in 1:dim(tfm)[1]){
    for (j in 1:dim(tfm)[2]){
      c1=cluVec[i,]
      c2=cluVec[j,]
      
      l1=Lab[i,]
      l2=Lab[j,]
      
      if(c1==c2&l1==l2){
        tfm[i,j]="TP"
      }
      
      else if(c1==c2&l1!=l2){
        tfm[i,j]="FP"
      }
      
      else if(c1!=c2&l1!=l2){
        tfm[i,j]="TN"
      }
      
      else if(c1!=c2&l1==l2){
        tfm[i,j]="FN"
      }
      
    }
  }
  
  return(tfm)
  
}



plotSeurat<-function(Seurat,Lab=NULL){
  
  if(!is.null(Lab)){
    plot(Seurat@tsne.rot,col=Lab[,1]*3,pch=as.numeric(Seurat@ident),xlab="TSNE1", ylab="TSNE2",main="Seurat")
  }
  
  else{
    plot(Seurat@tsne.rot,pch=as.numeric(Seurat@ident),xlab="TSNE1", ylab="TSNE2",main="Seurat")
  }
  
}


plotSC3<-function(SC3,Lab=NULL){
  
  p=plotPCA(SC3)
  
  pc1=p$data$PC1
  pc2=p$data$PC2
  
  cluNum=SC3@sc3$ks
  clu=paste("sc3_",toString(cluNum),"_clusters",sep="")
  
  
  
  if(!is.null(Lab)){
    plot(x=pc1,y=pc2,xlab="PC1", ylab="PC2",col=Lab[,1]*3,pch=as.numeric(SC3@phenoData@data[[21]]),main="SC3")
  }
  
  
  else{
    plot(x=pc1,y=pc2,xlab="PC1", ylab="PC2",pch=SC3@phenoData@data$sc3_4_clusters,main="SC3")
  }
  
}


plotCIDR<-function(CIDR,Lab=NULL){
  
  if(!is.null(Lab)){
    plot(CIDR@PC[,c(1,2)],col=Lab[,1]*3,pch=CIDR@clusters, main="CIDR", xlab="PC1", ylab="PC2")
  }
  
  else{
    plot(CIDR@PC[,c(1,2)],pch=CIDR@clusters, main="CIDR", xlab="PC1", ylab="PC2")
  }
  
}


plotSIMLR<-function(SIMLR,Lab=NULL){
  
  if(!is.null(Lab)){
    plot(SIMLR$ydata,col=Lab[,1]*3,xlab="SIMLR component 1",pch=SIMLR$y$cluster, ylab="SIMLR component 2",main="SIMILR")
  }
  
  else{
    plot(SIMLR$ydata,xlab="SIMLR component 1",pch=SIMLR$y$cluster, ylab="SIMLR component 2",main="SIMILR")
    
  }
}

plotTFM<-function(mat,main=NA,xlab=NA,ylab=NA){
  
  for(i in 1:nrow(mat)){
    for (j in 1:ncol(mat)){
      
      if(mat[i,j]=="TP"){
        mat[i,j]=1
      }
      
      else if(mat[i,j]=="TN"){
        mat[i,j]=2
      }
      
      else if(mat[i,j]=="FP"){
        mat[i,j]=3
      }
      
      else if(mat[i,j]=="FN"){
        mat[i,j]=4
      }
      
    }
  }
  
  nr <- dim(mat)[1]
  nc <- dim(mat)[2]
  
  mat=as.numeric(mat)
  
  mat=matrix(data=mat, ncol=nc, nrow=nr)
  
  #Request a square plot area:
  par(pty="s")
  #Plot an image using the row/column numbers as the x/y variables:
  image(x=1:nr,y=1:nc,mat,col=heat.colors(12),main=main,xlab=NA,ylab=NA)
  
  
}

