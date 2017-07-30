#clustering number
c=4

#Seurat
#pbmc <- Setup(pbmc, min.cells = 3, min.genes = 200, do.logNormalize = F, total.expr = 1e4, project = "10X_PBMC")
min.cells=3
min.genes=200

#pbmc <- RegressOut(pbmc, latent.vars = c("nUMI", "nGene"), do.scale = T, do.center = T)
do.scale = T
do.center = T

#pbmc <- MeanVarPlot(pbmc ,fxn.x = expMean, fxn.y = logVarDivMean, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)
x.low.cutoff = 0.0125
x.high.cutoff = 3
y.cutoff = 0.5

#pbmc <- JackStraw(pbmc, num.replicate = 100, do.print = FALSE)
num.replicate = 100

#pbmc <- FindClusters(pbmc, pc.use = 1:10, resolution = 1, print.output = 0, save.SNN = T)
pc.use = 1:10
resolution = 1
save.SNN = F

#pbmc <- RunTSNE(pbmc, dims.use = 1:10, do.fast = T)
dims.use = 1:10


#SC3
#sceset <- calculateQCMetrics(sceset, nmads = 5, pct_feature_controls_threshold = 80)
nmads = 5
pct_feature_controls_threshold = 80

#sceset <- sc3(sceset, ks = ClusterNum, biology = TRUE,
#pct_dropout_min = 10,pct_dropout_max = 90,d_region_min = 0.04, d_region_max = 0.07,svm_max = 5000)
pct_dropout_min = 10
pct_dropout_max = 90
d_region_min = 0.04
d_region_max = 0.07
svm_max = 5000
kmeans_nstart = 50
kmeans_iter_max = 50

#CIDR
#sData <- determineDropoutCandidates(sDatamin1 = 3, min2 = 8,
#N = 2000, alpha = 0.1)
min1 = 0.1
min2 = 0.1
N = 2000
alpha = 0.1

#sData <- wThreshold(sData)
cutoff = 0.5

#SIMLR
#res_example1 = SIMLR(X=Data,c=ClusterNum,k = 10, normalize = FALSE, cores.ratio = 1)
large_scale=TRUE
k=10
cores.ratio = 1