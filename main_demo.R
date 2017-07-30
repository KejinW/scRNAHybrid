setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./R/method.R")
source("./configure.R")

# Then transform the data in the format
# the user will always have a txt or csv file
# input file should be a string of the path to the file
# format should be “txt” or “csv” for defining the delimiter used in the input filed
input<-dataTrans("./Input/melanoma_data_with_labels.csv",Lab=TRUE)

# analyse the data with your hybrid method
# the configuration file should be the one that has all the parameters for the 
# different clustering methods
h<-Hybrid(input$data)

plotHybrid(h,PDF="plotHybrid.pdf")

# Plot anything you find interesting
h<-TPFNMatrix(h,input$label,plot=TRUE,PDF="plotTPFN.pdf")

h<-ARI(h,unlist(input$label),plot=TRUE,PDF="plotARI.pdf")

h<-ConsistentCluster(h,plot=TRUE,PDF="plotCCM.pdf")





