
library(SummarizedExperiment)
library(devtools)
library(caret)
library(mclust)
library(Rtsne)

##-----------------------------SAVER------------------------------- 
library(SAVER)
packageVersion("SAVER")
library(data.table)

# Need to remove first 10 rows of metadata
data.path <- "./data/yan_data.csv"
label.path <- "./data/yan_celldata.csv"
#class.label <- fread(label.path, header = F)
class.label <- as.matrix(fread(label.path, header = T))
class.label <- class.label[,-1]
ncls = unique(class.label)
raw.data <- fread(data.path, header = T)
cortex <- as.matrix(raw.data[, -1])
#cortex <- apply(raw.data[-1,],2,as.numeric)
cellnames <- colnames(raw.data)
colnames(cortex) <- cellnames[-1]
dim(cortex)     # dim = genes*cells
View(cortex)

#run SAVER
cortex.saver <- saver(cortex, ncores = 12)
str(cortex.saver)
#ARI
adjustedRandIndex(kmeans(Rtsne(t(as.matrix(cortex.saver[["estimate"]])))$Y, centers = 3)$cluster, class.label)
adjustedRandIndex(kmeans(Rtsne(t(as.matrix(cortex.saver[["se"]])))$Y, centers = 11)$cluster, class.label)
##---------------------------DrImpute-----------------------------
library(DrImpute)

data.path <- "./data/yan_data.csv"
#label.path <- "D:/bigdata/scIGANsaaa/all data/pollen/pollen_celldata.csv"
raw.data <- fread(data.path, header = T)
X <- as.matrix(raw.data[, -1]) #gene*cell

#X <- preprocess(X, min.expressed.gene = 0)
X.log <- log(X + 1)
#class.label <- fread(label.path, header = F)
#class.label <- as.matrix(fread(label.path, header = F))
#run DrImpute
set.seed(1)
X.imp <- DrImpute(X.log)

#ARI

# before imputation
adjustedRandIndex(kmeans(Rtsne(t(as.matrix(X.log)))$Y, centers = 11)$cluster, class.label)

# after imputation
adjustedRandIndex(kmeans(Rtsne(t(as.matrix(X.imp)))$Y, centers = 11)$cluster, class.label)

#------------------------DeepImpute-------------------------------
deep_path <- "..."
deep <- fread(deep_path, header = T)
deep = deep[,-1]
adjustedRandIndex(kmeans(Rtsne(as.matrix(deep),perplexity=10)$Y, centers = 5)$cluster, class.label)

#--------------------------MAGIC------------------------------
magic_path <- "..."
magic <- fread(magic_path, header = T)
magic = magic[,-1]
adjustedRandIndex(kmeans(Rtsne(as.matrix(magic))$Y, centers = 3)$cluster, class.label)

#---------------------------DCA--------------------------------
DCA_path <- "..."
DCA <- fread(DCA_path, header = T)
DCA <- DCA[,-1]
adjustedRandIndex(kmeans(Rtsne(as.matrix(DCA))$Y, centers = 3)$cluster, class.label)

#---------------------------scIGANs--------------------------------
scIGANs_path <- "..."
scIGANs <- fread(scIGANs_path, header = T)
adjustedRandIndex(kmeans(Rtsne(as.matrix(scIGANs))$Y, centers = 3)$cluster, class.label)
#---------------------------scGCT---------------------------------
scGCT_path <- "..."
scGCT <- fread(scGCT_path, header = T)
adjustedRandIndex(kmeans(Rtsne(as.matrix(scGCT))$Y, centers = 5)$cluster, class.label)

