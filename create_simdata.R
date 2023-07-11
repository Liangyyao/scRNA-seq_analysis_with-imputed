##---------------------模拟数据集-------------------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("splatter")
library(splatter)
library(arulesViz)
library(rpy2)
#install.packages("rpy2")
simulate <- function(nGroups=2, nGenes=225, batchCells=2000, dropout=3)
{
  if (nGroups > 1) method <- 'groups'
  else             method <- 'single'
  
  group.prob <- rep(1, nGroups) / nGroups
  
  # new splatter requires dropout.type
  if ('dropout.type' %in% slotNames(newSplatParams())) {
    if (dropout)
      dropout.type <- 'experiment'
    else
      dropout.type <- 'none'
    
    sim <- splatSimulate(group.prob=group.prob, nGenes=nGenes, batchCells=batchCells,
                         dropout.type=dropout.type, method=method,
                         seed=42, dropout.shape=-1, dropout.mid=dropout)
    
  } else {
    sim <- splatSimulate(group.prob=group.prob, nGenes=nGenes, batchCells=batchCells,
                         dropout.present=!dropout, method=method,
                         seed=42, dropout.shape=-1, dropout.mid=dropout)        
  }
  
  counts     <- as.data.frame(t(counts(sim)))
  truecounts <- as.data.frame(t(assays(sim)$TrueCounts))
  
  dropout    <- assays(sim)$Dropout
  #mode(dropout) <- 'integer'
  
  cellinfo   <- as.data.frame(colData(sim))
  geneinfo   <- as.data.frame(rowData(sim))
  
  list(counts=counts,
       cellinfo=cellinfo,
       geneinfo=geneinfo,
       truecounts=truecounts)
}

sim <- simulate()

counts <- sim$counts
geneinfo <- sim$geneinfo
cellinfo <- sim$cellinfo
truecounts <- sim$truecounts
write.csv(counts,"D:/bigdata/scIGANsaaa/sim2000.csv")
write.csv(geneinfo,"D:/bigdata/scIGANsaaa/sim2000_g.csv")
write.csv(cellinfo,"D:/bigdata/scIGANsaaa/sim2000_c.csv")
write.csv(truecounts,"D:/bigdata/scIGANsaaa/sim2000t.csv")


