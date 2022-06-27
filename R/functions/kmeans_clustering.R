library(ggplot2)
library(dplyr)
library(DESeq2)
library(reshape)
setwd('~/Desktop/HMH/rds')
theme_set(theme_bw())

## zscore function 
zscore <- function(input.data = input.data){
  input.data.rowsums <- rowSums(input.data)
  input.data.mean <- rowMeans(input.data)
  input.data.sd <- matrixStats::rowSds(as.matrix(input.data))
  names(input.data.sd) <- rownames(input.data)
  zscore <- (input.data-input.data.mean)/input.data.sd
  return(zscore)
}

sig.peak.norm <- read.csv('ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/atac_seq.21.11.29.sig.peaks.norm.csv',
                          row.names = 1)
sig.peak.norm %>% head()

fit <- kmeans(zscore(input.data), centers = 8, nstart = 10)
fit.cluster <- fit$cluster %>% data.frame()
colnames(fit.cluster) <- 'cluster'
fit.cluster <- fit.cluster %>% arrange(cluster)
fit.cluster$cluster <- factor(fit.cluster$cluster)
fit.cluster %>% head()
df.anno <- fit.cluster
fit.cluster$cluster %>% table()


pheatmap::pheatmap(zscore(log2(input.data[rownames(fit.cluster),]+1)), cluster_cols = F, 
                   cluster_rows = F, show_rownames = F, annotation_row = df.anno,
                   col = colorRampPalette(c("navy", "white", "red"))(1000),
                   fontsize_row = 6,fontsize_col = 10)

