library(ggplot2)
library(dplyr)
library(DESeq2)
setwd('~/Desktop/HMH/')
theme_set(theme_bw())


###################################################################
res.group <- data.frame(matrix(nrow = nrow(res), ncol = 5), row.names = rownames(res))
colnames(res.group) <- c('nai_te','nai_tc',
                         'wt_te_tc','wt_ko_te','wt_ko_tc')
res.group %>% head()

inc<- res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5) %>% rownames()
dec <-res %>% filter(padj <= 0.05) %>% filter(log2FoldChange <= -1.5) %>% rownames()
res.group[,'nai_te'] <- 'same'
res.group[inc,'nai_te'] <- 'incr'
res.group[dec,'nai_te'] <- 'decr'

inc<- res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5) %>% rownames()
dec <-res %>% filter(padj <= 0.05) %>% filter(log2FoldChange <= -1.5) %>% rownames()
res.group[,'nai_tc'] <- 'same'
res.group[inc,'nai_tc'] <- 'incr'
res.group[dec,'nai_tc'] <- 'decr'
#########################################

inc<- res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5) %>% rownames()
dec <-res %>% filter(padj <= 0.05) %>% filter(log2FoldChange <= -1.5) %>% rownames()
res.group[,'wt_te_tc'] <- 'same'
res.group[inc,'wt_te_tc'] <- 'incr'
res.group[dec,'wt_te_tc'] <- 'decr'


inc<- res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5) %>% rownames()
dec <-res %>% filter(padj <= 0.05) %>% filter(log2FoldChange <= -1.5) %>% rownames()
res.group[,'wt_ko_te'] <- 'same'
res.group[inc,'wt_ko_te'] <- 'incr'
res.group[dec,'wt_ko_te'] <- 'decr'


inc<- res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5) %>% rownames()
dec <-res %>% filter(padj <= 0.05) %>% filter(log2FoldChange <= -1.5) %>% rownames()
res.group[,'wt_ko_tc'] <- 'same'
res.group[inc,'wt_ko_tc'] <- 'incr'
res.group[dec,'wt_ko_tc'] <- 'decr'


res.group %>% head()
count.mtx %>% head()

write.csv(res.group, 'rds/ATAC_seq/atac_seq.readvalue.and.normalized/res_csv/res.group.2021.11.16.csv')
res.group <- read.csv('rds/ATAC_seq/atac_seq.readvalue.and.normalized/res_csv/res.group.2021.11.16.csv', row.names = 1)
res.group %>% dim()
res.group %>% head()
## select groups
res.group <- res.group[,grep('wt',colnames(res.group))]

g1 <- res.group %>% filter(wt_te_tc !='same') %>% rownames()
g2 <- res.group %>% filter(wt_ko_te !='same') %>% rownames()
g3 <- res.group %>% filter(wt_ko_tc !='same') %>% rownames()

union(union(g1,g2),g3) %>% length()
sig.at.least.one <- union(union(g1,g2),g3)


norm.peak[sig.at.least.one,] %>% head()
norm.peak[sig.at.least.one,] %>% dim()

write.csv(norm.peak[sig.at.least.one,], 
          'rds/ATAC_seq/atac_seq.readvalue.and.normalized/atac_seq.norm.peak.sig.at.least.one.2021.11.16.csv')

norm.peak %>% head()


sig.norm.peak <- norm.peak[sig.at.least.one,c(6:16)]
sig.norm.peak <- norm.peak[sig.at.least.one,c(6:7,9:13,15:16)]
sig.norm.peak %>% head()
sig.norm.peak %>% dim()
input.data <- sig.norm.peak
## zscore function 
zscore <- function(input.data = input.data){
  input.data.rowsums <- rowSums(input.data)
  input.data.mean <- rowMeans(input.data)
  input.data.sd <- matrixStats::rowSds(as.matrix(input.data))
  names(input.data.sd) <- rownames(input.data)
  zscore <- (input.data-input.data.mean)/input.data.sd
  return(zscore)
}



fit <- kmeans(zscore(input.data), centers = 7, nstart = 10)
fit.cluster <- fit$cluster %>% data.frame()
colnames(fit.cluster) <- 'cluster'
fit.cluster <- fit.cluster %>% arrange(cluster)
fit.cluster$cluster <- factor(fit.cluster$cluster)
fit.cluster %>% head()
df.anno <- fit.cluster

fit.cluster$cluster %>% table()
write.csv(fit.cluster, 'rds/ATAC_seq/atac_seq.readvalue.and.normalized/sample9.fit.cluster7.csv')

pheatmap::pheatmap(zscore(log2(input.data[rownames(fit.cluster),]+1)), cluster_cols = F, 
                   cluster_rows = F, show_rownames = F, annotation_row = df.anno,
                   col = colorRampPalette(c("navy", "white", "red"))(1000),
                   fontsize_row = 6,fontsize_col = 10)

pheatmap::pheatmap(zscore(input.data[rownames(fit.cluster),]), cluster_cols = F, 
                   cluster_rows = F, show_rownames = F,annotation_row = df.anno,
                   col = colorRampPalette(c("navy", "white", "red"))(1000),
                   fontsize_row = 6,fontsize_col = 10)



sig.norm.tmp <- read.csv('rds/ATAC_seq/atac_seq.readvalue.and.normalized/tmp.csv', header = T)
rownames(sig.norm.tmp) <- paste0('peak_', rownames(sig.norm.tmp))
sig.norm.tmp %>% dim()

sig.norm.tmp %>% head()
sig.norm.tmp <- sig.norm.tmp[sig.at.least.one,]
sig.norm.tmp[rownames(fit.cluster),'group'] %>% head()
sig.norm.tmp[rownames(fit.cluster),'group'] <- fit.cluster$cluster
write.csv(sig.norm.tmp[rownames(fit.cluster),], 'rds/ATAC_seq/atac_seq.readvalue.and.normalized/tmp2.csv')
