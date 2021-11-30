library(ggplot2)
library(dplyr)
library(DESeq2)
setwd('~/Desktop/HMH/')
theme_set(theme_bw())

norm.peak <- read.csv('rds/ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/atac_seq.21.11.29.read.norm.by.colmeans.with.naive.csv', 
                      row.names = 1)
norm.peak %>% head()
norm.peak %>% dim()
norm.peak.log2fc <- norm.peak

####################################################################
read.value <- read.csv('rds/ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/atac_seq.21.11.29.read.raw.with.naive.csv',
                       row.names = 5)
read.value %>% head()
read.value %>% dim()

#####################################################################

count.mtx <- read.value[,grep('naive|wt_te', colnames(read.value))]
group1 <- 'naive'
group2 <- 'wt_te'
count.mtx %>% head()

### information sheet
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- colnames(count.mtx)
info$cell_type <- substr(colnames(count.mtx),1,5)
info$mutation <- substr(colnames(count.mtx),1,2)
info$cell_type <- factor(info$cell_type, levels = c(group1, group2))
info

dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dds <- DESeq(dds)
res <- results(dds)
res %>% dim()
res <- data.frame(res)
res <- res[,c(2,5:6)]
colnames(res) <- paste0(paste0(group1,'_',group2), '_',colnames(res))
norm.peak.log2fc <- cbind(norm.peak.log2fc, res)
norm.peak.log2fc %>% head()


#######################################################################
count.mtx <- read.value[,grep('naive|wt_tc', colnames(read.value))]
group1 <- 'naive'
group2 <- 'wt_tc'
count.mtx %>% head()

### information sheet
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- colnames(count.mtx)
info$cell_type <- substr(colnames(count.mtx),1,5)
info$mutation <- substr(colnames(count.mtx),1,2)
info$cell_type <- factor(info$cell_type, levels = c(group1, group2))
info

dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dds <- DESeq(dds)
res <- results(dds)
res %>% dim()
res <- data.frame(res)
res <- res[,c(2,5:6)]
colnames(res) <- paste0(paste0(group1,'_',group2), '_',colnames(res))
norm.peak.log2fc <- cbind(norm.peak.log2fc, res)
norm.peak.log2fc %>% head()

#######################################################################
count.mtx <- read.value[,grep('wt_te|wt_tc', colnames(read.value))]
group1 <- 'wt_te'
group2 <- 'wt_tc'
count.mtx %>% head()

### information sheet
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- colnames(count.mtx)
info$cell_type <- substr(colnames(count.mtx),1,5)
info$mutation <- substr(colnames(count.mtx),1,2)
info$cell_type <- factor(info$cell_type, levels = c(group1, group2))
info

dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dds <- DESeq(dds)
res <- results(dds)
res %>% dim()
res <- data.frame(res)
res <- res[,c(2,5:6)]
colnames(res) <- paste0(paste0(group1,'_',group2), '_',colnames(res))
norm.peak.log2fc <- cbind(norm.peak.log2fc, res)
norm.peak.log2fc %>% head()


#######################################################################
count.mtx <- read.value[,grep('wt_te|ko_te', colnames(read.value))]
group1 <- 'wt_te'
group2 <- 'ko_te'
count.mtx %>% head()

### information sheet
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- colnames(count.mtx)
info$cell_type <- substr(colnames(count.mtx),1,5)
info$mutation <- substr(colnames(count.mtx),1,2)
info$cell_type <- factor(info$cell_type, levels = c(group1, group2))
info

dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dds <- DESeq(dds)
res <- results(dds)
res %>% dim()
res <- data.frame(res)
res <- res[,c(2,5:6)]
colnames(res) <- paste0(paste0(group1,'_',group2), '_',colnames(res))
norm.peak.log2fc <- cbind(norm.peak.log2fc, res)
norm.peak.log2fc %>% head()
res %>% head()
count.mtx %>% head()

#######################################################################
count.mtx <- read.value[,grep('wt_tc|ko_tc', colnames(read.value))]
group1 <- 'wt_tc'
group2 <- 'ko_tc'
count.mtx %>% head()

### information sheet
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- colnames(count.mtx)
info$cell_type <- substr(colnames(count.mtx),1,5)
info$mutation <- substr(colnames(count.mtx),1,2)
info$cell_type <- factor(info$cell_type, levels = c(group1, group2))
info

dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dds <- DESeq(dds)
res <- results(dds)
res %>% dim()
res <- data.frame(res)
res <- res[,c(2,5:6)]
colnames(res) <- paste0(paste0(group1,'_',group2), '_',colnames(res))
norm.peak.log2fc <- cbind(norm.peak.log2fc, res)
norm.peak.log2fc <- cbind(read.value[,c(2:4)],norm.peak.log2fc)
norm.peak.log2fc %>% head()

write.csv(norm.peak.log2fc,
          'rds/ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/atac_seq.21.11.30.norm.log2fc.all.csv')

