################################  DEGs ###################################
###########################################################################

#####################################
setwd('~/Desktop/HMH/')
library(dplyr)
library(ggplot2)
library(DESeq2)
library(edgeR)
theme <- theme_classic()
theme <- theme(panel.background = element_rect(fill = "white", colour = "black"))


######## DEG analysis using DESeq2
###### input data : count matrix

atac_raw_count <- readRDS('~/Desktop/HMH/rds/2021.09.20.atac_seq.combined.peaks.rds')
dim(atac_raw_count) ## 44167

## change the row names to peak row names
rownames(atac_raw_count) <- paste0('peak_',rownames(atac_raw_count))

### save the separated conditions (total 5) reads
naive <- atac_raw_count[,grepl('naive', colnames(atac_raw_count))]
wt_tem <- atac_raw_count[,grepl('wt_tem', colnames(atac_raw_count))]
wt_tcm <- atac_raw_count[,grepl('wt_tcm', colnames(atac_raw_count))]
ko_tem <- atac_raw_count[,grepl('ko_tem', colnames(atac_raw_count))]
ko_tcm <- atac_raw_count[,grepl('ko_tcm', colnames(atac_raw_count))]


#######################################################################
#######################################################################
### load libraries ###
library(dplyr)
library(ggplot2)
library(DESeq2)

## information sheet
## naive, wt_tem
count.mtx <- cbind(naive, wt_tem)
nofrow <- length(colnames(count.mtx))
info <- data.frame(matrix(nrow = nofrow, ncol = 3))
colnames(info) <- c('sample', 'type','mutation')
info$sample <- colnames(count.mtx)
info$type <- substr(info$sample,1,5)
info$mutation <- 'wt'
info

### convert character to factor in type information
info$type <- factor(info$type, levels = c('naive','wt_te')) ## wt_te/naive

## create dds object
dds <- DESeqDataSetFromMatrix(count.mtx, info, ~type)
dds <- DESeq(dds)
res <- results(dds)
res <- data.frame(res)

samples <- 'wt_tem_vs_naive'
write.csv(res[,c(2,5,6)], paste0('~/Desktop/HMH/rds/atac_seq.deseq.',samples,'.csv'))

## naive, wt_tcm
count.mtx <- cbind(naive, wt_tcm)
colnames(count.mtx)
nofrow <- length(colnames(count.mtx))
info <- data.frame(matrix(nrow = nofrow, ncol = 3))
colnames(info) <- c('sample', 'type','mutation')
info$sample <- colnames(count.mtx)

info$type <- substr(info$sample,1,5)
info$mutation <- 'wt'
info

### convert character to factor in type information
info$type <- factor(info$type, levels = c('naive','wt_tc'))


## create dds object
dds <- DESeqDataSetFromMatrix(count.mtx, info, ~type)
#dim(dds)
dds <- DESeq(dds)
res <- results(dds)
#dim(res)
res <- data.frame(res)
res[res$log2FoldChange > 2,][1:3,]
count.mtx[rownames(res[res$log2FoldChange > 2,][1:3,]),]

samples <- 'wt_tcm_vs_naive'
write.csv(res[,c(2,5,6)], paste0('~/Desktop/HMH/rds/atac_seq.deseq.',samples,'.csv'))
dim(res)

