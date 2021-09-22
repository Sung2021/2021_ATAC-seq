### This is a chunk of codes for 
### 2021 ATAC-seq peak Differential analysis
### 
#####################################
setwd('~/Desktop/HMH/')
library(dplyr)
library(ggplot2)
library(DESeq2)
library(edgeR)
library(GenomicRanges)
theme <- theme_classic()
### import peak count 
raw.peaks <- readRDS('rds/2021.09.20.atac_seq.combined.peaks.rds')

### count matrix
### info sheet
count.mtx <- raw.peaks[,c(4:length(colnames(raw.peaks)))]
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- colnames(count.mtx)
info$cell_type <- substr(info$sample, 4,6)
info$mutation <- substr(info$sample, 1,2)
info
info[c(1:2),2] <- 'naive'
info[c(1:2),3] <- 'naive'

info$sample <- factor(info$sample, levels = c('naive1','naive2',
                                              'wt_tem1','wt_tem2','wt_tem3',
                                              'wt_tcm1','wt_tcm2','wt_tcm3',
                                              'ko_tem1','ko_tem2',
                                              'ko_tcm1','ko_tcm2','ko_tcm3'))
dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dim(dds)
dds <- DESeq(dds)
vsd <- vst(dds,blind=TRUE)
plotPCA(vsd,intgroup= 'sample')
## PCA plot in ggplot
plot_data1 <- plotPCA(vsd,intgroup= 'sample',returnData=TRUE)
plot_data1$sample <- factor(plot_data1$sample, levels = c('naive1','naive2',
                                                    'wt_tem1','wt_tem2','wt_tem3',
                                                    'wt_tcm1','wt_tcm2','wt_tcm3',
                                                    'ko_tem1','ko_tem2',
                                                    'ko_tcm1','ko_tcm2','ko_tcm3'))
plot_data1$condition <- substr(plot_data1$sample,1,5)
plot_data1$condition <- factor(plot_data1$condition, levels = c('naive',
                                                                'wt_te','wt_tc',
                                                                'ko_te','ko_tc'))
theme <- theme(panel.background = element_rect(fill = "white", colour = "black"))

ggplot(plot_data1, aes(x = PC1,y=PC2, col=condition,label=condition)) + 
  geom_point(size=4, shape=factor(plot_data1$condition)) + theme + xlab('PC1: 70% variance') + ylab('PC2: 26% variance')

### save dds file for future analysis
saveRDS(dds, 'rds/2021.09.22.ATAC_seq_all.rds')
