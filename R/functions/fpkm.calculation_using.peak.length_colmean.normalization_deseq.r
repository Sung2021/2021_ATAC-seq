## generate DESeq2 object (dds) from count.mtx (raw count)
## generate GRange info from peak length information
## replace rowRanges(dds) with GRange object from peak length information
## use fpkm(dds) function in DESeq to calculate fpkm
## I couldn't find the fpkm formula to use it manually

#####################################
setwd('~/Desktop/HMH/rds')

library(ggplot2)
library(dplyr)
library(DESeq2)
library(GenomicRanges)
library(rtracklayer)

####################################################################
## 2022.07.05
## raw reads : zero expression filtered, blacklist filtered
df <- read.csv('ATAC_seq/Icos.2022.06/Icos.ATAC_seq.all.read.raw.only.bl.filtered.22.06.24.csv',
               row.names = 1)
df[1:3,]
df %>% dim()

## fpkm calculation 
count.mtx <- df[,c(4:8)]
se <- SummarizedExperiment(as.matrix(count.mtx), 
                           colData=DataFrame(sample=1:ncol(count.mtx)))
dds <- DESeqDataSet(se, ~ 1)
rowRanges(dds) <- makeGRangesFromDataFrame(df[,c(1:3)])
fpkm <- fpkm(dds)
fpkm %>% dim() ## 52241
fpkm[1:3,]
colMeans(fpkm)

## fpkm-colmean normalization
input.data <- fpkm
col.means <- input.data %>% colMeans()
tmp <- data.frame(matrix(nrow=nrow(input.data), 
                         ncol = ncol(input.data)))
colnames(tmp) <- colnames(input.data)
rownames(tmp) <- rownames(input.data)
for(i in 1:nrow(tmp)){
  tmp[i,] <- col.means
}
norm.peak <- cbind(df[,c(1:3)],(input.data/tmp))
norm.peak %>% dim()
fpkm.colmean.norm <- norm.peak
fpkm.colmean.norm %>% write.csv('ATAC_seq/Icos.2022.06/Icos.ATAC_seq.fpkm.colmean.norm.22.07.08.csv')

##########################################################################

#prcomp(t(input.data), scale. = T)$x
norm.mtx <- fpkm.colmean.norm.peak[,c(4:8)]
norm.mtx %>% dim()
input.data <- norm.mtx
pca = as.data.frame(prcomp(t(input.data), scale. = T)$x)
pca %>% str()
pca
PCA <- prcomp(t(input.data),scale.=T)
PCA <- (PCA$sdev/sum(PCA$sdev))[1:2]
PCA_percentage <- paste(round(PCA,4)*10^2,'%',sep = '')

ggplot(pca, aes(PC1,PC2,label=rownames(pca), 
                color= substr(rownames(pca),1,2))) + geom_point(size=4) +
  xlab(paste0('PC1: ',PCA_percentage[1])) +
  ylab(paste0('PC2: ',PCA_percentage[2])) + 
  theme(legend.position = 'none') +
  theme_bw()

ggplot(pca, aes(PC1,PC2,label=rownames(pca), 
                color= substr(rownames(pca),1,2))) + geom_point(size=4, alpha=0.5) +
  xlab(paste0('PC1: ',PCA_percentage[1])) +
  ylab(paste0('PC2: ',PCA_percentage[2])) + 
  theme(legend.position = 'none') + geom_text() +
  theme_bw()


#####################################################################

## DEPeak 
## information sheet
count.input <- df[,c(4:ncol(df))] ## raw integer input
info <- data.frame(matrix(nrow = ncol(count.input), ncol = 2))
colnames(info) <- c('sample', 'cond')
info$sample <- colnames(count.input)
info$cond <- substr(info$sample, 1,2)
info$cond <- factor(info$cond, levels = c('WT','KO'))
info 

dds <- DESeqDataSetFromMatrix(count.input, info, ~ cond)
dim(dds)
dds <- DESeq(dds)
res <- results(dds)
dim(res)
res <- data.frame(res)
res[1:3,]

res %>% ggplot(aes(log2FoldChange, -log10(padj))) + geom_point(size=1, alpha=0.5) +
  geom_vline(xintercept = c(-log2(2), log2(2)), color='red') +
  geom_hline(yintercept = -log10(0.05), color='blue') + theme_classic()

fc <- 2
pval <- 0.05
qval <- 0.05
res$color <- 'NA'
red <- res %>% filter(pvalue < pval) %>%  filter(log2FoldChange > log2(fc)) %>% rownames()
blue <- res %>% filter(pvalue < pval) %>%  filter(log2FoldChange < -log2(fc)) %>% rownames()
res[red,]$color <- 'red'
res[blue,]$color <- 'blue'
res$size <- 'small'
res[(res %>% filter(color %in% c('red','blue')) %>% rownames()),]$size <- 'normal'
res %>% ggplot(aes(log2FoldChange, -log10(pvalue),
                   color=color,size=size)) + geom_point(alpha=0.5) +
  geom_vline(xintercept = c(-log2(fc), log2(fc)), color='red') +
  geom_hline(yintercept = -log10(pval), color='blue') + 
  scale_color_manual(values = c('blue','grey','red')) +
  scale_size_manual(values = c(1.5,0.5)) +
  theme_classic() +  theme(legend.position = 'none')


fc <- 1.5
res %>% filter(padj < 0.05) %>%  filter(log2FoldChange > log2(fc)) %>% dim()
res %>% filter(padj < 0.05) %>%  filter(log2FoldChange < -log2(fc)) %>% dim()
res %>% filter(pvalue < 0.05) %>%  filter(log2FoldChange > log2(fc)) %>% dim()
res %>% filter(pvalue < 0.05) %>%  filter(log2FoldChange < -log2(fc)) %>% dim()


