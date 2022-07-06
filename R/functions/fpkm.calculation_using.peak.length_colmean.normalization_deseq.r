

### This is a chunk of codes for 
### ATAC-seq analysis
#####################################
setwd('~/Desktop/HMH/rds')

library(ggplot2)
library(dplyr)
library(DESeq2)
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

####################################################################
## 2022.07.05
## raw reads : zero expression filtered, blacklist filtered
df <- read.csv('ATAC_seq/Icos.2022.06/Icos.ATAC_seq.all.read.raw.only.bl.filtered.22.06.24.csv',
               row.names = 1)
df[1:3,]
df %>% dim()


count.mtx <- df[,c(4:8)]
se <- SummarizedExperiment(as.matrix(count.mtx), 
                           colData=DataFrame(sample=1:ncol(count.mtx)))
dds <- DESeqDataSet(se, ~ 1)
rowRanges(dds) <- makeGRangesFromDataFrame(df[,c(1:3)])
fpkm <- fpkm(dds)
fpkm %>% dim() ## 52241
fpkm[1:3,]
colMeans(fpkm)

## colmean normalization
input.data <- fpkm
input.data[1:3,]
col.means <- input.data %>% colMeans()
col.means 
tmp <- data.frame(matrix(nrow=nrow(input.data), ncol = ncol(input.data)))
colnames(tmp) <- colnames(input.data)
rownames(tmp) <- rownames(input.data)
tmp[1:3,]
for(i in 1:nrow(tmp)){
  tmp[i,] <- col.means
}
tmp[1:10,]
norm.peak <- cbind(df[,c(1:3)],(input.data/tmp))
norm.peak[1:3,]
norm.peak %>% dim()

fpkm.colmean.norm.peak <- norm.peak

fpkm.colmean.norm.peak %>% 
  write.csv('ATAC_seq/Icos.2022.06/fpkm.colmean.norm.22.07.06/Icos.ATAC.fpkm.colmean.norm.peak.22.07.06.csv')


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

## information sheet
count.input <- df[,c(4:ncol(df))]
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

res %>% ggplot(aes(log2FoldChange, -log10(padj))) + geom_point(alpha=0.5)+
  geom_vline(xintercept = c(-log2(2), log2(2)), color='red') +
  geom_hline(yintercept = -log10(0.05), color='blue') + theme_bw()

fc <- 1.5
res %>% filter(padj < 0.05) %>%  filter(log2FoldChange > log2(fc)) %>% dim()
res %>% filter(padj < 0.05) %>%  filter(log2FoldChange < -log2(fc)) %>% dim()
res %>% filter(pvalue < 0.05) %>%  filter(log2FoldChange > log2(fc)) %>% dim()
res %>% filter(pvalue < 0.05) %>%  filter(log2FoldChange < -log2(fc)) %>% dim()
