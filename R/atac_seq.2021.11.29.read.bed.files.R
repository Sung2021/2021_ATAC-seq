library(ggplot2)
library(dplyr)
library(DESeq2)
library(GenomicRanges)
setwd('~/Desktop/HMH/')
theme_set(theme_bw())
## read input merged bed from macs2 call
dir <- 'rds/ATAC_seq/atac_seq.bed.21.11.29'
list.files(dir, full.names = T)

### when reading coverage bed
for(i in 1:length(list.files(dir, full.names = T))){
  bed <- read.table(list.files(dir, full.names = T)[i])
  colnames(bed) <- c('chr','start','end',
                     'peak','score','strand',
                     'signalvalue','pvalue','qvalue','peak.score',
                     'read','frg_len','ref_len','coverage')
  bed$sample <- substr(list.files(dir, full.names = T)[i], 36,42)
  assign(substr(list.files(dir, full.names = T)[i], 36,42),bed)
}
ls(pattern = 'wt_')
ls(pattern = 'ko_')
ls(pattern = 'na_')
wt_tcm1 %>% head()
wt_tcm1 %>% dim()

wt <- list(wt_tcm1, wt_tcm2, wt_tcm3,wt_tem1, wt_tem2)
ko <- list(ko_tcm1, ko_tcm2, ko_tcm3,ko_tem1, ko_tem2)
na <- list(naive1., naive2.)

read.raw <- wt[[1]][,c(1:4)]
for(i in 1:5){
  read.raw[,paste0(wt[[i]][,'sample'][1],'_read')] <- wt[[i]][11]
}
read.raw %>% head()
for(i in 1:5){
  read.raw[,paste0(ko[[i]][,'sample'][1],'_read')] <- ko[[i]][11]
}
read.raw %>% head()
for(i in 1:2){
  read.raw[,paste0(na[[i]][,'sample'][1],'_read')] <- na[[i]][11]
}
read.raw %>% head()
read.raw[,c(1:4,15:16,8:9,5:7,13:14,10:12)] %>% head()
read.raw <- read.raw[,c(1:4,15:16,8:9,5:7,13:14,10:12)]
write.csv(read.raw, 
          'rds/ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/atac_seq.21.11.29.read.raw.with.naive.csv')

read.raw[,c(5:ncol(read.raw))] %>% head()
table(read.raw[,c(5:ncol(read.raw))] %>% rowMeans() == 0 )
norm.factor <- read.raw[,c(5:ncol(read.raw))] %>% colMeans()
norm.colmeans <- data.frame(matrix(nrow = nrow(read.raw), ncol=12)) ## ncol=10 without naive
for(i in 1:12){
  norm.colmeans[,i] <- norm.factor[i]
}
norm.colmeans %>% head()
norm.peak.atac <- read.raw[,c(5:ncol(read.raw))]/norm.colmeans
rownames(norm.peak.atac) <- read.raw$peak
norm.peak.atac %>% head()
write.csv(norm.peak.atac, 
          'rds/ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/atac_seq.21.11.29.read.norm.by.colmeans.with.naive.csv')

