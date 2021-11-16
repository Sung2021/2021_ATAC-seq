library(ggplot2)
library(dplyr)
library(DESeq2)
setwd('~/Desktop/HMH/')
theme_set(theme_bw())

########## read the samples ########
dir <- 'rds/ATAC_seq/atac_seq.bed'
feature.list <- list.files(dir, pattern = 'bed$', full.names = T)
group <- 'naive'
group <- 'wt'
group <- 'ko'

feature.list[grep(group, feature.list)]

for(i in 1:6){
  group <- 'wt' #1 assign group
  sample.name <-substr(feature.list[grep(group, feature.list)][i], 27,33)
  bed <- read.table(feature.list[grep(group, feature.list)][i])
  colnames(bed) <- c('chr','start','end',
                     'peak','score','strand',
                     'signalvalue','pvalue','qvalue','peak.score',
                     'read','frg_len','ref_len','coverage')
  colnames(bed)[11] <- paste0(sample.name,'.read')
  assign(sample.name, bed)
}
wt_tcm1 %>% head()
wt_tem1 %>% head()

## ko has 5 samples
for(i in 1:5){
  group <- 'ko' #1 assign group
  sample.name <-substr(feature.list[grep(group, feature.list)][i], 27,33)
  bed <- read.table(feature.list[grep(group, feature.list)][i])
  colnames(bed) <- c('chr','start','end',
                     'peak','score','strand',
                     'signalvalue','pvalue','qvalue','peak.score',
                     'read','frg_len','ref_len','coverage')
  colnames(bed)[11] <- paste0(sample.name,'.read')
  assign(sample.name, bed)
}
ko_tcm1 %>% head()
ko_tem1 %>% head()

## naive has 2 samples
for(i in 1:2){
  group <- 'naive' #1 assign group
  sample.name <-substr(feature.list[grep(group, feature.list)][i], 27,33)
  bed <- read.table(feature.list[grep(group, feature.list)][i])
  colnames(bed) <- c('chr','start','end',
                     'peak','score','strand',
                     'signalvalue','pvalue','qvalue','peak.score',
                     'read','frg_len','ref_len','coverage')
  colnames(bed)[11] <- paste0(sample.name,'.read')
  assign(sample.name, bed)
}
naive01 %>% head()
naive02 %>% head()

############################################################################

tmp <- readRDS('rds/ATAC_seq/2021.09.20.atac_seq.combined.peaks.rds')
tmp %>% head()
read.value <- tmp



############ reference ##########


col.means <- read.value[,c(4:16)] %>% colMeans()
col.means.df <- data.frame(matrix(nrow = nrow(read.value), ncol=(16-4+1)))
colnames(col.means.df) <- colnames(read.value[,c(4:16)])
for(i in 1:(16-4+1)){
  col.means.df[,i] <- col.means[i]
}
col.means.df

normalization_factor <- col.means

norm.peak <- read.value[,c(4:16)]/col.means.df
norm.peak <- cbind(read.value[,c(1:3)], norm.peak)
saveRDS(norm.peak, 'rds/ATAC_seq/atac_seq.2021.11.16.normalized.by.colmean.rds')
write.csv(norm.peak, 'rds/ATAC_seq/atac_seq.2021.11.16.normalized.by.colmean.csv')

tmp <- readRDS('rds/ATAC_seq/atac_seq.readvalue.and.normalized/atac_seq.2021.11.16.normalized.by.colmean.rds')
tmp %>% dim()
tmp %>% head()
