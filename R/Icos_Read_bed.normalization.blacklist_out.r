### Icos WT/KO 72hr ATAC-seq data analysis
#####################################
setwd('~/Desktop/HMH/rds')

library(ggplot2)
library(dplyr)
library(DESeq2)
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

theme_set(theme_bw())
## read input merged bed from macs2 call
dir <- 'ATAC_seq/Icos.2022.06/bed.cov'
list.files(dir, full.names = T)

### when reading coverage bed
samples <- c(paste0('Icos_','KO',1:3),
             paste0('Icos_','WT',1:3))
for(i in 1:length(list.files(dir, full.names = T))){
  bed <- read.table(list.files(dir, full.names = T)[i])
  colnames(bed) <- c('chr','start','end',
                     'peak','score','strand',
                     'signalvalue','pvalue','qvalue','peak.score',
                     'read','frg_len','ref_len','coverage')
  bed$sample <- samples[i] 
  assign(samples[i],bed)
}
ls(pattern = 'Icos_')

## assign new name
## sample, chr, start, end
samples.bed <- c(paste0('Icos_','WT',1:3, '.bed'),
             paste0('Icos_','KO',1:3, '.bed'))
bed.list <- list(Icos_WT1,Icos_WT2, Icos_WT3,
                Icos_KO1,Icos_KO2, Icos_KO3)
for(i in 1:6){
  bed.tmp <- bed.list[[i]]
  rownames(bed.tmp) <- paste0(bed.tmp$sample, '.',bed.tmp$chr,'_',bed.tmp$start,'_', bed.tmp$end)
  assign(samples.bed[i],bed.tmp)
}
ls(pattern = 'Icos')

Icos.all.bed <- rbind(Icos_WT1.bed,
                      Icos_WT2.bed,
                      Icos_WT3.bed,
                      Icos_KO1.bed,
                      Icos_KO2.bed,
                      Icos_KO3.bed)
Icos.all.bed %>% dim()
## save the bed files as csv : bed file including all sample data in on format.
Icos.all.bed %>% write.csv('ATAC_seq/Icos.2022.06/Icos.ATAC_seq.bed.all.22.06.09.csv')

levels(as.factor(Icos.all.bed$sample))
df1 <- Icos.all.bed %>% filter(sample == levels(as.factor(Icos.all.bed$sample))[1]) %>% select(peak, read)
df2 <- Icos.all.bed %>% filter(sample == levels(as.factor(Icos.all.bed$sample))[2]) %>% select(peak, read)
df3 <- Icos.all.bed %>% filter(sample == levels(as.factor(Icos.all.bed$sample))[3]) %>% select(peak, read)
df4 <- Icos.all.bed %>% filter(sample == levels(as.factor(Icos.all.bed$sample))[4]) %>% select(peak, read)
df5 <- Icos.all.bed %>% filter(sample == levels(as.factor(Icos.all.bed$sample))[5]) %>% select(peak, read)
df6 <- Icos.all.bed %>% filter(sample == levels(as.factor(Icos.all.bed$sample))[6]) %>% select(peak, read)

df <- cbind(df1,df2[,2],df3[2],df4[,2],df5[,2],df6[,2])
rownames(df) <- paste0(Icos_KO1$chr,'.',Icos_KO1$start,'_',Icos_KO1$end)
df <- df[,-1] # removal of merge peak name
colnames(df) <- levels(as.factor(Icos.all.bed$sample))
df <- cbind(Icos_WT1[,c(1:3)],df) ## genomic information
## save the raw read file as the csv format
df %>% write.csv('ATAC_seq/Icos.2022.06/Icos.ATAC_seq.all.read.raw.only.22.06.09.csv')

## blacklist removal
## blacklist
bl <- read.table('~/Desktop/ref/blacklist/mm10.blacklist.bed')
colnames(bl) <- c('chr','start','end')
bl %>% dim() ## check how many blacklist regions in the file
bl$chr %>% table() ## mm10 blacklist only contains regions from chr1-19

## sebsetting out the bl regions from the pool
## using GRange format
df.gr <- makeGRangesFromDataFrame(df[,c(1:3)])
bl.gr <- makeGRangesFromDataFrame(bl[,c(1:3)])
##
setdiff(df.gr, bl.gr)  # check how many regions remaining 
df %>% dim() # 

df.filtered.gr <- setdiff(df.gr, bl.gr) %>% as.data.frame()
rownames(df.filtered.gr) <- paste0(df.filtered.gr$seqnames,'.',
                       df.filtered.gr$start,'_',df.filtered.gr$end) ## new name
df.filtered.gr[1:3,]
df <- df %>% filter(rownames(df) %in% rownames(df.filtered.gr)) ## 75065 regions
df <- df %>% filter(chr %in% c(paste0('chr',1:19),
                               paste0('chr', c('X','Y')))) ## 75063 regions
## keep sex chromosome region : donor and recipient both are male.
df %>% dim()
## save data with raw, bl filtered
df %>% write.csv('ATAC_seq/Icos.2022.06/Icos.ATAC_seq.all.read.raw.only.bl.filtered.22.06.09.csv')
df$chr %>% table()


## colmean normalization
df[1:3,]
col.means <- df[,c(4:9)] %>% colMeans()
col.means 
tmp <- data.frame(matrix(nrow=nrow(df), ncol = ncol(df)-3))
colnames(tmp) <- colnames(df[,4:9])
rownames(tmp) <- rownames(df)
for(i in 1:nrow(tmp)){
  tmp[i,] <- col.means[i]
}
norm.peak <- cbind(df[,1:3],df[,4:9]/col.means)
norm.peak[1:3,]
## save data with normalized value
norm.peak %>% 
  write.csv('ATAC_seq/Icos.2022.06/Icos.ATAC_seq.all.read.raw.only.bl.filtered.colmean.norm.22.06.09.csv')
norm.peak %>% dim()


### analyze data
norm.peak <- read.csv('ATAC_seq/Icos.2022.06/Icos.ATAC_seq.all.read.raw.only.bl.filtered.colmean.norm.22.06.09.csv',
                      row.names = 1)
norm.peak[1:3,]
