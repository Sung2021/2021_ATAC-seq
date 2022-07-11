### This is a chunk of codes for 
### ATAC-seq analysis
#####################################
setwd('~/Desktop/HMH/rds')

library(ggplot2)
library(dplyr)
library(DESeq2)
library(GenomicRanges)
library(rtracklayer)

gene.info <- read.csv('RNA_seq/gene_genomic_info.csv',row.names = 1)
tss.info <- read.csv('ATAC_seq/gene.tss.info.22.06.30.csv',row.names = 1)
dis.info <- read.csv('ATAC_seq/gene.distal.info.22.06.30.csv', row.names = 1)

gene.info <- gene.info[,c(2:5,7)]
gene.info[1:3,]
gene.info$location_name <- paste0(gene.info$seqnames,'_',
                                  gene.info$start, '_',
                                  gene.info$end,'_',
                                  gene.info$SYMBOL,'_',
                                  'gene_body')
tss.info <- tss.info[,c(1:4,6)]
tss.info[1:3,]
tss.info$location_name <- paste0(tss.info$seqnames,'_',
                                  tss.info$start, '_',
                                  tss.info$end,'_',
                                  tss.info$SYMBOL,'_',
                                  'promoter')
dis.info[1:3,]
tmp <- dis.info[,c(1:2,3:4)]
colnames(tmp)[3:4] <- c('start','end')
tmp$location_name <- paste0(tmp$seqnames,'_',
                            tmp$start, '_',
                            tmp$end,'_',
                            tmp$SYMBOL,'_',
                            'dis.up')
tmp[1:3,]
tmp.up <- tmp

tmp <- dis.info[,c(1:2,5:6)]
colnames(tmp)[3:4] <- c('start','end')
tmp$location_name <- paste0(tmp$seqnames,'_',
                            tmp$start, '_',
                            tmp$end,'_',
                            tmp$SYMBOL,'_',
                            'dis.up')
tmp[1:3,]
tmp.dn <- tmp

location.info <- rbind(gene.info,tss.info,tmp.up, tmp.dn)
location.info %>% dim()
location.info %>% write.csv('ATAC_seq/location.info.22.07.11.csv')

location.info <- read.csv('ATAC_seq/location.info.22.07.11.csv',
                          row.names = 1)
rownames(location.info) <- location.info$location_name
location.info <- location.info[,c(2:4)]
location.info[1:3,]
location.info %>% write.csv('ATAC_seq/location.info.22.07.11.csv')
location.info <- read.csv('ATAC_seq/location.info.22.07.11.csv',
                          row.names = 1)
location.info[1:3,]
