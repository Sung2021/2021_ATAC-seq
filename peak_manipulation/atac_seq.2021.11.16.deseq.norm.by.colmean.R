library(ggplot2)
library(dplyr)
library(DESeq2)
setwd('~/Desktop/HMH/')
theme_set(theme_bw())

norm.peak <- readRDS('rds/ATAC_seq/atac_seq.readvalue.and.normalized/atac_seq.2021.11.16.normalized.by.colmean.rds')
norm.peak %>% head()
read.value <- readRDS('rds/ATAC_seq/atac_seq.readvalue.and.normalized/2021.09.20.atac_seq.combined.peaks.rds')
saveRDS(read.value,'rds/ATAC_seq/atac_seq.readvalue.and.normalized/2021.11.16.atac_seq.combined.peaks.rds')
rownames(read.value) <- paste0('peak_', rownames(read.value))
read.value %>% head()
read.value %>% dim()


colnames(read.value)

###################################################################
######## DESeq by groups #######
count.mtx <- read.value[,grep('naive|wt_te', colnames(read.value))]
group1 <- 'naive'
group2 <- 'wt_te'
count.mtx %>% head()


### information sheet
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- colnames(count.mtx)
info$cell_type <- substr(colnames(count.mtx),1,5)
info$mutation <- 'no'
info$cell_type <- factor(info$cell_type, levels = c(group1, group2))
info


dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dds <- DESeq(dds)
res <- results(dds)
res %>% dim()
res <- data.frame(res)
res %>% filter(padj <= 0.05) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 2 | log2FoldChange <= -2) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5 | log2FoldChange <= -1.5) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>% dim()
res %>% head()
res <- res[,c(2,5,6)]
write.csv(res,paste0('rds/ATAC_seq/atac_seq.readvalue.and.normalized/res/', group1,'_',group2, '.csv'))



#####################################################################
count.mtx <- read.value[,grep('naive|wt_tc', colnames(read.value))]
group1 <- 'naive'
group2 <- 'wt_tc'
count.mtx %>% head()


### information sheet
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- colnames(count.mtx)
info$cell_type <- substr(colnames(count.mtx),1,5)
info$mutation <- 'no'
info$cell_type <- factor(info$cell_type, levels = c(group1, group2))
info


dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dds <- DESeq(dds)
res <- results(dds)
res %>% dim()
res <- data.frame(res)
res %>% filter(padj <= 0.05) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 2 | log2FoldChange <= -2) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5 | log2FoldChange <= -1.5) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>% dim()
res %>% head()
res <- res[,c(2,5,6)]
write.csv(res,paste0('rds/ATAC_seq/atac_seq.readvalue.and.normalized/res/', group1,'_',group2, '.csv'))



#####################################################################
count.mtx <- read.value[,grep('wt_te|wt_tc', colnames(read.value))]
group1 <- 'wt_te'
group2 <- 'wt_tc'
count.mtx %>% head()


### information sheet
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- colnames(count.mtx)
info$cell_type <- substr(colnames(count.mtx),1,5)
info$mutation <- 'no'
info$cell_type <- factor(info$cell_type, levels = c(group1, group2))
info


dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dds <- DESeq(dds)
res <- results(dds)
res %>% dim()
res <- data.frame(res)
res %>% filter(padj <= 0.05) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 2 | log2FoldChange <= -2) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5 | log2FoldChange <= -1.5) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>% dim()
res %>% head()
res <- res[,c(2,5,6)]
write.csv(res,paste0('rds/ATAC_seq/atac_seq.readvalue.and.normalized/res/', group1,'_',group2, '.csv'))



#####################################################################
count.mtx <- read.value[,grep('wt_te|ko_te', colnames(read.value))]
group1 <- 'wt_te'
group2 <- 'ko_te'
count.mtx %>% head()
mutation <- factor(c(rep('wt',3), rep('ko',2)))

### information sheet
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- colnames(count.mtx)
info$cell_type <- substr(colnames(count.mtx),1,5)
info$mutation <- mutation
info$cell_type <- factor(info$cell_type, levels = c(group1, group2))

info


dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dds <- DESeq(dds)
res <- results(dds)
res %>% dim()
res <- data.frame(res)
res %>% filter(padj <= 0.05) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 2 | log2FoldChange <= -2) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5 | log2FoldChange <= -1.5) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>% dim()
res %>% head()
res <- res[,c(2,5,6)]
write.csv(res,paste0('rds/ATAC_seq/atac_seq.readvalue.and.normalized/res/', group1,'_',group2, '.csv'))


#####################################################################
count.mtx <- read.value[,grep('wt_tc|ko_tc', colnames(read.value))]
group1 <- 'wt_tc'
group2 <- 'ko_tc'
count.mtx %>% head()
mutation <- factor(c(rep('wt',3), rep('ko',2)))

### information sheet
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- colnames(count.mtx)
info$cell_type <- substr(colnames(count.mtx),1,5)
info$mutation <- mutation
info$cell_type <- factor(info$cell_type, levels = c(group1, group2))

info


dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dds <- DESeq(dds)
res <- results(dds)
res %>% dim()
res <- data.frame(res)
res %>% filter(padj <= 0.05) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 2 | log2FoldChange <= -2) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5 | log2FoldChange <= -1.5) %>% dim()
res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1) %>% dim()
res %>% head()
res <- res[,c(2,5,6)]
write.csv(res,paste0('rds/ATAC_seq/atac_seq.readvalue.and.normalized/res/', group1,'_',group2, '.csv'))

########################################################
dir <- 'rds/ATAC_seq/atac_seq.readvalue.and.normalized/res_csv'
list.files(dir, full.names = T)
i <- 1
sample.name <-substr(list.files(dir, full.names = T)[i], nchar(list.files(dir, full.names = T)[i])-14,
                     nchar(list.files(dir, full.names = T)[i])-4)
sample.name 

for(i in 1:5){
  sample.name <-substr(list.files(dir, full.names = T)[i], 
                       nchar(list.files(dir, full.names = T)[i])-14,
                       nchar(list.files(dir, full.names = T)[i])-4)
  csv <- read.csv(list.files(dir, full.names = T)[i], header = T, row.names = 1)
  assign(sample.name, csv)
  print(sample.name)
}

wt_tc_ko_tc %>% head()







###################################################################
res.group <- data.frame(matrix(nrow = nrow(res), ncol = 5), row.names = rownames(res))
colnames(res.group) <- c('nai_te','nai_tc',
                         'wt_te_tc','wt_ko_te','wt_ko_tc')
res.group %>% head()


count.mtx %>% head()
res %>% head()
#res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 1.5) %>% dim()
#res %>% filter(padj <= 0.05) %>% filter(log2FoldChange <= -1.5) %>% dim()
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
