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
###### information sheet
atac_norm <- readRDS('~/Desktop/HMH/rds/2021.09.13.atac_seq.normalized.peaks.rds')
atac_raw <- readRDS('~/Desktop/HMH/rds/2021.09.13.atac_seq.raw.peaks.rds')
atac_norm[1:3,]
dim(atac_norm)
colnames(atac_norm)
## change the row names to peak row names
rownames(atac_norm) <- paste0('peak_',rownames(atac_norm))
rownames(atac_raw) <- paste0('peak_',rownames(atac_raw))

## information sheet
## wt_tem vs wt_tcm
count.mtx <- atac_raw[, grep('wt_',colnames(atac_raw))]
info <- data.frame(matrix(nrow = 6, ncol = 3))
colnames(info) <- c('sample', 'type','mutation')
info$sample <- colnames(count.mtx)
info$type <- c(rep('tem',3), rep('tcm',3))
info$type <- factor(info$type, levels = c('tcm','tem'))
info$mutation <- 'wt'
info
str(info)
## create dds object
dds <- DESeqDataSetFromMatrix(count.mtx, info, ~type)
#dim(dds)
dds <- DESeq(dds)
res <- results(dds)
#dim(res)
res <- data.frame(res)
colnames(res) <- paste0('wt_tcm_vs_wt_tem_', colnames(res))
#write.csv(res, '~/Desktop/HMH/rds/2021.09.13.wt_tcm_vs_wt_tem_deg.csv')
## save data for later use
res.wt_tcm_wt_tem <- res

### check the details of the DEG output

res[c('peak_8460','peak_1021'),]
count.mtx['peak_8460',grep('wt_',colnames(count.mtx))] ## Tcf7 : tcm more
count.mtx['peak_1021',grep('wt_',colnames(count.mtx))] ## Ctla4 : tcm more



## information sheet
## wt_tem vs ko_tem
count.mtx <- atac_raw[, grep('_tem',colnames(atac_raw))]
count.mtx[1:3,]
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'type','mutation')
info$sample <- colnames(count.mtx)
info$type <-'tem'
info$mutation <- c(rep('ko',2), rep('wt',3))
info$mutation <- factor(info$mutation,levels = c('wt','ko'))
info
str(info)
## create dds object
dds <- DESeqDataSetFromMatrix(count.mtx, info, ~mutation)
#dim(dds)
dds <- DESeq(dds)
res <- results(dds)
#dim(res)
res <- data.frame(res)
colnames(res) <- paste0('wt_tem_vs_ko_tem_', colnames(res))

# write.csv(res, '~/Desktop/HMH/rds/2021.09.13.wt_tem_ko_tem_deg.csv')


res['peak_8458',]
count.mtx['peak_8458',grep('_tem',colnames(count.mtx))]

## save data for later use
res.wt_tem_ko_tem <- res




## information sheet
## wt_tcm vs ko_tcm
count.mtx <- atac_raw[, grep('_tcm',colnames(atac_raw))]
count.mtx[1:3,]
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'type','mutation')
info$sample <- colnames(count.mtx)
info$type <-'tcm'
info$mutation <- c(rep('ko',3), rep('wt',3))
info$mutation <- factor(info$mutation,levels = c('wt','ko'))
info
str(info)
## create dds object
dds <- DESeqDataSetFromMatrix(count.mtx, info, ~mutation)
#dim(dds)
dds <- DESeq(dds)
res <- results(dds)
#dim(res)
res <- data.frame(res)
colnames(res) <- paste0('wt_tcm_vs_ko_tcm_', colnames(res))

# write.csv(res, '~/Desktop/HMH/rds/2021.09.13.wt_tcm_ko_tcm_deg.csv')


res['peak_8458',]
count.mtx['peak_8458',grep('_tcm',colnames(count.mtx))]

## save data for later use
res.wt_tcm_ko_tcm <- res




## information sheet
## naive vs wt_tem
colnames(atac_raw)
count.mtx <- atac_raw[,c('naive1', 'naive2',"wt_tem1" ,"wt_tem2", "wt_tem3")]
count.mtx[1:3,]
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'type','mutation')
info$sample <- colnames(count.mtx)
info$type <- c(rep('naive',2),
               rep('tem',3))
info$mutation <- 'wt'
info$type <- factor(info$type, levels = c('naive','tem'))
info
str(info)
## create dds object
dds <- DESeqDataSetFromMatrix(count.mtx, info, ~type)
#dim(dds)
dds <- DESeq(dds)
res <- results(dds)
#dim(res)
res <- data.frame(res)
colnames(res) <- paste0('naive_vs_wt_tem_', colnames(res))

# write.csv(res, '~/Desktop/HMH/rds/2021.09.13.naive_wt_tem_deg.csv')


res['peak_8458',]
count.mtx['peak_8458',grep('_tcm',colnames(count.mtx))]

## save data for later use
res.naive_tem <- res


## information sheet
## naive vs wt_tcm
colnames(atac_raw)
count.mtx <- atac_raw[,c('naive1', 'naive2',"wt_tcm1" ,"wt_tcm2", "wt_tcm3")]
count.mtx[1:3,]
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'type','mutation')
info$sample <- colnames(count.mtx)
info$type <- c(rep('naive',2),
               rep('tcm',3))
info$mutation <- 'wt'
info$type <- factor(info$type, levels = c('naive','tcm'))
info
str(info)
## create dds object
dds <- DESeqDataSetFromMatrix(count.mtx, info, ~type)
#dim(dds)
dds <- DESeq(dds)
res <- results(dds)
#dim(res)
res <- data.frame(res)
colnames(res) <- paste0('naive_vs_wt_tcm_', colnames(res))

# write.csv(res, '~/Desktop/HMH/rds/2021.09.13.naive_wt_tcm_deg.csv')

res['peak_8458',]
count.mtx['peak_8458',]

## save data for later use
res.naive_tcm <- res

###########################################################################
###########################################################################
###########################################################################

res.wt_tcm_wt_tem[1:3,]
res.wt_tem_ko_tem[1:3,]
res.wt_tcm_ko_tcm[1:3,]
res.naive_tem[1:3,]
res.naive_tcm[1:3,]

res.combined <- cbind(res.wt_tcm_wt_tem,
                      res.wt_tem_ko_tem,
                      res.wt_tcm_ko_tcm,
                      res.naive_tem,
                      res.naive_tcm)
##### save combined data as csv
write.csv(res.combined, '~/Desktop/HMH/rds/2021.09.13.atac_seq_deg.combined.csv')
##### save log2FC only
write.csv(res.combined[,grep('log2', colnames(res.combined))], '~/Desktop/HMH/rds/2021.09.13.atac_seq_deg.combined.log2.csv')

# res.combined[,grep('log2', colnames(res.combined))]['peak_8458',]
# res.combined[1:3,]
