library(ggplot2)
library(dplyr)
library(DESeq2)
library(reshape)
setwd('~/Desktop/HMH/rds')
theme_set(theme_bw())

##############################################################
## zscore function 
zscore <- function(input.data = input.data){
  input.data.rowsums <- rowSums(input.data)
  input.data.mean <- rowMeans(input.data)
  input.data.sd <- matrixStats::rowSds(as.matrix(input.data))
  names(input.data.sd) <- rownames(input.data)
  zscore <- (input.data-input.data.mean)/input.data.sd
  return(zscore)
}


##############################################################

norm.peak <- read.csv('rds/ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/atac_seq.21.11.29.read.norm.by.colmeans.with.naive.csv', 
                      row.names = 1)
norm.peak %>% head()
norm.peak %>% dim()


####################################################################
read.value <- read.csv('ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/atac_seq.21.11.29.read.raw.with.naive.csv',
                       row.names = 5)
read.value %>% head()
read.value %>% dim()

####################################################################
norm.peak.log2fc <- read.csv('rds/ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/atac_seq.21.11.30.norm.log2fc.all.csv',
                             row.names = 1)
norm.peak.log2fc %>% head()
norm.peak.log2fc %>% dim()

####################################################################

sig <- read.csv('rds/ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/atac_seq.21.11.29.sig.peaks.norm.csv',
                row.names = 1)
sig %>% head()
sig %>% dim()

sig.norm <- norm.peak.log2fc[rownames(sig),]
sig.norm %>% head()
sig.norm %>% dim()



count.mtx <- read.value[,grep('wt_te|wt_tc', colnames(read.value))]
group1 <- 'wt_te'
group2 <- 'wt_tc'
count.mtx %>% head()

### information sheet
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- colnames(count.mtx)
info$cell_type <- substr(colnames(count.mtx),1,5)
info$mutation <- substr(colnames(count.mtx),1,2)
info$cell_type <- factor(info$cell_type, levels = c(group1, group2))
info

dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dds <- DESeq(dds)
res <- results(dds)
res %>% dim()
res <- data.frame(res)
res <- res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 2 | log2FoldChange <= -2)
res <- res[,c(2,5:6)]
#colnames(res) <- paste0(paste0(group1,'_' ,group2,'_'), colnames(res))
res$comparison <- paste0(group1,'_' ,group2)
res %>% head()
res %>% dim()

res.3 <- res

#############################################################################
count.mtx <- read.value[,grep('wt_te|ko_te', colnames(read.value))]
group1 <- 'wt_te'
group2 <- 'ko_te'
count.mtx %>% head()

### information sheet
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- colnames(count.mtx)
info$cell_type <- substr(colnames(count.mtx),1,5)
info$mutation <- substr(colnames(count.mtx),1,2)
info$cell_type <- factor(info$cell_type, levels = c(group1, group2))
info

dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dds <- DESeq(dds)
res <- results(dds)
res %>% dim()
res <- data.frame(res)
res <- res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 2 | log2FoldChange <= -2)
res <- res[,c(2,5:6)]
#colnames(res) <- paste0(paste0(group1,'_' ,group2,'_'), colnames(res))
res$comparison <- paste0(group1,'_' ,group2)
res %>% head()
res.4 <- res


res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= log2(2) | log2FoldChange <= -log2(2)) %>% dim()

#####################################################################
count.mtx <- read.value[,grep('wt_tc|ko_tc', colnames(read.value))]
group1 <- 'wt_tc'
group2 <- 'ko_tc'
count.mtx %>% head()

### information sheet
info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- colnames(count.mtx)
info$cell_type <- substr(colnames(count.mtx),1,5)
info$mutation <- substr(colnames(count.mtx),1,2)
info$cell_type <- factor(info$cell_type, levels = c(group1, group2))
info

dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dds <- DESeq(dds)
res <- results(dds)
res %>% dim()
res <- data.frame(res)
res <- res %>% filter(padj <= 0.05) %>% filter(log2FoldChange >= 2 | log2FoldChange <= -2)
res <- res[,c(2,5:6)]
#colnames(res) <- paste0(paste0(group1,'_' ,group2,'_'), colnames(res))
res$comparison <- paste0(group1,'_' ,group2)
res %>% head()
res.5 <- res

res.3 %>% dim()
res.4 %>% dim()
res.5 %>% dim()


loci.3 <- res.3 %>% rownames()
loci.4 <- res.4 %>% rownames()
loci.5 <- res.5 %>% rownames()
loci.all <- list(loci.3=loci.3,
                 loci.4=loci.4,
                 loci.5=loci.5)
union(loci.3,union(loci.4,loci.5)) %>% length()
sig.peaks <- union(loci.3,union(loci.4,loci.5))

library(ggvenn)
ggvenn(loci.all)
ggvenn(loci.all,stroke_size = 0.5, set_name_size = 4, show_percentage = F)


#############################################################################

sig.peaks.norm <- norm.peak[sig.peaks,]
sig.peaks.norm %>% dim()
sig.peaks.norm %>% head()
## zscore function 
zscore <- function(input.data = input.data){
  input.data.rowsums <- rowSums(input.data)
  input.data.mean <- rowMeans(input.data)
  input.data.sd <- matrixStats::rowSds(as.matrix(input.data))
  names(input.data.sd) <- rownames(input.data)
  zscore <- (input.data-input.data.mean)/input.data.sd
  return(zscore)
}
input.data <- sig.peaks.norm[,-c(1:2)]
input.data %>% dim()

fit <- kmeans(zscore(input.data), centers = 8, nstart = 10)
fit.cluster <- fit$cluster %>% data.frame()
colnames(fit.cluster) <- 'cluster'
fit.cluster <- fit.cluster %>% arrange(cluster)
fit.cluster$cluster <- factor(fit.cluster$cluster)
fit.cluster %>% head()
df.anno <- fit.cluster
fit.cluster$cluster %>% table()


pheatmap::pheatmap(zscore(log2(input.data[rownames(fit.cluster),]+1)), cluster_cols = F, 
                   cluster_rows = F, show_rownames = F, annotation_row = df.anno,
                   col = colorRampPalette(c("navy", "white", "red"))(1000),
                   fontsize_row = 6,fontsize_col = 10)




sig.peak.norm.cluster <- sig.peaks.norm[rownames(fit.cluster),] %>% cbind(df.anno)
write.csv(sig.peak.norm.cluster, 'rds/ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/atac_seq.21.11.30.sig.peak.cluster.csv')
input.data[rownames(fit.cluster),] %>% head()


################################################################################


sig.peak.norm.cluster <- read.csv('rds/ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/atac_seq.21.11.30.sig.peak.cluster.csv',
                                  row.names = 1)
sig.peak.norm.cluster %>% head()
sig.peak.norm.cluster %>% select(cluster) %>% table()


gene.exp.func <- function(gene,input.data=fpkm){
  gene <- as.character(gene)
  if(gene %in% rownames(input.data)){
    df <- input.data[gene,] %>% melt() %>% data.frame()
    df$condition <- substr(df$variable,1,5)
    df2 <- df %>% group_by(condition) %>% summarize(avg=mean(value), sd=sd(value))
    p <- ggplot(df2, aes(condition,avg)) + geom_col()+ geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd)) + 
      ylab(gene) + theme_classic()
    print(p)
  }
  else{
    message('no gene in this dataset')
  }
}

locus <- 'union_peak_1140'
locus.enr.func <- function(locus, input.data=sig.peak.norm.cluster[,-c(13:14)]){
  locus <- as.character(locus)
  input.data <- sig.peak.norm.cluster[,-c(13:14)]
  df <- input.data[locus,] %>% melt() %>% data.frame()
  df$condition <- substr(df$variable,1,5)
  df$condition <- factor(df$condition,levels = c('naive','wt_te','wt_tc','ko_te','ko_tc'))
  df2 <- df %>% group_by(condition) %>% summarize(avg=mean(value), sd=sd(value))
  p <- ggplot(df2, aes(condition,avg)) + geom_col()+ geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd)) + 
    ylab(locus) + theme_classic()+ ggtitle(levels(sig.peak.norm.cluster$group)[i])
  print(p)
}
locus <- as.character(locus)
input.data <- sig.peak.norm.cluster[,-c(13:14)]
df <- input.data[locus,] %>% melt() %>% data.frame()
df$condition <- substr(df$variable,1,5)
df$condition <- factor(df$condition,levels = c('naive','wt_te','wt_tc','ko_te','ko_tc'))
df2 <- df %>% group_by(condition) %>% summarize(avg=mean(value), sd=sd(value))
p <- ggplot(df2, aes(condition,avg)) + geom_col()+ geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd)) + 
  ylab(locus) + theme_classic()
print(p)



sig.peak.norm.cluster <- read.csv('rds/ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/atac_seq.21.11.30.sig.peak.cluster.csv',
                                  row.names = 1)
sig.peak.norm.cluster %>% head()
sig.peak.norm.cluster %>% dim()
sig.peak.norm.cluster %>% select(cluster) %>% table()
sig.peak.norm.cluster[sig.peak.norm.cluster$cluster == '1','group'] <- 'ko_tem_enr'
sig.peak.norm.cluster[sig.peak.norm.cluster$cluster == '2','group'] <- 'ko_enr_2'
sig.peak.norm.cluster[sig.peak.norm.cluster$cluster == '3','group'] <- 'wt_tcm_enr'
sig.peak.norm.cluster[sig.peak.norm.cluster$cluster == '4','group'] <- 'ko_enr_1'
sig.peak.norm.cluster[sig.peak.norm.cluster$cluster == '5','group'] <- 'ko_enr_3'
sig.peak.norm.cluster[sig.peak.norm.cluster$cluster == '6','group'] <- 'wt_enr_2'
sig.peak.norm.cluster[sig.peak.norm.cluster$cluster == '7','group'] <- 'wt_tem_enr'
sig.peak.norm.cluster[sig.peak.norm.cluster$cluster == '8','group'] <- 'wt_enr_1'
sig.peak.norm.cluster$group <- factor(sig.peak.norm.cluster$group, 
                                      levels = c('wt_enr_1',
                                                 'wt_enr_2',
                                                 'ko_enr_1',
                                                 'ko_enr_2',
                                                 'ko_enr_3',
                                                 'wt_tem_enr',
                                                 'wt_tcm_enr',
                                                 'ko_tem_enr'))
sig.peak.norm.cluster %>% select(group) %>% table()



input.data <- sig.peak.norm.cluster[,-c(1:2,13:14)]
input.data %>% dim()
input.data %>% head()


df.anno <- sig.peak.norm.cluster[,c(14,13)]
df.anno %>% head()

pheatmap::pheatmap(zscore(log2(input.data+1)), cluster_cols = F, 
                   cluster_rows = F, show_rownames = F, annotation_row = df.anno,
                   col = colorRampPalette(c("navy", "white", "red"))(1000),
                   fontsize_row = 6,fontsize_col = 10)

###########################################################################
########### 2021.12.09 ############

sig.peak.norm <- read.csv('ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/atac_seq.21.11.29.sig.peaks.norm.csv',
                          row.names = 1)
sig.peak.norm %>% head()

input.data <- sig.peak.norm
numberofcluster <- 7
fit <- kmeans(zscore(input.data), centers = numberofcluster, nstart = 10)
fit.cluster <- fit$cluster %>% data.frame()
colnames(fit.cluster) <- 'cluster'
fit.cluster <- fit.cluster %>% arrange(cluster)
fit.cluster$cluster <- factor(fit.cluster$cluster)
fit.cluster %>% head()
df.anno <- fit.cluster
fit.cluster$cluster %>% table()


pheatmap::pheatmap(zscore(log2(input.data[rownames(fit.cluster),]+1)), cluster_cols = F, 
                   cluster_rows = F, show_rownames = F, annotation_row = df.anno,
                   col = colorRampPalette(c("navy", "white", "red"))(1000),
                   fontsize_row = 6,fontsize_col = 10)
write.csv(fit.cluster, 'ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/sig.peak.norm.cluster7.2021.12.09.csv')



##############################################################################
##### 2021.12.09 atac_seq log2FC log2(2) peaks k-means clustering analysis

norm.peak.log2fc <- read.csv('ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/atac_seq.21.11.30.norm.log2fc.all.csv', 
                             row.names = 1)
norm.peak.log2fc %>% head()
norm.peak.log2fc %>% dim()

peak.tmp <- norm.peak.log2fc[,grep('wt_te_wt_tc', colnames(norm.peak.log2fc))]
peak.tmp %>% dim()
peak.tmp %>% head()
peak.tmp %>% filter(wt_te_wt_tc_padj < 0.05) %>% 
  filter(wt_te_wt_tc_log2FoldChange > log2(2)|wt_te_wt_tc_log2FoldChange < -log2(2)) %>% dim()
loci.3 <- peak.tmp %>% filter(wt_te_wt_tc_padj < 0.05) %>% 
  filter(wt_te_wt_tc_log2FoldChange > log2(2)|wt_te_wt_tc_log2FoldChange < -log2(2)) %>% rownames()

peak.tmp <- norm.peak.log2fc[,grep('wt_te_ko_te', colnames(norm.peak.log2fc))]
peak.tmp %>% dim()
peak.tmp %>% head()
peak.tmp %>% filter(wt_te_ko_te_padj < 0.05) %>% 
  filter(wt_te_ko_te_log2FoldChange > log2(2)|wt_te_ko_te_log2FoldChange < -log2(2)) %>% dim()
loci.4 <- peak.tmp %>% filter(wt_te_ko_te_padj < 0.05) %>% 
  filter(wt_te_ko_te_log2FoldChange > log2(2)|wt_te_ko_te_log2FoldChange < -log2(2)) %>% rownames()

peak.tmp <- norm.peak.log2fc[,grep('wt_tc_ko_tc', colnames(norm.peak.log2fc))]
peak.tmp %>% dim()
peak.tmp %>% head()
peak.tmp %>% filter(wt_tc_ko_tc_padj < 0.05) %>% 
  filter(wt_tc_ko_tc_log2FoldChange > log2(2)|wt_tc_ko_tc_log2FoldChange < -log2(2)) %>% dim()
loci.5 <- peak.tmp %>% filter(wt_tc_ko_tc_padj < 0.05) %>% 
  filter(wt_tc_ko_tc_log2FoldChange > log2(2)|wt_tc_ko_tc_log2FoldChange < -log2(2)) %>% rownames()

loci.all <- list(loci.3=loci.3,
                 loci.4=loci.4,
                 loci.5=loci.5)
union(loci.3,union(loci.4,loci.5)) %>% length()
sig.peaks <- union(loci.3,union(loci.4,loci.5))

library(ggvenn)
ggvenn(loci.all)
ggvenn(loci.all,stroke_size = 0.5, set_name_size = 4, show_percentage = F)


norm.peak.log2fc <- read.csv('ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/atac_seq.21.11.30.norm.log2fc.all.csv', 
                             row.names = 1)
norm.peak.log2fc %>% head()
norm.peak <- norm.peak.log2fc[,grep('read', colnames(norm.peak.log2fc))]
norm.peak <- norm.peak[,-c(1:2)]
norm.peak %>% head()

input.data <- norm.peak[sig.peaks,]
set.seed(123)
numberofcluster <- 5
fit <- kmeans(zscore(input.data), centers = numberofcluster, nstart = 10)
fit.cluster <- fit$cluster %>% data.frame()
colnames(fit.cluster) <- 'cluster'
fit.cluster <- fit.cluster %>% arrange(cluster)
fit.cluster$cluster <- factor(fit.cluster$cluster)
fit.cluster %>% head()
df.anno <- fit.cluster
fit.cluster$cluster %>% table()
pheatmap::pheatmap(zscore(log2(input.data[rownames(fit.cluster),]+1)), cluster_cols = F, 
                   cluster_rows = F, show_rownames = F, annotation_row = df.anno,
                   col = colorRampPalette(c("navy", "white", "red"))(1000),
                   fontsize_row = 6,fontsize_col = 10)
