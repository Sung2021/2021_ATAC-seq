

## volcano plot
## padj not NA : 52105 peaks
res %>% filter(!(padj == 'NA')) %>% 
  ggplot(aes(log2FoldChange, -log10(padj))) + geom_point(size=0.5) +
  geom_vline(xintercept = c(-log2(1.2),log2(1.2)), color='orange') +
  geom_vline(xintercept = c(-log2(1.5),log2(1.5)), color='red') +
  geom_vline(xintercept = c(-log2(2),log2(2)), color='pink') +
  geom_hline(yintercept = c(-log10(0.05)), color='red') +
  theme_classic()


res %>% filter(!(padj == 'NA')) %>% 
  ggplot(aes(log2FoldChange, -log10(padj))) + geom_point(size=0.5) +
  geom_vline(xintercept = c(-log2(1.2),log2(1.2)), color='orange') +
  geom_vline(xintercept = c(-log2(1.5),log2(1.5)), color='red') +
  geom_vline(xintercept = c(-log2(2),log2(2)), color='pink') +
  geom_hline(yintercept = c(-log10(0.05)), color='red') +
  xlim(c(-3,3)) + ylim(c(-0.5,50)) +
  theme_classic()

res %>% filter(!(pvalue == 'NA')) %>% 
  ggplot(aes(log2FoldChange, -log10(pvalue))) + geom_point(size=0.5) +
  geom_vline(xintercept = c(-log2(1.5),log2(1.5)), color='red') +
  geom_vline(xintercept = c(-log2(2),log2(2)), color='pink') +
  geom_hline(yintercept = c(-log10(0.05)), color='red') +
  theme_classic()


c1 <- res %>% filter(pvalue < 0.05 & log2FoldChange > log2(1.5)) %>% rownames()
c2 <- res %>% filter(pvalue < 0.05 & log2FoldChange < -log2(1.5)) %>% rownames()
c3 <- res %>% filter(pvalue < 0.05 & log2FoldChange > log2(2)) %>% rownames()
c4 <- res %>% filter(pvalue < 0.05 & log2FoldChange < -log2(2)) %>% rownames()
for(i in 1:4){
  list(c1,c2,c3,c4)[[i]] %>% write.csv(paste0('RNA_seq/Icos/gene_list.c',i,'.csv'))
}


# add gene name information
res$gene <- '' ## genes unwanted to plot
genes <- c('Ctla4','Icos','Sell','Il2ra','Runx3', 'Tnfrsf4','Nkg7',
           'Il21','Cxcr5','Nr4a1')
res %>% filter(padj < 0.05 & log2FoldChange > log2(2)) %>% rownames()
genes <- c('Ifi204','Il2ra','Ifitm3',
           'Ccr4','Igfbp4','Cd7',
           'Ifi211','Sox12',
           'Plac8','Cish','Pdlim4',
           'Itgb3')
genes <- c('Nrn1')
genes <- c('Il17ra','Ulk1','Il4ra','Lag3',
           'Hdac6','Ccr7','Id2','Bcl2l1')
genes <- c('Ctla4','Icos','Sell','Runx3',
           'Selplg','Tnfrsf4','Ptgir','Irf4')
res$gene <- ''
res[rownames(res) %in% genes,]$gene <- genes ## genes wanted to plot
res$color <- 'grey'
res[(res %>% filter(pvalue < 0.05 & log2FoldChange > log2(1.5)) %>% 
       rownames()),]$color <- 'sig'
res[(res %>% filter(pvalue < 0.05 & log2FoldChange < -log2(1.5)) %>% 
       rownames()),]$color <- 'sig'
## draw the plot
res %>% ggplot(aes(log2FoldChange, -log10(pvalue), 
                   label=gene)) + geom_point(size=0.5, alpha=0.5, aes(color=color)) +
  scale_color_manual(values = c('grey','black')) +
  geom_vline(xintercept = c(-log2(1.5),log2(1.5)), color='red') +
  geom_hline(yintercept = c(-log10(0.05)), color='red')  +
  geom_text(hjust = -0.05, nudge_x = 0.05) +
  xlim(c(-3,3)) + ylim(c(-0.5,50)) +
  theme_classic() + theme(legend.position = 'none')

test <- read.csv('ATAC_seq/Icos.2022.06/modified.bed/test.csv')
test[1:3,]
table(test$X.1, test$DEG)

test %>% filter(X.1 =='B') %>% filter(DEG == 'up.fc.2') %>% select(SYMBOL) %>% unique()
test %>% filter(X.1 =='B') %>% filter(DEG == 'dn.fc.2') %>% 
  select(SYMBOL,genomic.info) %>% table()
test %>% write.csv('ATAC_seq/Icos.2022.06/modified.bed/A.B.genomic.DEG.22.07.12.csv')

genes ='Plac8'
g3
g4
test[test$SYMBOL %in% genes,]

test %>% filter(X.1 =='A') %>% filter(DEG == 'dn.fc.1.5')
test %>% filter(X.1 =='B') %>% filter(DEG == 'dn.fc.1.5')

rbind(test %>% filter(X.1 =='A') %>% filter(DEG == 'dn.fc.1.5'),
      test %>% filter(X.1 =='B') %>% filter(DEG == 'dn.fc.1.5')) %>% 
  write.csv('ATAC_seq/Icos.2022.06/modified.bed/test1.csv')

rbind(test %>% filter(X.1 =='A') %>% filter(DEG == 'dn.fc.2'),
      test %>% filter(X.1 =='B') %>% filter(DEG == 'dn.fc.2')) %>% 
  write.csv('ATAC_seq/Icos.2022.06/modified.bed/test2.csv')

rbind(test %>% filter(X.1 =='A') %>% filter(DEG == 'up.fc.2'),
      test %>% filter(X.1 =='B') %>% filter(DEG == 'up.fc.2')) %>% 
  write.csv('ATAC_seq/Icos.2022.06/modified.bed/test4.csv')

rbind(test %>% filter(X.1 =='A') %>% filter(DEG == 'up.fc.1.5'),
      test %>% filter(X.1 =='B') %>% filter(DEG == 'up.fc.1.5')) %>% 
  write.csv('ATAC_seq/Icos.2022.06/modified.bed/test3.csv')





bed <- read.table('ATAC_seq/Icos.2022.06/modified.bed/Icos.ATAC_seq.B.bed', 
                  row.names = 4)
colnames(bed)[1:3] <- c('seqnames','start','end')
rownames(bed) <- paste0(bed$seqnames,'_',
                        bed$start,'_',
                        bed$end)
bed$location.name <- rownames(bed)
bed[1:3,]
bed %>% dim()


location.info <- read.csv('ATAC_seq/location.info.22.07.11.csv',
                          row.names = 1)
location.info[1:3,]

## findoverlaps
peaks.gr <- bed[,c(1:3)] %>% makeGRangesFromDataFrame(keep.extra.columns = T)
location.gr <- location.info %>% makeGRangesFromDataFrame(keep.extra.columns = T)

minimum = 50 
ol <- findOverlaps(location.gr, peaks.gr, minoverlap = minimum) %>% data.frame()
ol[1:3,]
location.info[ol$queryHits,][1:3,]
bed[ol$subjectHits,][1:3,]
tmp <- cbind(location.info[ol$queryHits,],
             bed[ol$subjectHits,])
tmp <- tmp[order(tmp$location.name),]
tmp %>% dim()
tmp$genomic.info %>% table()
tmp[1:3,]
tmp <- tmp[order(tmp$location.name),]
tmp %>% write.csv('ATAC_seq/Icos.2022.06/modified.bed/Icos.ATAC_seq.B.bed.annotation.csv')

a =rle(sort(tmp$SYMBOL))
count.info <- data.frame(a$values,
                         a$lengths)
count.info %>% filter(a.lengths > 3)
count.info %>% ggplot(aes(a.lengths)) + geom_histogram(binwidth = 1)
count.info %>% top_n(20, a.lengths)
count.info %>% dim()
count.info %>% write.csv('ATAC_seq/Icos.2022.06/modified.bed/B.bed.count.info.csv')


rna.seq <- read.csv('RNA_seq/Icos/Icos.72hr.RNA_seq.6378.genes.res.csv', row.names = 1)
rna.seq[1:3,]

tmp <- read.csv('ATAC_seq/Icos.2022.06/modified.bed/Icos.ATAC_seq.B.bed.annotation.csv',
                row.names = 1)
tmp[1:3,]
tmp$RNA_seq_log2FC <- ''
tmp[tmp$SYMBOL %in% rownames(rna.seq),]
i <- 1
for(i in 1:nrow(tmp)){if(tmp$SYMBOL[i] %in% rownames(rna.seq)){
  print(tmp$SYMBOL[i])
  gene <- tmp$SYMBOL[i]
  fc <- rna.seq[gene,]$log2FoldChange
  tmp[i,]$RNA_seq_log2FC <-fc
}
  else{
    print('no')
  }
}
tmp[1:10,]
tmp %>% write.csv('ATAC_seq/Icos.2022.06/modified.bed/test.csv')      



###########################################################

### This is a chunk of codes for 
### ATAC-seq analysis
#####################################
setwd('~/Desktop/HMH/rds')

library(ggplot2)
library(dplyr)
library(DESeq2)
library(GenomicRanges)
library(rtracklayer)

bed <- read.table('ATAC_seq/Icos.2022.06/modified.bed/Icos.ATAC_seq.A.bed', 
                  row.names = 4)
bed[1:3,]
bed %>% dim()

gene.info <- read.csv('RNA_seq/gene_genomic_info.csv',row.names = 1)
tss.info <- read.csv('ATAC_seq/gene.tss.info.22.06.30.csv',row.names = 1)
dis.info <- read.csv('ATAC_seq/gene.distal.info.22.06.30.csv', row.names = 1)

gene.info[1:3,]
gene.info$location_name <- paste0(gene.info$seqnames,'_',
                                  gene.info$start, '_',
                                  gene.info$end,'_',
                                  gene.info$SYMBOL,'_',
                                  'gene_body')
gene.info <- gene.info[,c(2:5,7)]

tss.info <- tss.info[,c(1:4)]
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
                            'dis.dn')
tmp[1:3,]
tmp.dn <- tmp


gene.info$genomic.info <- 'gene_body'
tss.info$genomic.info <-'promoter'
tmp.up$genomic.info <- 'dis.up'
tmp.dn$genomic.info <- 'dis.dn'
gene.info[1:3,]
tss.info[1:3,]
tmp.up[1:3,]
tmp.dn[1:3,]



location.info <- rbind(gene.info,tss.info,tmp.up, tmp.dn)
location.info[1:10,]
location.info$genomic.info %>% table()
rownames(location.info) <- location.info$location_name
location.info[1:3,]
location.info %>% dim()
location.info %>% write.csv('ATAC_seq/location.info.22.07.11.csv')
location.info <- read.csv('ATAC_seq/location.info.22.07.11.csv', row.names = 1)
location.info[1:3,]

location.info <- location.info[,c(2:5)]
location.info[1:3,]
location.info %>% write.csv('ATAC_seq/location.info.22.07.11.csv')
location.info <- read.csv('ATAC_seq/location.info.22.07.11.csv',
                          row.names = 1)
location.info[1:3,]

anno.info <- 'gene_body'
anno.info <- 'dis.up'
anno.info <- 'dis.dn'
anno.info <- 'promoter'

rows <- grep(anno.info, location.info$location_name)
location.info[rows,]

anno.list <- read.csv('ATAC_seq/Icos.2022.06/modified.bed/Icos.ATAC_seq.A.bed.annotation.csv',
                      row.names = 1)
anno.list[1:3,]
anno.list$genomic.info <- ''
anno.info <- c('gene_body','dis.up','dis.dn','promoter')
for(i in 1:4){
  rows <- grep(anno.info[i], rownames(anno.list))
  anno.list[rows,]$genomic.info <- as.character(anno.info[i])
}
anno.list$genomic.info %>% table()



##################################################

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

#################################################################
setwd('~/Desktop/HMH/rds/')
library(dplyr)
library(ggplot2)
library(DESeq2)
library(reshape)

read <- read.csv('RNA_seq/Icos/Icos.72hr.RNA_seq.featurecounts.40bp.trimming.raw.22.06.29.csv',
                 row.names = 1)
read[1:3,]

fpkm <- read.csv('RNA_seq/Icos/Icos.72hr.RNA_seq.featurecounts.40bp.trimming.fpkm.22.06.29.csv',
                 row.names = 1)
fpkm[1:3,]
fpkm %>% dim()

fpkm <- read.csv('RNA_seq/Icos/Icos.72hr.RNA_seq.fpkm.mean.over1.6378.genes.22.06.30.csv',
                 row.names = 1)
fpkm[1:3,]
fpkm %>% dim()


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
tmp[1:3,]
norm.mtx <- fpkm/tmp
norm.mtx[1:3,]
           

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
