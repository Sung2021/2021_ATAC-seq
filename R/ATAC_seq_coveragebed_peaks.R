library(ggplot2)

setwd('~/Desktop/Sung_work/sam/')
cov.beds <- list.files(path = '~/Desktop/Sung_work/sam/bed', 
                       pattern = 'coverage.bed$', full.names = T)
## total 13 samples
condition1 <- read.table(cov.beds[[1]])
condition2 <- read.table(cov.beds[[2]])
condition3 <- read.table(cov.beds[[3]])
condition4 <- read.table(cov.beds[[4]])
condition5 <- read.table(cov.beds[[5]])
condition6 <- read.table(cov.beds[[6]])
condition7 <- read.table(cov.beds[[7]])
condition8 <- read.table(cov.beds[[8]])
condition9 <- read.table(cov.beds[[9]])
condition10 <- read.table(cov.beds[[10]])
condition11 <- read.table(cov.beds[[11]])
condition12 <- read.table(cov.beds[[12]])
condition13 <- read.table(cov.beds[[13]])


### generate the combined peak matrix
### merged peak information chr, start, end
combined <- condition1[,c(1:3)]
colnames(combined) <- c('chr','start','end')

cov.beds[[1]]
## condition1,2,3 score
combined[,'c1'] <- condition1[,11]
combined[,'c2'] <- condition2[,11]
combined[,'c3'] <- condition3[,11]
combined[,'c4'] <- condition4[,11]
combined[,'c5'] <- condition5[,11]
combined[,'c6'] <- condition6[,11]
combined[,'c7'] <- condition7[,11]
combined[,'c8'] <- condition8[,11]
combined[,'c9'] <- condition9[,11]
combined[,'c10'] <- condition10[,11]
combined[,'c11'] <- condition11[,11]
combined[,'c12'] <- condition12[,11]
combined[,'c13'] <- condition13[,11]

combined[1:3,]
cov.beds
colnames(combined)[4:16] <- c('ko_tcm1','ko_tcm2','ko_tcm3',
                              'ko_tem1','ko_tem2',
                              'naive1','naive2',
                              'wt_tcm1','wt_tcm2','wt_tcm3',
                              'wt_tem1','wt_tem2','wt_tem3')
## reordered colnames

combined <- combined[,c(1:3,9:10,14:16,11:13,7:8,4:6)]

## condition1,2,3 coverage
#combined[,'c1_cov'] <- condition1[,8]
#combined[,'c2_cov'] <- condition2[,8]
#combined[,'c3_cov'] <- condition3[,8]

######## choose only the peaks on chr1-chr20 (mouse)
chrs <- paste0('chr',1:20)
combined <- combined[combined$chr %in% chrs,]

combined[1:3,]
boxplot(rowMeans(combined[,4:16]))
combined[rowMeans(combined[,4:16]) >= 20000,]
## 19136 chr17 34000039 34002278
combined <- combined[-19136,]
dim(combined)
boxplot(rowMeans(combined[,4:16]))
summary(rowMeans(combined[,4:16]))

######## inspect the combined peaks for unwanted chromosomes
combined[grep('_', combined$chr),]
combined[grep('M', combined$chr),]
combined[grep('X', combined$chr),]
combined[grep('Y', combined$chr),]
combined[1:3,]

## intermediate matrix for the normalized score
## median of score of each condition
tmp.md <- data.frame(matrix(nrow = length(rownames(combined)), 
                            ncol = length(colnames(combined))-3))
colnames(tmp.md) <- colnames(combined)[4:16]
for(i in 1:13){
  tmp.md[,i] <- median(combined[,3+i])
}
tmp.md[1:3,]

## matrix to be normalized
tmp.norm <- data.frame(matrix(nrow = length(rownames(combined)), 
                            ncol = length(colnames(combined))-3))
colnames(tmp.norm) <- colnames(combined)[4:16]
tmp.norm <- combined[,4:16]
tmp.norm[1:3,]
tmp.norm <- tmp.norm/tmp.md
tmp.norm[1:3,]

combined.norm <- cbind(combined[,1:3], tmp.norm)
combined.norm[1:3,]
saveRDS(combined.norm, '~/Desktop/HMH/rds/2021.09.20.atac_seq.normalized.peaks.rds')
write.csv(combined.norm, '~/Desktop/HMH/rds/2021.09.20.atac_seq.normalized.peaks.csv')

### combined is the count-matrix for DESeq2
combined[1:3,]
saveRDS(combined, 
        '~/Desktop/HMH/rds/2021.09.20.atac_seq.combined.peaks.rds')
write.csv(combined,
          '~/Desktop/HMH/rds/2021.09.20.atac_seq.combined.peaks.csv')

combined[44165:44167,]



colnames(combined.norm)
### inspect the distribution of the normalized peaks
ggplot(combined.norm, aes(naive1,naive2)) + geom_point() + ggpubr::stat_cor(method = 'pearson')
ggplot(combined.norm, aes(wt_tem1,wt_tcm1)) + geom_point() + ggpubr::stat_cor(method = 'pearson')


###################### analysis of p value cut off #####################
setwd('~/Desktop/Sung_work/sam/macs2/')
cutoff.analysis <- list.files(path = '~/Desktop/Sung_work/sam/macs2', 
                       pattern = '_cutoff_analysis.txt$', full.names = T)
naive <- read.table('naive_cutoff_analysis.txt', header = T)
wt_tem <- read.table('wt_tem_cutoff_analysis.txt', header = T)
wt_tcm <- read.table('wt_tcm_cutoff_analysis.txt', header = T)
ko_tem <- read.table('ko_tem_cutoff_analysis.txt', header = T)
ko_tcm <- read.table('ko_tcm_cutoff_analysis.txt', header = T)

rownames(naive) <- paste0('naive_', rownames(naive))
rownames(wt_tem) <- paste0('wt_tem_', rownames(wt_tem))
rownames(wt_tcm) <- paste0('wt_tcm_', rownames(wt_tcm))
rownames(ko_tem) <- paste0('ko_tem_', rownames(ko_tem))
rownames(ko_tcm) <- paste0('ko_tcm_', rownames(ko_tcm))
tmp <- rbind(naive, wt_tem, wt_tcm,
             ko_tem,ko_tcm)

tmp$sample <- 'NA'
tmp[,'sample'] <- substr(rownames(tmp),1,6)
tmp$sample <- factor(tmp$sample, levels = c('naive_',
                                            'wt_tem','wt_tcm',
                                            'ko_tem','ko_tcm'))
p1 <- ggplot(tmp, aes(factor(pscore), npeaks, color=sample)) + geom_line(aes(group=sample)) +
  theme(axis.text.x= element_text(angle=45)) + ggtitle('raw number of peaks on y axis') + geom_vline(xintercept = 6)
p2 <- ggplot(tmp, aes(factor(pscore), npeaks, color=sample)) + geom_line(aes(group=sample)) +
  theme(axis.text.x= element_text(angle=45)) + scale_y_log10() + ggtitle('log10 scaled number of peaks on y axis') +
  geom_vline(xintercept = 12)
cowplot::plot_grid(p1,p2, ncol = 1)

tmp[tmp$pscore %in% c(3.6),]
tmp[tmp$npeaks <= 40000,]


ggplot(tmp, aes(factor(pscore), npeaks, color=sample)) + geom_line(aes(group=sample)) +
  theme(axis.text.x= element_text(angle=45, size=5)) +facet_wrap(.~sample, ncol = 1) +
  geom_hline(yintercept = c(30000,40000,50000)) + scale_y_log10() + geom_vline(xintercept = c(6,12,20), color='red')

