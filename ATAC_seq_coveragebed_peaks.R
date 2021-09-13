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
combined[,'c1'] <- condition1[,5]
combined[,'c2'] <- condition2[,5]
combined[,'c3'] <- condition3[,5]
combined[,'c4'] <- condition4[,5]
combined[,'c5'] <- condition5[,5]
combined[,'c6'] <- condition6[,5]
combined[,'c7'] <- condition7[,5]
combined[,'c8'] <- condition8[,5]
combined[,'c9'] <- condition9[,5]
combined[,'c10'] <- condition10[,5]
combined[,'c11'] <- condition11[,5]
combined[,'c12'] <- condition12[,5]
combined[,'c13'] <- condition13[,5]

combined[1:3,]
cov.beds
colnames(combined)[4:16] <- c('ko_tcm1','ko_tcm2','ko_tcm3',
                              'ko_tem1','ko_tem2',
                              'naive1','naive2',
                              'wt_tcm1','wt_tcm2','wt_tcm3',
                              'wt_tem1','wt_tem2','wt_tem3')
## condition1,2,3 coverage
#combined[,'c1_cov'] <- condition1[,8]
#combined[,'c2_cov'] <- condition2[,8]
#combined[,'c3_cov'] <- condition3[,8]

######## choose only the peaks on chr1-chr20 (mouse)
chrs <- paste0('chr',1:20)
combined <- combined[combined$chr %in% chrs,]
######## inspect the combined peaks for unwanted chromosomes
combined[grep('_', combined$chr),]
combined[grep('M', combined$chr),]
combined[grep('X', combined$chr),]
combined[grep('Y', combined$chr),]

## intermediate matrix for the normalized score
## median of score of each condition
tmp.md <- data.frame(matrix(nrow = length(rownames(combined)), 
                            ncol = length(colnames(combined))-3))
colnames(tmp.md) <- colnames(combined)[4:16]
for(i in 1:13){
  tmp.md[,i] <- median(combined[,3+i])
}
tmp.norm <- data.frame(matrix(nrow = length(rownames(combined)), 
                            ncol = length(colnames(combined))-3))
colnames(tmp.norm) <- colnames(combined)[4:16]
tmp.norm <- combined[,4:16]
tmp.norm <- tmp.norm/tmp.md
tmp.norm <- tmp.norm[,c(6:7,11:13,8:10,4:5,1:3)]
combined.norm <- cbind(combined[,1:3], tmp.norm)

saveRDS(combined.norm, '~/Desktop/HMH/rds/2021.09.13.atac_seq.normalized.peaks.rds')
write.csv(combined.norm, '~/Desktop/HMH/rds/2021.09.13.atac_seq.normalized.peaks.csv')

### inspect the distribution of the normalized peaks
rowMeans(combined.norm[,4:16])
hist(rowMeans(combined.norm[,4:16]), breaks = 100)
summary(rowMeans(combined.norm[,4:16]))
boxplot(rowMeans(combined.norm[,4:16]))
combined.norm[rowMeans(combined.norm[,4:16]) > 200,]

##### remove one expression outlier chr17 34000246 34002280 : which 16825
##### remove one expression outlier 25471 chr17 34000185 34002280
combined.norm <- combined.norm[-25471,] ### 58693 left
dim(combined.norm)

## need to renormalize again
combined[25471,]
combined <- combined[-25471,]

## intermediate matrix for the normalized score
## median of score of each condition
tmp.md <- data.frame(matrix(nrow = length(rownames(combined)), 
                            ncol = length(colnames(combined))-3))
colnames(tmp.md) <- colnames(combined)[4:16]
for(i in 1:13){
  tmp.md[,i] <- median(combined[,3+i])
}
tmp.norm <- data.frame(matrix(nrow = length(rownames(combined)), 
                              ncol = length(colnames(combined))-3))
colnames(tmp.norm) <- colnames(combined)[4:16]
tmp.norm <- combined[,4:16]
tmp.norm <- tmp.norm/tmp.md
tmp.norm <- tmp.norm[,c(6:7,11:13,8:10,4:5,1:3)]
combined.norm <- cbind(combined[,1:3], tmp.norm)

colnames(combined.norm)
### inspect the distribution of the normalized peaks
ggplot(combined.norm, aes(naive1,naive2)) + geom_point() + ggpubr::stat_cor(method = 'pearson')
ggplot(combined.norm, aes(wt_tem1,wt_tcm1)) + geom_point() + ggpubr::stat_cor(method = 'pearson')

saveRDS(combined.norm, '~/Desktop/HMH/rds/2021.09.13.atac_seq.normalized.peaks.rds')
write.csv(combined.norm, '~/Desktop/HMH/rds/2021.09.13.atac_seq.normalized.peaks.csv')
