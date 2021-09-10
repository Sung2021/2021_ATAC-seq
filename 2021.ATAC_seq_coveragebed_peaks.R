library(ggplot2)

setwd('~/Desktop/Sung_work/sam/')
cov.beds <- list.files(path = '~/Desktop/Sung_work/sam', 
                       pattern = 'coverage.bed$', full.names = T)
condition1 <- read.table(cov.beds[[1]])
condition2 <- read.table(cov.beds[[2]])
condition3 <- read.table(cov.beds[[3]])

### generate the combined peak matrix
### merged peak information chr, start, end
combined <- condition1[,c(1:3)]
colnames(combined) <- c('chr','start','end')
## condition1,2,3 score
combined[,'c1'] <- condition1[,5]
combined[,'c2'] <- condition2[,5]
combined[,'c3'] <- condition3[,5]
## condition1,2,3 coverage
combined[,'c1_cov'] <- condition1[,8]
combined[,'c2_cov'] <- condition2[,8]
combined[,'c3_cov'] <- condition3[,8]

######## choose only the peaks on chr1-chr20 (mouse)
chrs <- paste0('chr',1:20)
combined <- combined[combined$chr %in% chrs,]
######## inspect the combined peaks for unwanted chromosomes
combined[grep('_', combined$chr),]
combined[grep('M', combined$chr),]
combined[grep('X', combined$chr),]
combined[grep('Y', combined$chr),]

## median of score of each condition
combined$c1_md <- median(combined$c1)
combined$c2_md <- median(combined$c2)
combined$c3_md <- median(combined$c3)

## intermediate matrix for the normalized score
tmp.df <- combined[,c(4:6, 10:12)]
tmp.df[,c(7:9)] <- 'NA'
colnames(tmp.df)[7:9] <- c(paste0('normalized_','c',1:3))
tmp.df[,7] <- tmp.df$c1/tmp.df$c1_md ## each raw score divided by median of whole raw score
tmp.df[,8] <- tmp.df$c2/tmp.df$c2_md ## each raw score divided by median of whole raw score
tmp.df[,9] <- tmp.df$c3/tmp.df$c3_md ## each raw score divided by median of whole raw score
## add normalized file to combined file
combined <- cbind(combined, tmp.df[,c(7:9)])

##### remove one expression outlier chr17 34000246 34002280 : which 16825
rownames(combined[combined$normalized_c1 > 100,])
combined <- combined[-16825,] ### 38851 peaks left

### inspect the distribution of the normalized peaks
ggplot(combined, aes(normalized_c1, normalized_c2)) + geom_point()
ggplot(combined, aes(normalized_c1, normalized_c3)) + geom_point()
ggplot(combined, aes(normalized_c2, normalized_c3)) + geom_point()

hist(combined$normalized_c1, breaks = 100)
hist(combined$normalized_c2, breaks = 100)
hist(combined$normalized_c3, breaks = 100)

### save combined peaks to csv 
write.csv(combined[,c(1:3,13:15,7:9)], '~/Desktop/HMH/rds/2021.09.10.wt_tem_peaks.normalized.csv')
saveRDS(combined, '~/Desktop/HMH/rds/2021.09.10.wt_tem_peaks.normalized.rds' )
