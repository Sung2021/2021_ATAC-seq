#############################################################################
#############################################################################

library(DiffBind)
library(tidyverse)
library(rtracklayer)

setwd('~/Desktop/Sung_work/sam/ATAC_seq/')
list.files('./bam/', pattern = c('ko_tcm','rmd.bam$'))[c(1,3,5)]

### sample names
pattern1 ='wt_tcm'
pattern2 ='ko_tcm'
sample_raw <- c(list.files('./bam/', pattern = c('wt_tcm','rmd.bam$'))[c(1,3,5)],
                list.files('./bam/', pattern = c('ko_tcm','rmd.bam$'))[c(1,3,5)])
sample_raw2 <- c(list.files('./bam/macs2/narrowPeak/', pattern = pattern1), 
                 list.files('./bam/macs2/narrowPeak/', pattern = pattern2))
sample_raw
sample_raw2

### generate sample sheet
### "/Users/sung/Desktop/Sung_work/sam/ATAC_seq" samples.csv
samples <- data.frame(matrix(nrow=6,ncol=8)) ## 2 samples, 8 metadata
names(samples) <- c('SampleID','cell_type','mutation', 'Factor',
                    'Replicate','bamReads','Peaks','PeakCaller')
samples$SampleID <- substr(sample_raw,1,7)
samples$cell_type <- substr(sample_raw, 4,6)
samples$mutation <- substr(sample_raw, 1,2)
### Factor : the criteria to compare groups
samples$Factor <- samples$mutation
samples$Replicate <- c(1:3,1:3)
samples$bamReads <- paste0('bam/', sample_raw)
samples$Peaks <- paste0('bam/macs2/narrowPeak/', sample_raw2)
samples$PeakCaller <- 'bed'
samples
write.csv(samples, 'samples.csv')

DEP <- dba(sampleSheet = samples)
DEP <- dba.count(DEP)
DEP
plot(DEP)
saveRDS(DEP, 'rds/2021.10.11.Tcm_wt_ko.DEP.rds')
info <- dba.show(DEP)
write.csv(info,'rds/2021.10.11.Tcm_wt_ko.DEP.info.csv')


DEP.con <- dba.contrast(DEP) ### count the reads
DEP.con <- dba.analyze(DEP.con) ### This is DESeq2 analysis
dba.show(DEP.con, bContrasts = T)


DEP.con.db <- dba.report(DEP.con)
db <- data.frame(DEP.con.db)
db <- db %>% arrange(seqnames, start)
db[1:3,]
write.csv(db, 'DEP.con.db.csv')

dba.plotVenn(DEP.con, contrast = 1,bDB=TRUE,bGain=TRUE, bLoss=TRUE, bAll=FALSE)
dba.plotPCA(DEP.con)
dba.plotVolcano(DEP.con)
dba.plotBox(DEP.con)
dba.plotHeatmap(DEP.con, contrast = 1, correlations = FALSE, scale="row")


peak_reads <- data.frame(matrix(nrow=nrow(db), ncol=3))
colnames(peak_reads) <- c('Chr','Start','End')
peak_reads$Chr <- db$seqnames
peak_reads$Start <- db$start
peak_reads$End <- db$end

gr = makeGRangesFromDataFrame(peak_reads, keep.extra.columns=T)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPseeker)
peakAnno.peak_reads <- annotatePeak(gr, tssRegion=c(-3000, 3000), 
                                   TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                   annoDb="org.Mm.eg.db")
peakAnno.peak_reads@anno
peakAnno.peak_reads@detailGenomicAnnotation
peakAnno.peak_reads@annoStat
peakAnno.peak_reads@peakNum


anno.db <- data.frame(peakAnno.peak_reads@anno)
anno.db[1:3,]
saveRDS(anno.db, 'rds/2021.10.11.Tcm_wt_ko.DEP.anno.db.rds')

##### anno.db statistics ######
anno.db$feature <- 'NA'
anno.db[grep('Promoter', anno.db$annotation),'feature'] <- 'Promoter'
anno.db[grep('Intron', anno.db$annotation),'feature'] <- 'Intron'
anno.db[grep('Exon', anno.db$annotation),'feature'] <- 'Exon'
anno.db[grep('3\' UTR', anno.db$annotation),'feature'] <- '3\' UTR'
anno.db[grep('5\' UTR', anno.db$annotation),'feature'] <- '5\' UTR'
anno.db[grep('Downstream', anno.db$annotation),'feature'] <- 'Downstream'
anno.db[grep('Distal Intergenic', anno.db$annotation),'feature'] <- 'Distal Intergenic'
anno.db$feature <- factor(anno.db$feature)
## bar plot-feature
ggplot(anno.db, aes(feature, fill=feature)) + geom_bar() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
df <- data.frame(table(anno.db$feature))
df$perc <- df$Freq/sum(df$Freq)*100
## pie plot
ggplot(df, aes(x="", y=perc, fill=Var1)) + geom_col() + coord_polar(theta = 'y') +
  geom_text(aes(label = round(perc, digits = 1)),
            position = position_stack(vjust = 0.5))
## bar plot-chromosome
ggplot(anno.db, aes(geneChr, fill=geneChr)) + geom_bar() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

table(anno.db$geneChr)
hist((anno.db$distanceToTSS), breaks = 500)

ggplot(anno.db, aes(feature, distanceToTSS)) + geom_jitter(size=0.5, alpha=0.5) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
