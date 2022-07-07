## 
final.bam 

## call peak by macs2
macs2 callpeak -t final.bam \
               -c no input \ ## no input uses genomic background
               -f BAMPE \
               -g mm \
               -n name \
               -B \
               -q 0.05 
               
macs2 callpeak -t $path/${i} -c -f BAMPE -g mm -n ${i} -B -q 0.05  $path/${i}_macs2 

## among outputs, use treat_pileup.bdg file be used
## sort by chromosomal order (optional?)
sort -k1,1 -k2,2n $path/${input}${i}_treat_pileup.bdg > $path/${input}${i}.sorted.bdg

## covert the sorted pileup files to bw for the IGV use
bedGraphToBigWig $path/${input}${i}.sorted.bdg ~/Desktop/software/mm10.chrom.sizes.txt $path/${input}${i}.sorted.bdg.bw 

## output
pileup.treat.bw  ## with the score column


## sum of total raw counts in each library
raw counts calculated by coverage between each bed and merged(union) peaks
merge peaks :filtered by fold change > 4 (column 7)
and merged peaks from WT and KO

samtools merge -o ICOS_KO.merge.bam Icos_KO1.final.bam Icos_KO2.final.bam
samtools merge -o ICOS_WT.merge.bam Icos_WT1.final.bam Icos_WT2.final.bam Icos_WT3.final.bam
 
cat ICOS_WT_merged_peaks.narrowPeak | awk '{if ($7 >= 4) print }' > ICOS_WT.bed
cat ICOS_KO_merged_peaks.narrowPeak | awk '{if ($7 >= 4) print }' > ICOS_KO.bed

bedtools merge -i ICOS_WT.bed -i ICOS_KO.bed > merge.bed

## output 
merge.bed

## each sample cov.bed compared to merge.bed
## bam to bed first

mkdir pe.bed
for i in WT1 WT2 WT3 KO1 KO2; do bedtools bamtobed -i Icos_${i}.final.bam > pe.bed/${i}.pe.bed ;done


## calculate the coverage
for input in *.pe.bed; do bedtools coverage -a $path/merge.bed -b $path/${input} > $path/${input}.cov.bed;done

cov.bed files : raw reads of ATAC-seq

#####################################
## read cov.bed in r to generated raw count mtx
library(ggplot2)
library(dplyr)
library(DESeq2)
library(GenomicRanges)
library(rtracklayer)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

theme_set(theme_bw())
## read input merged bed from macs2 call
dir <- 'ATAC_seq/Icos.2022.06/cov.bed.22.07.07'
list.files(dir, full.names = T, pattern = 'cov.bed')
## Shaoqi's method using union peaks from merged narrowPeak
## coverage output : 'chr','start','end','read','frg_len','ref_len','coverage'
### when reading coverage bed
samples <- c(paste0('Icos_','KO',1:2),
             paste0('Icos_','WT',1:3))
tmp <- read.table(list.files(dir, full.names = T, pattern = 'cov.bed')[1])
for(i in 1:length(list.files(dir, full.names = T))){
  bed <- read.table(list.files(dir, full.names = T)[i])
  colnames(bed) <- c('chr','start','end','read','frg_len','ref_len','coverage')
  bed$sample <- samples[i] 
  rownames(bed) <- paste0(bed$chr,'_',
                          bed$start,'_',
                          bed$end)
  assign(samples[i],bed)
}

bed.list <- list(Icos_WT1,Icos_WT2,Icos_WT3,
     Icos_KO1,Icos_KO2)

bed.tmp <- data.frame(matrix(nrow = nrow(bed.list[[1]]),
                             ncol = 5))
bed.tmp <- cbind(bed.list[[1]][,c(1:3)], bed.tmp)
rownames(bed.tmp) <- rownames(bed.list[[1]])
colnames(bed.tmp)[4:ncol(bed.tmp)] <- c(paste0('WT',1:3), paste0('KO',1:2))

for(i in 1:5){
  bed.tmp[,3+i] <- bed.list[[i]][,4]
}

bed.tmp %>% dim()
bed.tmp %>% write.csv('ATAC_seq/Icos.2022.06/Icos.ATAC_seq.raw.read.from.bed.22.07.07.csv')
## bed.tmp is the raw count mtx

## scaling factor : total raw count of each library
scaling_factor <- bed.tmp[,4:ncol(bed.tmp)] %>% colSums() %>% as.vector()

   WT1      WT2      WT3      KO1      KO2 
10327503 10517737 10698878 11423816  8870570 

## scaling bw in R
scaling_factor <- bed.tmp[,4:ncol(bed.tmp)] %>% colSums() %>% as.vector()

library(rtracklayer)
dir <- 'ATAC_seq/Icos.2022.06/bw.22.07.07'
files <-list.files(dir, full.names = T, pattern = 'bw')
files <- files[c(3:5,1:2)] ## reorder the file names
for(i in 1:5){
  names <- c(paste0('WT',1:3), paste0('KO',1:2))
  bw <- import(files[i])
  bw$score <- (bw$score/scaling_factor[i])*10^6
  export(bw, paste0(dir,'/',names[i],'.scaled.bw'))
}

## outputs : scaled.bw
