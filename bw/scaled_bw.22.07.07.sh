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

## scaling bw in R
library(rtracklayer)
bw <- import('ATAC_seq/Icos.2022.06/bw.22.07.06/WT1.pileup.bw')
bw$score <- (bw$score/(sum of total raw counts in each library)*10^6
export(bw, 'scaled.bw)

## sum of total raw counts in each library
raw counts calculated by coverage between each bed and merged(union) peaks
merge peaks :filtered by fold change > 4 (column 7)
and merged peaks from WT and KO

cat ICOS_WT_merged_peaks.narrowPeak | awk '{if ($7 >= 4) print }' > ICOS_WT.bed
cat ICOS_KO_merged_peaks.narrowPeak | awk '{if ($7 >= 4) print }' > ICOS_KO.bed

bedtools merge -i ICOS_WT.bed -i ICOS_KO.bed > merge.bed

## calculate the coverage
for input in *.pe.bed; do bedtools coverage -a $path/merge.bed -b $path/${input} > $path/${input}.cov.bed;done

cov.bed files : raw reads of ATAC-seq
#############################################
