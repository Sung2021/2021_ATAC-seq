#!/bin/bash

path=~/Desktop/Sung_work/sam
cd ~/Desktop/Sung_work/sam

### merge bams from one condition (combining replicates in one condition)
for input in naive ko_tem ; do samtools merge $path/${input}.bam $path/${input}1.36.q30.sorted2.rmd.bam $path/${input}2.36.q30.sorted2.rmd.bam ; done

for input in wt_tem wt_tcm ko_tcm ; do samtools merge $path/${input}.bam $path/${input}1.36.q30.sorted2.rmd.bam $path/${input}2.36.q30.sorted2.rmd.bam $path/${input}3.36.q30.sorted2.rmd.bam ; done

### macs2 callpeak with cutoff-analysis output
for input in naive wt_tem wt_tcm ko_tem ko_tcm; do macs2 callpeak --cutoff-analysis -t $path/${input}.bam -f BAMPE -g mm -n ${input} \
-B -q 0.05 --SPMR --outdir $path/macs2; done




#########################################################################
#########################################################################
####### codes before 2021.09.13

#!/bin/bash


input=wt_tem
path=~/Desktop/Sung_work/sam

for j in 1 2 3; do macs2 callpeak -t $path/${input}${j}.36.q30.sorted2.rmd.bam -f BAMPE -g mm -n ${input}${j} -B -q 0.05 --SPMR --outdir $path/macs2; done

##### generate merged.bed as the first input for bedtools coverage
cd $path/macs2/
cat *Peak | sort -k1,1 -k2,2n | bedtools merge -i stdin | awk 'BEGIN{OFS="\t"} $0=$0"\t"NR' > ${input}.merged.bed

##### generate each paired end bed files from bam files
##### this is also the second input file for bedtools coverage
for k in 1 2 3; do bedtools bamtobed -i $path/${input}${k}.36.q30.sorted2.rmd.bam > $path/${input}${k}.pe.bed; done

##### bedtools coverage: how much of the genome does my data cover?
##### generate each coverage bed files comparing merged.bed and each PE bed files
for m in 1 2 3; do bedtools coverage -a $path/macs2/merged.bed -b $path/${input}${m}.pe.bed > $path/${input}${m}_coverage.bed; done


##### we only need column 5 in coverage.bed for coverage information (cut -f 5)
cd ~/Desktop/Sung_work/sam
for i in *coverage.bed; do cat ${i} | cut -f5 > ${i}.coverage;done
