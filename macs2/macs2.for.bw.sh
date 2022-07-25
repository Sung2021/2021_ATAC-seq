## call peaks by macs2
macs2 callpeak -t # input bam file 
               -c mock data # no antibody sample
               -f BAMPE # paired end bam)
               -g mm
               -n prefix ## name that represents the sample
               -B ## to store -log10 pvalue and qvalue in the output
               --SPMR # to save signal per million read
               -q 0.05 # qvalue cutoff
               -- outdir ## to designate output folder

### no SPMR
### example

#!/bin/bash
path=.
mkdir $path/macs2.out
for i in *.bam;do mkdir $path/macs2.out/${i}
macs2 callpeak -t $path/${i} -c -f BAMPE -g mm -n ${i} -B -q 0.05 --outdir $path/macs2.out/${i} ;done

#!/bin/bash
path=.
for i in WT1 WT2 WT3 KO1 KO2; do bedGraphToBigWig $path/Icos_${i}.final.bam.sorted.bam/Icos_${i}.final.bam.sorted.bam_treat_pileup.bdg \
~/Desktop/ref/mm10/mm10.chrom.sizes.txt ${i}.pileup.bw ; done 

## in R
library(rtracklayer)
bw <- import('~/Desktop/Sung_work/sam/ATAC_seq/Tle3/macs2_conditional/test.bw')
bw[1:3,]
bw$score/sum(coverage)*10^6
