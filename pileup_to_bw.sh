#!/bin/bash

##############
##### this step can be done after macs2 peak calling 
##### this step use pileup files from macs2 peak calling
##### only for visualization in IGV 
##### replicates are handled in for loop for the convenience

path=~/Desktop/Sung_work/sam/macs2
## set input name of the replicates
input=wt_tem

## sort the pileup files (from macs2 output) by chromosomal order
for i in 1 2 3; do sort -k1,1 -k2,2n $path/${input}${i}_treat_pileup.bdg > $path/${input}${i}.sorted.bdg ; done

## covert the sorted pileup files to bw for the IGV use
for i in 1 2 3; do bedGraphToBigWig $path/${input}${i}.sorted.bdg ~/Desktop/software/mm10.chrom.sizes.txt $path/${input}${i}.sorted.bdg.bw ; done

## add delete intermediate files command below


#############################################
######## updated 2021.09.13 ############

#!/bin/bash

path=~/Desktop/Sung_work/sam

##for input in naive1 naive2 wt_tem1 wt_tem2 wt_tem3 wt_tcm1 wt_tcm2 wt_tcm3 ko_tem1 ko_tem2 ko_tcm1 ko_tcm2 ko_tcm3 ; do echo ${input}; cat ${input}_coverage.bed | awk '{if($4 =='8459') print $0}' ; done

## sort the pileup files (from macs2 output) by chromosomal order
#for input in naive1 naive2 wt_tem1 wt_tem2 wt_tem3 wt_tcm1 wt_tcm2 wt_tcm3 ko_tem1 ko_tem2 ko_tcm1 ko_tcm2 ko_tcm3; do sort -k1,1 -k2,2n $path/macs2/${input}_treat_pileup.bdg > $path/macs2/${input}.sorted.bdg ; done

## covert the sorted pileup files to bw for the IGV use
for input in naive1 naive2 wt_tem1 wt_tem2 wt_tem3 wt_tcm1 wt_tcm2 wt_tcm3 ko_tem1 ko_tem2 ko_tcm1 ko_tcm2 ko_tcm3; do bedGraphToBigWig $path/macs2/${input}.sorted.bdg ~/Desktop/software/mm10.chrom.sizes.txt $path/macs2/${input}.sorted.bdg.bw ; done
