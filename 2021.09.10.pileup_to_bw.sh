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
