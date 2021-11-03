#!/bin/bash

path=~/Desktop/Sung_work/sam/ATAC_seq/macs2/_narrowPeak/

for input in naive wt_tcm ko_tem ko_tcm; do mkdir $path/homer/${input} ; done

for input in naive wt_tcm ko_tem ko_tcm; do perl /Users/sung/opt/miniconda3/share/homer-4.10-0/bin/findMotifsGenome.pl $path/${input}_peaks.narrowPeak mm10 /Users/sung/Desktop/Sung_work/sam/ATAC_seq/macs2/_narrowPeak/homer/${input} ; done 
