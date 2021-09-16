#!/bin/bash

##############
##### this step can be done after macs2 peak calling 

##### path to the location which bw files are stored
path=~/Desktop/Sung_work/sam/macs2/tmp

##### multibigwigsummary 
##### input bigwig files : *.bdg.bw 
##### naive wt_tem wt_tcm ko_tem ko_tcm : total 5 kinds

multiBigwigSummary bins -b $(ls *.bdg.bw) -o results.npz --outRawCounts scores_per_bin.tab
