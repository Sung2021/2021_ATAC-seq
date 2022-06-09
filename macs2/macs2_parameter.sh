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

### example

#!/bin/bash

path=.

mkdir $path/macs2
macs2 callpeak -t $path/merged.bam -f BAMPE -g mm -n merged -B -q 0.05 --SPMR --outdir $path/macs2

