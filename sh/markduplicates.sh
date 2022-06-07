#!/bin/bash

path=.
crop=36

#picard markduplicates
for input in Icos_WT1 Icos_WT2 Icos_WT3 Icos_KO1 Icos_KO2 Icos_KO3 ; do 
java -jar ~/Desktop/software/picard.jar MarkDuplicates \
I= $path/${input}.${crop}.q30.sorted.bam  \
O= $path/${input}.${crop}.q30.sorted.rmd.bam \
M= $path/${input}.${crop}.q30.sorted.rmd.metrics.txt \
ASSUME_SORTED=TRUE  REMOVE_DUPLICATES=true CREATE_INDEX=true QUIET=true ; done
