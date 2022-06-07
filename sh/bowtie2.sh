#!/bin/bash

path=.
crop=36

### Alignment of samples : bowtie2 for ATAC-seq
### bowtie2 index location : ~/Desktop/ref/mm10/mm10.index
### output: bam

for input in Icos_WT1 Icos_WT3 Icos_KO1 Icos_KO2 Icos_KO3; do mkdir $path/${input}_sam
bowtie2 --mm \
-p 10 \
--no-unal \
--non-deterministic \
--no-discordant \
-x ~/Desktop/ref/mm10/mm10.index \
-1 $path/${input}.${crop}.R1.paired.fastq.gz \
-2 $path/${input}.${crop}.R2.paired.fastq.gz \
-S $path/${input}_sam/${input}.${crop}.sam 2> $path/${input}_sam/${input}.${crop}.bowtie2.log ; done
