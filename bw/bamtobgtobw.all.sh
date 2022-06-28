#!/bin/bash

path=.
mkdir bw
for input in *.rmd.bam ; do samtools sort -@ 10 ${input} -o ${input}.sorted.bam
bedtools genomecov -ibam ${input}.sorted.bam -bg > ${input}.bedgraph
sort -k1,1 -k2,2n ${input}.bedgraph > ${input}.sorted.bedgraph
bedGraphToBigWig ${input}.sorted.bedgraph ~/Desktop/ref/mm10/mm10.chrom.sizes.txt ${input}.bw; done
