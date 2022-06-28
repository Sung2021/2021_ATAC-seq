#!/bin/bash 

samtools sort -@ 10 ICOS_KO.merge.bam -o ICOS_KO.sorted.bam
bedtools genomecov -ibam ICOS_KO.sorted.bam -bg > ICOS_KO.merge.bedgraph
sort -k1,1 -k2,2n ICOS_KO.merge.bedgraph > ICOS_KO.sorted.bedgraph
bedGraphToBigWig ICOS_KO.sorted.bedgraph ~/Desktop/ref/mm10/mm10.chrom.sizes.txt ICOS_KO.bw

