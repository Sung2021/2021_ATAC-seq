#!/bin/bash
path=.

samtools merge -o d05.wt.merge.bam d05.wt1.100.final.bam d05.wt2.100.final.bam d05.wt3.100.final.bam
samtools merge -o d05.ko.merge.bam d05.ko1.100.final.bam d05.ko2.100.final.bam d05.ko3.100.final.bam
samtools merge -o d12.wt.merge.bam d12.wt1.100.final.bam d12.wt2.100.final.bam d12.wt3.100.final.bam
samtools merge -o d12.ko.merge.bam d12.ko1.100.final.bam d12.ko2.100.final.bam d12.ko3.100.final.bam

