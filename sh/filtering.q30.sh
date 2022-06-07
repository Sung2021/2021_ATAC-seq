#!/bin/bash

path=.
crop=36

### dropping read under q30 
### samtools 

for input in Icos_WT1 Icos_WT2 Icos_WT3 Icos_KO1 Icos_KO2 Icos_KO3 ; 
do samtools view -@ 10 -q 30 -b $path/${input}_sam/${input}.${crop}.sam > $path/${input}.${crop}.q30.bam ; done
