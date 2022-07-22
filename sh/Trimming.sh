#!/bin/bash

path=.

### Trimming to keep 36bp
## run
crop=36

for input in Icos_WT2 Icos_WT3 Icos_KO1 Icos_KO2 Icos_KO3; do 
java -jar ~/Desktop/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 \
${input}_R1_001.fastq.gz ${input}_R2_001.fastq.gz \
$path/${input}.${crop}.R1.paired.fastq.gz $path/${input}.${crop}.R1.unpaired.fastq.gz \
$path/${input}.${crop}.R2.paired.fastq.gz $path/${input}.${crop}.R2.unpaired.fastq.gz \
CROP:${crop}; done
