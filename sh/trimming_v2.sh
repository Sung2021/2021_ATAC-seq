#!/bin/bash

path=.

### Trimming to keep 36bp
## run
crop=36
day=d05

for input in wt1 wt2 wt3 ko1 ko2 ko3; do 
java -jar ~/Desktop/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 \
${day}.${input}.R1.fastq.gz ${day}.${input}.R2.fastq.gz \
$path/trimmed/${day}.${input}.${crop}.R1.paired.fastq.gz $path/trimmed/${day}.${input}.${crop}.R1.unpaired.fastq.gz \
$path/trimmed/${day}.${input}.${crop}.R2.paired.fastq.gz $path/trimmed/${day}.${input}.${crop}.R2.unpaired.fastq.gz \
CROP:${crop}; done

