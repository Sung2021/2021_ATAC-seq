#!/bin/bash

path=.
mkdir $path/homer
for i in *.bed; do mkdir $path/homer/${i}
perl /Users/sung/opt/miniconda3/share/homer-4.10-0/bin/findMotifsGenome.pl \
$path/${i} mm10 $path/homer/${i}; done
