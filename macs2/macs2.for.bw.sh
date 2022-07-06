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

### no SPMR
### example

#!/bin/bash
path=.
mkdir $path/macs2.out
for i in *.bam;do mkdir $path/macs2.out/${i}
macs2 callpeak -t $path/${i} -f BAMPE -g mm -n ${i} -B -q 0.05 --outdir $path/macs2.out/${i} ;done
