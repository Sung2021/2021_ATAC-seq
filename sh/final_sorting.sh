for input in  Icos_WT1 Icos_WT2 Icos_WT3 Icos_KO1 Icos_KO2 Icos_KO3; do samtools sort -@ 10 -n ${input}.36.q30.sorted.rmd.bam -o ${input}.final.bam; done
