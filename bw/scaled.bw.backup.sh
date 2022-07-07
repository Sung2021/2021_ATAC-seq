for i in WT1 WT2 WT3 KO1 KO2; do macs2 callpeak -t Icos_${i}.final.bam  -c -f BAMPE -g mm -n ${i} -B -q 0.05  --outdir macs2/${i} ;done 

for i in WT1 WT2 WT3 KO1 KO2; do sort -k1,1 -k2,2n macs2/${i}/${i}_treat_pileup.bdg > bdg/${i}.sorted.bdg ; done

for i in WT1 WT2 WT3 KO1 KO2; do bedGraphToBigWig bdg/${i}.sorted.bdg ~/Desktop/software/mm10.chrom.sizes.txt bw/${i}.pileup.treat.bw ; done

samtools merge -o ICOS_KO.merge.bam Icos_KO1.final.bam Icos_KO2.final.bam
samtools merge -o ICOS_WT.merge.bam Icos_WT1.final.bam Icos_WT2.final.bam Icos_WT3.final.bam

for i in WT KO; do mkdir macs2/${i}_merged/ ;done 
for i in WT KO; do macs2 callpeak -t merged_bam/Icos_${i}.merge.bam  -c -f BAMPE -g mm -n ${i} -B -q 0.05  --outdir macs2/${i}_merged/ ;done 

for i in WT KO; do cat macs2/${i}_merged/${i}_peaks.narrowPeak | awk '{if ($7 >= 4) print }' > ${i}.fc.filtered.bed ; done

bedtools merge -i ICOS_WT.bed -i ICOS_KO.bed > merge.bed

mkdir pe.bed
for i in WT1 WT2 WT3 KO1 KO2; do bedtools bamtobed -i Icos_${i}.final.bam > pe.bed/${i}.pe.bed ;done

mkdir cov.bed
for i in WT1 WT2 WT3 KO1 KO2; do bedtools coverage -a merge.bed -b pe.bed/${i}.pe.bed > cov.bed/${i}.cov.bed ;done
