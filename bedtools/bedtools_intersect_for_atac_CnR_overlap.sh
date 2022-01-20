 ##### -a control -b test 
 ##### bedtools intersect 
 ##### -wo : including both peak info and overlapping info
 
 bedtools intersect -a atac.peak.sig.8069.bed -b CnR.less.dyn.peak.bed -wo > atac.CnR.wo.bed
 bedtools intersect -a atac.peak.sig.8069.bed -b CnR.less.dyn.peak.bed -wo | wc -l  
