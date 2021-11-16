##############################################################
##############################################################
norm.peaks <- readRDS('rds/2021.09.20.atac_seq.normalized.peaks.rds')
rownames(norm.peaks) <- paste0('peak_', rownames(norm.peaks))
norm.peaks <- norm.peaks[,c(grep('naive',colnames(norm.peaks)), 
                            grep('wt_tcm', colnames(norm.peaks)))]

raw.peaks <- readRDS('rds/2021.09.20.atac_seq.combined.peaks.rds')
rownames(raw.peaks) <- paste0('peak_', rownames(raw.peaks))
raw.peaks[1:3,][,c(4:length(colnames(raw.peaks)))]
raw.peaks[1:3,]
### count matrix
### info sheet
count.mtx <- raw.peaks[,c(4:length(colnames(raw.peaks)))]
count.mtx <- count.mtx[,c(grep('wt_tem',colnames(count.mtx)), 
                          grep('wt_tcm', colnames(count.mtx)))]
count.mtx[1:3,]

info <- data.frame(matrix(nrow = length(colnames(count.mtx)), ncol = 3))
colnames(info) <- c('sample', 'cell_type','mutation')
info$sample <- colnames(count.mtx)
info$cell_type <- substr(info$sample, 4,6)
info$mutation <- substr(info$sample, 1,2)
info

info$cell_type <- factor(info$cell_type, levels = c('tem','tcm'))

### DESeq object 
dds <- DESeqDataSetFromMatrix(count.mtx, info, ~cell_type)
dim(dds)
dds <- DESeq(dds)
res <- results(dds)
dim(res)
res <- data.frame(res)
res[1:3,]

###################################################
### volcano plot
ggplot(res, aes(log2FoldChange, -log10(padj))) + 
  geom_point(size=0.5, alpha=0.5) +
  theme +
  xlab('Log2(fold_changes in WT Tcm/WT Tem)') +
  geom_vline(xintercept = c(-1,1),col='dark green', linetype=1) +
  geom_hline(yintercept = -log10(0.05), col='dark green', linetype=1)
  
### when slicing, 'NA' are generated and are included in the dataset
### to avoid 'NA's, using dplyr piping 
dim(res) ## total peaks
res.inc <- res %>% filter(log2FoldChange >= log2(2)) %>% filter(padj < 0.05)
res.des <- res %>% filter(log2FoldChange <= -log2(2)) %>% filter(padj < 0.05)
res.sig <- rbind(res.inc, res.des)
dim(res.sig)
dim(res.inc)
dim(res.des)

count.mtx[rownames(res[res$log2FoldChange >= log2(2) & res$padj < 0.05,][1:3,]),]
norm.peaks[rownames(res[res$log2FoldChange >= log2(2) & res$padj < 0.05,][1:3,]),]



