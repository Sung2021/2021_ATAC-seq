ggvenn(loci.all,stroke_size = 0.5, set_name_size = 4, 
       show_percentage = F, show_elements = F)

loci.ven <- ggvenn::list_to_data_frame(loci.all) %>% data.frame(row.names = 1)
loci.ven %>% head()
loci.overlap <- loci.ven[loci.ven %>% rowSums() == 3, ] %>% rownames()

loci.ven.tmp <- loci.ven[,c(1:2)]
loci.ven.tmp %>% head()
loci.ven.tmp[loci.ven.tmp %>% rowSums()== 2,] %>% rownames() %>% length()
loci1 <- loci.ven.tmp[loci.ven.tmp %>% rowSums()== 2,] %>% rownames()
loci2 <- setdiff(loci.3, loci1)
loci3 <- setdiff(loci.4, loci1)
setdiff(loci.3, loci1) %>% length()
setdiff(loci.4, loci1) %>% length()

loci.ven.tmp <- loci.ven[,c(2:3)]
loci.ven.tmp %>% head()
loci.ven.tmp[loci.ven.tmp %>% rowSums()== 2,] %>% rownames() %>% length()
loci4 <- loci.ven.tmp[loci.ven.tmp %>% rowSums()== 2,] %>% rownames()
loci5 <- setdiff(loci.4, loci1)
loci6 <- setdiff(loci.5, loci1)
setdiff(loci.4, loci4) %>% length()
setdiff(loci.5, loci4) %>% length()


loci.ven.tmp <- loci.ven[,c(1,3)]
loci.ven.tmp %>% head()
loci.ven.tmp[loci.ven.tmp %>% rowSums()== 2,] %>% rownames() %>% length()
loci7 <- loci.ven.tmp[loci.ven.tmp %>% rowSums()== 2,] %>% rownames()
loci8 <- setdiff(loci.3, loci7)
loci9 <- setdiff(loci.5, loci7)
setdiff(loci.3, loci7) %>% length()
setdiff(loci.5, loci7) %>% length()


ggvenn(list(loci.3=loci.all[[1]] ,
            loci.4=loci.all[[2]]),
       stroke_size = 0.5, set_name_size = 4, 
       show_percentage = F, show_elements = F)

ggvenn(list(loci.3=loci.all[[1]] ,
            loci.5=loci.all[[3]]),
       stroke_size = 0.5, set_name_size = 4, 
       show_percentage = F, show_elements = F)

ggvenn(list(loci.4=loci.all[[2]] ,
            loci.5=loci.all[[3]]),
       stroke_size = 0.5, set_name_size = 4, 
       show_percentage = F, show_elements = F)


loci1 <- loci.ven.tmp[loci.ven.tmp %>% rowSums()== 2,] %>% rownames()
loci4 <- loci.ven.tmp[loci.ven.tmp %>% rowSums()== 2,] %>% rownames()
loci7 <- loci.ven.tmp[loci.ven.tmp %>% rowSums()== 2,] %>% rownames()
loci.overlap

loci9 <- setdiff(loci1,loci.overlap)
loci10 <- setdiff(loci4,loci.overlap)
loci11 <- setdiff(loci7,loci.overlap)


loci9 %>% length()
loci10 %>% length()
loci11 %>% length()

setdiff(loci.3, union(loci1,loci7)) %>% length()
setdiff(loci.4, union(loci1,loci4)) %>% length()
setdiff(loci.5, union(loci4,loci7)) %>% length()

loci12 <- setdiff(loci.3, union(loci1,loci7))
loci13 <- setdiff(loci.4, union(loci1,loci4))
loci14 <- setdiff(loci.5, union(loci4,loci7))

loci.venn.all <- list(l1=loci1,
                      l2=loci2,
                      l3=loci3,
                      l4=loci4,
                      l5=loci5,
                      l6=loci6,
                      l7=loci7,
                      l8=loci8,
                      l9=loci9,
                      l10=loci10,
                      l11=loci11,
                      l12=loci12,
                      l13=loci13,
                      l14=loci14)
saveRDS(loci.venn.all, 'ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/log2fc.1.loci.info.rds')
norm.peak %>% head()
norm.peak[,'venn'] <- 'not_significant'
norm.peak[loci9,'venn'] <- 'loci9'
norm.peak[loci10,'venn'] <- 'loci10'
norm.peak[loci11,'venn'] <- 'loci11'
norm.peak[loci12,'venn'] <- 'loci12'
norm.peak[loci13,'venn'] <- 'loci13'
norm.peak[loci14,'venn'] <- 'loci14'
norm.peak %>% head()
norm.peak$venn %>% table()
write.csv(norm.peak, 'ATAC_seq/atac_seq.readvalue.and.normalized_21.11.29/log2fc.1.loci.info.csv')
