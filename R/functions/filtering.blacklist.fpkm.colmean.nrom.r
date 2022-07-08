## blacklist
bl <- read.table('~/Desktop/ref/blacklist/mm10.blacklist.bed')
colnames(bl) <- c('chr','start','end')
bl[1:3,]
bl %>% dim() ## 164 regions
bl$chr %>% table() ## mm10 blacklist only contains regions from chr1-19
