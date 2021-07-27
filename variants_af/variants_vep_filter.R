suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

print('Preparing AF list')
af <- fread('../input/multisample_20210716.dv.bcfnorm.filt.unrelated.vcf.gz.AFlist.txt',showProgress=T) %>% filter(V3 > 0 & V4 > 0)
colnames(af) <- c('Uploaded_variation','PL_AF','PL_AN','PL_AC')

print('Loading VEP file')

vep <- fread('/genom_polaka/input/multisample_20210716.dv.bcfnorm.filt.unrelated.tsv.gz',skip=113,showProgress=T)

colnames(vep)[1] <- 'Uploaded_variation'


muko <- vep %>% 
  filter(Uploaded_variation == 'chr7_117559590_ATCT_A') %>%
  left_join(af,by='Uploaded_variation')
muko %>% write.table('muko_vep.tsv',col.names = T, row.names = F,quote = F,sep='\t')

nbs <- vep %>% 
  filter(Uploaded_variation == 'chr8_89971213_ATTTGT_A') %>%
  left_join(af,by='Uploaded_variation')
nbs %>% write.table('nbs.tsv',col.names = T, row.names = F,quote = F,sep='\t')

