suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

print('Preparing AF list')
af <- fread('../input/multisample_20210519.dv.bcfnorm.filtered.AFlist.txt',showProgress=T) %>% filter(V3 > 0 & V4 > 0)
colnames(af) <- c('Uploaded_variation','PL_AF','PL_AN','PL_AC')

print('Loading VEP file')

vep <- fread('/genom_polaka/input/multisample_20210519.dv.bcfnorm.filtered.vep_CLINVAR.tsv.gz',skip=120,showProgress=T) %>% filter(PICK == 1)

colnames(vep)[1] <- 'Uploaded_variation'


maoa <- vep %>% 
  filter(SYMBOL == 'CFTR' | Location == 'chr7:117498313-117519392') %>%
  left_join(af,by='Uploaded_variation')
maoa %>% write.table('CFTR.tsv',col.names = T, row.names = F,quote = F,sep='\t')

