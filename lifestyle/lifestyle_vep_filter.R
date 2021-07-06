suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

print('Preparing AF list')
af <- fread('../input/multisample_20210519.dv.bcfnorm.filtered.AFlist.txt',showProgress=T) %>% filter(V3 > 0 & V4 > 0)
colnames(af) <- c('Uploaded_variation','PL_AF','PL_AN','PL_AC')

print('Loading VEP file')

vep <- fread('/genom_polaka/input/multisample_20210519.dv.bcfnorm.filtered.vep_CLINVAR.tsv.gz',skip=120,showProgress=T) %>% filter(PICK == 1)

colnames(vep)[1] <- 'Uploaded_variation'

vep <- vep %>% left_join(af, by = 'Uploaded_variation') 

maoa <- vep %>% filter(SYMBOL == 'MAOA' | SYMBOL == 'DRD4') %>% 
  left_join(af,by='Uploaded_variation')
maoa %>% write.table('dopamine_genes.tsv',col.names = T,    row.names = F,quote = F,sep='\t')

lifestyle <- read.csv('lifestyle_vep_pos.csv',header=T,sep='\t')

vep %>% filter(Existing_variation %in% lifestyle$Uploaded_variation) %>% 
  #left_join(af,by='Uploaded_variation') %>% 
  write.table('lifestyle_genes.tsv',col.names = T, row.names = F,quote = F,sep='\t')


alko_pos <- read.csv('alko_rsid_pos.csv',header=T,sep='\t')
vep %>% filter(Existing_variation %in% alko_pos$Uploaded_variation) %>% 
  write.table('2alko_rsid_vep_out.tsv',col.names = T,
                                                       row.names = F,quote = F,sep='\t')

alko_genes <- read.table('alko_genes.txt',header=T)
vep %>% filter(SYMBOL %in% alko_genes$SYMBOL) %>% write.table('alko_genes_vep_out.tsv',col.names = T,
                                                       row.names = F,quote = F,sep='\t')

