library(tidyverse)
library(data.table)

af <- fread('input/AFlist.txt')

af_single <- af %>% filter(V2-V3==0)
af_multi <- af %>% filter(V2-V3!=0)

af_single$Location <- paste(af_single$V1,af_single$V2,sep=':')

af_multi$Location <- paste(af_multi$V1,af_multi$V2,sep=':')
af_multi$Location <- paste(af_multi$Location,af_multi$V3,sep='-')

new_af <- rbind(af_single,af_multi) %>% select(Location,V4)
colnames(new_af) <- c('Location','PL_AF')

acmg <- read.csv('input/diseases/acmg_filtered.tsv',skip=101,sep='\t')
acmg <- acmg %>% left_join(new_af, by = 'Location')

putative <- read.csv('input/diseases/putative_filtered.tsv',skip=101,header=T,sep='\t')
putative <- putative %>% left_join(new_af, by = 'Location')

clin_sig <- read.csv('input/diseases/clinsig_filtered.tsv',skip=101,header=T,sep='\t')
clin_sig <- clin_sig %>% left_join(new_af, by = 'Location')

write.table(acmg, 'input/diseases/acmg_ready.tsv',sep='\t',
            row.names = F,
            col.names = T,
            quote = F)

write.table(putative, 'input/diseases/putative_ready.tsv',sep='\t',
            row.names = F,
            col.names = T,
            quote = F)

write.table(clin_sig, 'input/diseases/clin_sig_ready.tsv',sep='\t',
            row.names = F,
            col.names = T,
            quote = F)


### Filter by AF and IMPACT

common <- fread('input/diseases/common.tsv',skip=101)
common <- common %>% left_join(new_af, by = 'Location')

common %>% 
  select(VARIANT_CLASS,IMPACT) %>%
  group_by(VARIANT_CLASS,IMPACT) %>%
  summarise(n=n()) %>% 
  mutate(AF = '>0.5%') %>%
  write.table('output/common_summary.tsv', sep='\t', quote = F, row.names = F)

mediumrare <- fread('input/diseases/medium_rare.tsv',skip=101)
mediumrare <- mediumrare %>% left_join(new_af, by = 'Location')

mediumrare %>% 
  select(VARIANT_CLASS,IMPACT) %>%
  group_by(VARIANT_CLASS,IMPACT) %>%
  summarise(n=n()) %>% 
  mutate(AF = '0.1%-0.5%') %>%
  write.table('output/mediumrare_summary.tsv', sep='\t', quote = F, row.names = F)

rare <- fread('input/diseases/ultra_rare.tsv',skip=101)
rare <- rare %>% left_join(new_af, by = 'Location')

rare %>% 
  select(VARIANT_CLASS,IMPACT) %>%
  group_by(VARIANT_CLASS,IMPACT) %>%
  summarise(n=n()) %>% 
  mutate(AF = '<0.1%') %>%
  write.table('output/rare_summary.tsv', sep='\t', quote = F, row.names = F)

library(gridExtra)

impact_all_af <- rbind(common,mediumrare,rare) %>% 
  select(IMPACT, PL_AF) %>% 
  ggplot(aes(log(PL_AF))) +
  geom_histogram(fill='#48C095',col='#27384A',bins=30) +
  ylab('log Allele Frequency') +
  xlab('Variant impact') +
  ggtitle('log Allele frequency in functional categories')+
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5)) +
  facet_wrap(~IMPACT,nrow=4,scales = 'free') +
  ggsave('diseases/diseases_files/figure-gfm/impact_all_af.png')
