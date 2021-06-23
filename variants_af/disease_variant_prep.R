suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

print('Preparing AF list')
af <- fread('../input/multisample_20210519.dv.bcfnorm.filtered.AFlist.txt',showProgress=T) %>% filter(V3 > 0 & V4 > 0)
colnames(af) <- c('Uploaded_variation','PL_AF','PL_AN','PL_AC')

print('Loading VEP file')
# vep <- fread('/genom_polaka/input/multisample_20210519.dv.bcfnorm.filtered.vep.tsv.gz',
#              skip=101,showProgress=T)S

vep <- fread('/genom_polaka/input/multisample_20210519.dv.bcfnorm.filtered.vep_CLINVAR.tsv.gz',skip=120,showProgress=T) %>% filter(PICK == 1)

colnames(vep)[1] <- 'Uploaded_variation'

vep <- vep %>% left_join(af, by = 'Uploaded_variation') 

af_vep <- vep %>% select(PL_AF,PL_AC,IMPACT,Consequence,EXON,INTRON) 
af_vep <- af_vep %>%
  mutate(group = case_when(PL_AF > 0.005~ ">0.5%",
                          PL_AF >0 & PL_AF < 0.001 ~ "<0.1%" ,
                           !PL_AF < 0.001 & !PL_AF > 0.005 ~ "0.1-0.5%"
  ))

af_vep$group <- factor(af_vep$group)
af_vep$group <- ordered(af_vep$group, levels = c('<0.1%','0.1-0.5%','>0.5%'))


#af_vep$group <- ifelse(af_vep$PL_AF == 0, 'missing', af_vep$group)

af_list <- c('HIGH','MODERATE','LOW')
group.colors <- c(HIGH = "#27384A", MODERATE ="#48C095", LOW = "#B6B6B6")

af_vep$group <- factor(af_vep$group)
af_vep$group <- ordered(af_vep$group, levels = c("<0.1%", "0.1-0.5%", ">0.5%"))


af_vep %>% 
  filter(IMPACT != 'MODIFIER' & group != 'missing') %>%
  group_by(IMPACT,group) %>% 
  summarise(n=n()) %>% 
  write.table('../input/diseases/impact_stacked_ready.tsv',sep='\t',
              row.names = F,
              col.names = T,
              quote = F)
  # ggplot(aes(fill=IMPACT,y=n,x=group)) +
  # geom_bar(position="fill", stat="identity") +
  # theme_classic() +
  # scale_fill_manual(values = group.colors) +
  # xlab('Allele frequencies') +
  # ylab('% of variants') +
  # # ggtitle('% of variants with known impact') +
  # # theme(plot.title = element_text(hjust = 0.5),
  # #      plot.subtitle = element_text(hjust = 0.5)) +
  # ggsave('variants_af_files/figure-gfm/stacked_af_impact.png',dpi=320)

class.colors <- c(exonic = "#27384A", intronic ="#48C095", noncoding = "#B6B6B6")

consequence_vep <- af_vep %>% select(group,Consequence,EXON,INTRON,IMPACT)
consequence_vep$type <- NA
consequence_vep$type <- ifelse(consequence_vep$EXON != '-','Exonic',consequence_vep$type)
consequence_vep$type <- ifelse(consequence_vep$INTRON != '-','Intronic',consequence_vep$type)
consequence_vep$type <- ifelse(consequence_vep$INTRON == '-' & consequence_vep$EXON == '-',
                               'Noncoding',consequence_vep$type)

consequence_vep %>% 
  group_by(group,type,IMPACT,Consequence) %>% 
  summarise(n=n()) %>%
  write.table('../input/diseases/consequence_ready.tsv',sep='\t',
              row.names = F,
              col.names = T,
              quote = F)
  # ggplot(aes(fill=class,y=n,x=group)) + 
  # geom_bar(position="fill", stat="identity") + 
  # theme_classic() +
  # scale_fill_manual(values = class.colors) +
  # xlab('Allele frequencies') +
  # ylab('% of variants') +
  # # ggtitle('% of variants with known impact') +
  # # theme(plot.title = element_text(hjust = 0.5), 
  # #      plot.subtitle = element_text(hjust = 0.5)) +
  # ggsave('variants_af_files/figure-gfm/consequence.png',dpi=320)


verticals <- c(log10(0.001),log10(0.005))

af_vep %>% 
  #filter(group != 'missing' & PL_AF != 0) %>%
  ggplot(aes(log10(PL_AF))) +
  geom_histogram(fill='#48C095',col='#27384A', bins=30) +
  geom_vline(xintercept = verticals, linetype='dashed',col='#BC0020') +
  ylab('Variant count') + xlab('log10 AF') + theme_classic() +
  # ggtitle('Allele frequencies distribution') +
  # theme(plot.title = element_text(hjust = 0.5), 
  #       plot.subtitle = element_text(hjust = 0.5)) +
  ggsave("variants_af_files/figure-gfm/afhist.png",dpi=320)

af_vep %>% 
  #filter(group != 'missing') %>%
  ggplot(aes(PL_AC)) +
  geom_histogram(fill='#48C095',col='#27384A', bins=30) +
  #geom_vline(xintercept = verticals, linetype='dashed',col='#BC0020') +
  ylab('Variant count') + xlab('AC') + theme_classic() +
  # ggtitle('Allele counts distribution') +
  # theme(plot.title = element_text(hjust = 0.5), 
  #       plot.subtitle = element_text(hjust = 0.5)) +
  ggsave("variants_af_files/figure-gfm/ac_hist.png",dpi=320)


print('Filtering ACMG variants')
acmg_list <- read.table('../input/diseases/Ensembl_ACMG_v3.txt',sep='\t')
vep %>% filter(Gene %in% acmg_list$V1 & IMPACT == 'HIGH') %>% 
  write.table('../input/diseases/acmg_ready.tsv',sep='\t',
              row.names = F,
              col.names = T,
              quote = F)

print('Filtering putative variants')
vep %>% 
  filter(PL_AF < 0.001) %>% 
  filter(IMPACT == 'HIGH' | IMPACT == 'MODERATE') %>%
  write.table('../input/diseases/putative_ready.tsv',sep='\t',
              row.names = F,
              col.names = T,
              quote = F)

print('Filtering pathogenic variants')  
vep %>% filter(grepl('pathogenic',CLIN_SIG) & PL_AF < 0.001) %>% 
  write.table('../input/diseases/clin_sig_ready.tsv',sep='\t',
              row.names = F,
              col.names = T,
              quote = F)

### Filter by AF and IMPACT

print('Plotting AF per IMPACT')

vep %>% 
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
  ggsave('variants_af_files/figure-gfm/impact_all_af.png')

print('Filtering AF common summary')

vep %>% 
  filter(PL_AF > 0.005) %>%
  select(VARIANT_CLASS,IMPACT) %>%
  group_by(VARIANT_CLASS,IMPACT) %>%
  summarise(n=n()) %>% 
  mutate(AF = '>0.5%') %>%
  write.table('../output/common_summary.tsv', sep='\t', quote = F, row.names = F)

print('Filtering AF mediumrare summary')
vep %>% filter(PL_AF > 0.001 & PL_AF < 0.005) %>%
  select(VARIANT_CLASS,IMPACT) %>%
  group_by(VARIANT_CLASS,IMPACT) %>%
  summarise(n=n()) %>% 
  mutate(AF = '0.1-0.5%') %>%
  write.table('../output/mediumrare_summary.tsv', sep='\t', quote = F, row.names = F)

print('Filtering AF rare summary')
vep %>% filter(PL_AF < 0.001) %>%
  select(VARIANT_CLASS,IMPACT) %>%
  group_by(VARIANT_CLASS,IMPACT) %>%
  summarise(n=n()) %>% 
  mutate(AF = '<0.1%') %>%
  write.table('../output/rare_summary.tsv', sep='\t', quote = F, row.names = F)

# print('Filtering AF 0 summary')
# vep %>% filter(PL_AF == 0) %>%
#   select(VARIANT_CLASS,IMPACT) %>%
#   group_by(VARIANT_CLASS,IMPACT) %>%
#   summarise(n=n()) %>% 
#   mutate(AF = '0%') %>%
#   write.table('../output/missing_summary.tsv', sep='\t', quote = F, row.names = F)

gc()
print('Done')
