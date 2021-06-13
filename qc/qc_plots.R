suppressMessages(library(tidyverse))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))
suppressMessages(library(data.table))

varmiss <-  fread('../output/plink_stats/plink.lmiss')
varmiss %>% 
  ggplot(aes(log10(F_MISS))) + 
  geom_histogram(fill='#48C095',col='#27384A',bins=30) +
  geom_vline(xintercept = log10(0.1), col='red',linetype='dashed') +
  theme_classic() + 
  ggtitle('Variants missing call rate') + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab('log10 Call rate missing') +
  ggsave("qc_files/figure-gfm/var_miss.png",dpi=320)


varmiss <-  fread('../output/plink_stats/multisample_filtered.lmiss')
varmiss %>% 
  ggplot(aes(log10(F_MISS))) + 
  geom_histogram(fill='#48C095',col='#27384A',bins=30) +
  theme_classic() + 
  ggtitle('Filtered variants with missing call rate < 10%') + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab('log10 Call rate missing') +
  ggsave("qc_files/figure-gfm/filtered_var_miss.png",dpi=320)

imiss <- fread('../output/plink_stats/plink.imiss') %>% as.data.frame()
colnames(imiss)[2] <- 'sample'
imiss %>% ggplot(aes(F_MISS)) + 
  geom_histogram(fill='#48C095',col='#27384A') +
  xlab('Call rate missing') + 
  ylab('count') +
  ggtitle('Missing call-rate per genome') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))  +
  ggsave("qc_files/figure-gfm/samples_miss.png",dpi=320)

imiss <- fread('../output/plink_stats/multisample_filtered.imiss') %>% as.data.frame()
colnames(imiss)[2] <- 'sample'
imiss %>% ggplot(aes(F_MISS)) + 
  geom_histogram(fill='#48C095',col='#27384A') +
  xlab('Call rate missing') + 
  ylab('count') +
  ggtitle('Filtered missing call-rate per genome < 10%') +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))  +
  ggsave("qc_files/figure-gfm/filtered_samples_miss.png",dpi=320)

gc()