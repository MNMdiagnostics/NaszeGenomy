library(tidyverse)
library(data.table)

imiss <-  fread('output/plink_stats/plink.imiss')
misshap <- fread('output/plink_stats/plink.missing.hap')
varmiss <-  fread('output/plink_stats/plink.lmiss')


imiss %>% 
  ggplot(aes(F_MISS)) +geom_histogram(fill='#48C095',col='#27384A') +
  theme_classic() + ylab('Samples') + xlab('Missing call rate')

varmiss %>% filter(F_MISS > 0.005) %>% 
  ggplot(aes(F_MISS)) +geom_histogram(fill='#48C095',col='#27384A') +
  theme_classic() + xlim(c(0.005,1)) + ggtitle('N variants with missing call rate > 0.5%')

varmiss %>%
  ggplot(aes(F_MISS)) + 
  geom_boxplot() + 
  facet_wrap(~factor(CHR), nrow = 5,ncol = 5)


varmiss_NA <- varmiss %>% filter(is.na(F_MISS))
