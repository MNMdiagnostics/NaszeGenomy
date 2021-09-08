suppressMessages(library(tidyverse))
suppressMessages(library(data.table))


af_check <- fread('af_check.tsv',header=T) %>% 
  mutate(gnomADg_AF_NFE = as.numeric(gnomADg_AF_NFE)) %>%
  mutate(gnomADg_AF_NFE = ifelse(is.na(gnomADg_AF_NFE)==T, 1e-5, gnomADg_AF_NFE))



af_check <- af_check %>% select(Uploaded_variation,PL_AF,gnomADg_AF_NFE) %>% 
  filter(PL_AF < 0.01)
  
sum(is.na(af_check$gnomADg_AF_NFE))

af_check %>% ggplot(aes(x=PL_AF,y=gnomADg_AF_NFE)) + 
  #geom_bin2d(bins = 70) +
  #scale_fill_continuous(type = "viridis")  +
  stat_binhex()
  geom_smooth(formula= y ~ x, method='lm') +
  ggsave('af_check.png')
