suppressMessages(library(tidyverse))
suppressMessages(library(data.table))


af_check <- fread('variants_af/af_check.tsv',header=T)

#af_check <- af_check %>% filter(gnomADg_AF_NFE == '-')

af_check <- af_check %>% select(Uploaded_variation,PL_AF,gnomADg_AF_NFE) %>% 
  mutate(gnomADg_AF_NFE = ifelse(gnomADg_AF_NFE == '-', 0,af_check$gnomADg_AF_NFE)) %>%
  mutate(gnomADg_AF_NFE = as.numeric(af_check$gnomADg_AF_NFE)) %>%
  filter(PL_AF < 0.01)
sum(is.na(af_check$gnomADg_AF_NFE))

af_check %>% ggplot(aes(x=PL_AF,y=log10(gnomADg_AF_NFE))) + 
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis")  +
  theme_minimal()
