suppressMessages(library(tidyverse))
suppressMessages(library(gridExtra))

ped <- read.table('input/popinfo.tsv',sep='\t',header=T)
feno <- read.table('input/polincluded.tsv',sep='\t',header=T) %>% select(-Covid_severity)

ped <- ped %>% left_join(feno) %>% filter(Polish_Genome_included == 1) 
ped$Sex <- ifelse(ped$Sex == 1,'M','F')
ped <- na.omit(ped)

ped$range <- NA
ped$range <- ifelse(ped$Age < 10, '<10',ped$range)
ped$range <- ifelse(ped$Age < 20 & ped$Age > 10 , '10-20',ped$range)
ped$range <- ifelse(ped$Age >=20 & ped$Age < 30 , '20-30',ped$range)
ped$range <- ifelse(ped$Age >=30 & ped$Age < 40 , '30-40',ped$range)
ped$range <- ifelse(ped$Age >=40 & ped$Age < 50 , '40-50',ped$range)
ped$range <- ifelse(ped$Age >=50 & ped$Age < 60 , '50-60',ped$range)
ped$range <- ifelse(ped$Age >=60 & ped$Age < 70 , '60-70',ped$range)
ped$range <- ifelse(ped$Age >=70 & ped$Age < 80 , '70-80',ped$range)
ped$range <- ifelse(ped$Age >=80 & ped$Age <= 90 , '80-90',ped$range)
ped$range <- ifelse(ped$Age >90, '90+',ped$range)

ped %>% 
  ggplot(aes(x = range, fill = Sex)) + 
  geom_bar(data = filter(ped, Sex == "F"),col='black') + 
  geom_bar(data = filter(ped, Sex == "M"), aes(y = ..count.. * (-1)),col='black') + 
  scale_y_continuous(breaks = seq(-100, 100, 20), labels = abs(seq(-100, 100, 20))) + 
  xlab('Age (years)') + ylab('Number of individuals') + coord_flip() + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_fill_manual(values = c('#27384A','#B6B6B6'))

ggsave('age_pyramid.tiff')








