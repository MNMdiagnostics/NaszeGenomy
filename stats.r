suppressMessages(library(tidyverse))

flagstat_path = input/multiqc_data/multiqc_samtools_flagstat.txt
depth_path = input/depth_concat.txt

flagstat <- read.table(flagstat_path,header=T,sep='\t')
flagstat <- flagstat %>% select(Sample,total_passed,properly.paired_passed,flagstat_total,mapped_passed) %>%
    mutate(total_passed = total_passed/1e+06, properly.paired_passed = properly.paired_passed/1e+06, flagstat_total = flagstat_total/1e+06,
    mapped_passed = mapped_passed/1e+06)

depth <- read.table(depth_path,header=T,sep='\t')

stats <- flagstat %>% left_join(depth, by= c('Sample'='sample'))

x <- stats %>% summarise(total_mean = mean(total_passed), total_sd = sd(total_passed), paired_mean = mean(properly.paired_passed), 
    paired_sd = sd(properly.paired_passed), depth_mean = mean(average_depth), depth_sd = sd(average_depth))

write.table(x, file = "output/library_stats.tsv",  sep = "\t")     
print('library_stats.tsv saved')

