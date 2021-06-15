library(detectRUNS)
library(tidyverse)

runs <- readExternalRuns(inputFile = 'output/roh/roh.txt',program='BCFtools')
runs_froh <- runs %>% filter(chrom != 'chrX' & chrom != 'chrY')

froh <- Froh_inbreeding(runs_froh,mapFile = 'input/roh_data.map')
write.table('output/roh/froh.output',sep='\t', quote = F, row.names = F)

#froh_wide <- Froh_inbreeding(runs,mapFile = 'input/roh_data.map', genome_wide=F)


# consecutiveRuns <- consecutiveRUNS.run(
#   genotypeFile ='input/roh_data.ped',
#   mapFile = 'input/roh_data.map',
#   minSNP = 20,
#   ROHet = FALSE,
#   maxGap = 10^6,
#   minLengthBps = 25000,
#   maxOppRun = 1,
#   maxMissRun = 1
# )
# 
# 
# write.table(consecutiveRuns,'output/roh/consec_runs.tsv', sep='\t', quote = F, row.names = F)