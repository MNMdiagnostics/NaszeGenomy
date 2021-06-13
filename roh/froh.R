library(detectRUNS)
library(data.table)


mapFile <- fread('input/roh_data.map')
mapFile <- as.data.frame(mapFile)

runs <- 'output/roh/roh.txt'
runs <- readExternalRuns(inputFile = runs, program = 'BCFtools')

froh <- Froh_inbreeding(runs = runs, mapFile = 'input/roh_data.map')
froh_chr <- Froh_inbreeding(runs = runs, mapFile = 'input/roh_data.map', genome_wide = F)



write.table(froh,'output/roh/froh.output',sep='\t')

plot_Runs(runs = runs)
