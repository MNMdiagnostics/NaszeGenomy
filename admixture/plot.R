library(tidyverse)
source("barNaming.R")
K <- 6
name <- "Eu_thin_LD"
tblfile <- paste(name, K, "Q", sep =".")
famfile <- paste(name, "fam", sep =".")

tbl <- read.table(tblfile)
indTable<-read.table(famfile, col.names=c("Pop","Sample","c1","c2","c3","c4"))
merged=cbind(tbl, indTable)
ordered=merged[order(merged$Pop), ]

pdfname <- paste(name, K, "pdf", sep = ".")
pdf(file=pdfname, 20, 4)
mp <- barplot(t(as.matrix(ordered[, 1:K])), col=rainbow(K),border=0, space=0, axes=F, axisname=F)
text(mp, par("usr")[3], labels = barNaming(ordered$Pop), srt = 90, adj = c(1.1,1.1), xpd = TRUE, cex=.14)
dev.off()

