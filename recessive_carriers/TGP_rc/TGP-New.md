## Libraries and options



```R
library(data.table)
library(glue)
library(dplyr)
library(plotly)
library(ggplot2)
library(grid)
library(gridExtra)
options(scipen=999)
```

## Load and prepare data


```R
# Read input VCF
tgp_vcf <- fread("/home/matdawidziuk/tgp_pl/multisample_20210519.dv.bcfnorm.filtered.nogt.vcf.gz", skip="#CHROM", sep="\t")
```


```R
# Check dim
dim(tgp_vcf)
head(tgp_vcf)
```


    Error in eval(expr, envir, enclos): object 'tgp_vcf' not found
    Traceback:




```R
# Read TGP vep annotated file
tgp_vep <- fread("/data/NGS/annotations/test_design/pl_genomes/run_1/output/multisample_20210519.dv.bcfnorm.filtered.nogt.header.format.dummygt_VEP.tsv", sep="\t")
```


```R
# Check dim
dim(tgp_vep)
head(tgp_vep)
```


```R
#  Split INFO
tgp_vcf[, c("AF", "AQ", "AN", "AC") := tstrsplit(INFO, ";", fixed=TRUE)]
```


```R
# Clean new columns, remove AQ (duplicates QUAL)
tgp_vcf[, AF := gsub("[^0-9.-]", "", AF)][, AN := gsub("[^0-9.-]", "", AN)][, AC := gsub("[^0-9.-]", "", AC)][, AQ := NULL]
```


```R
head(tgp_vcf)
```


```R
# Create key - further used to merge with VCF file
tgp_vep[, key2 := gsub(">", "_", gsub(":", "_", key))]
setkey(tgp_vep, key2)
setkey(tgp_vcf, ID)
```


```R
# Copy selected columns from VCF to TSV
selected_columns <- c("REF", "ALT", "QUAL", "AF", "AN", "AC")
tgp_vep[tgp_vcf, (selected_columns) := mget(paste0('i.', selected_columns))]
head(tgp_vep)
```


```R
#save(tgp_vep, file="~/tgp_vep_2021_06_26.RData")
```


```R
dim(tgp_vep)
```


<ol class=list-inline>
	<li>37318044</li>
	<li>154</li>
</ol>




```R
# Read AR genes from OMIM
omim_recessive <- fread("/srv/data/TGP-notebooks/recessive_genes.bed", sep="\t", col.names=c("chr", "start", "end", "symbol"))
omim_recessive <- unique(omim_recessive, by = "symbol")
```


```R
tgp_recessive <- tgp_vep[which((Symbol %in% omim_recessive$symbol) & AC!=0)]
```


```R
tgp_vep_nonsyn_rare <- tgp_recessive[which( gnomAD_AF < 0.01 & AF < 0.01 & Impact %in% c("MODERATE", "HIGH"))]

```


```R
#save(tgp_vep_nonsyn_rare,file="~/tgp_vep_nonsyn_rare_2021_06_27.RData" )
#save(tgp_recessive,file="~/tgp_recessive_2021_06_27.RData" )
```


```R
gnomad_vep <- fread("/data/NGS/annotations/test_design/pl_genomes/run_1/output/gnomad_genome_v3.0_gt_subset-gatk-haplotype.vcf.gz.annotated.tsv", skip="#Uploaded_variation")
```


```R
gnomad_vep[, gnomAD_AF := as.numeric(gnomad_vep$gnomad_genome.vcf.gz_AF)]
gnomad_vep[, gnomAD_AF_NFE := as.numeric(gnomad_vep$gnomad_genome.vcf.gz_AF_nfe)]
gnomad_vep[, gnomAD_AF:=ifelse(is.na(gnomAD_AF),0,gnomAD_AF)]
gnomad_vep[, gnomAD_AF_NFE:=ifelse(is.na(gnomAD_AF_NFE),0,gnomAD_AF_NFE)]
```

    Warning message in eval(jsub, SDenv, parent.frame()):
    “NAs introduced by coercion”


```R
gnomad_nonsyn_rare <- gnomad_vep[which( gnomAD_AF < 0.01 & IMPACT %in% c("MODERATE", "HIGH"))]
```


```R
idx_duplicated <- which(duplicated(gnomad_nonsyn_rare$key))
gnomad_nonsyn_rare <- gnomad_nonsyn_rare[-idx_duplicated,]
dim(gnomad_nonsyn_rare)
#save(gnomad_nonsyn_rare,file="~/gnomad_nonsyn_rare_2021_06_27.RData" )
```


```R
#colnames(gnomad_nonsyn_rare)
```

## Start from here


```R
#load("~/tgp_vep_nonsyn_rare_2021_06_27.RData")
load("/home/jupyter-tomgambin/gnomad_nonsyn_rare_2021_06_27.RData" )
load("/home/jupyter-tomgambin/tgp_recessive_2021_06_27.RData" )
tgp_vep_nonsyn_rare <- tgp_recessive[which( gnomAD_AF < 0.01 & Impact %in% c("MODERATE", "HIGH"))]
```


```R
nrow(tgp_vep_nonsyn_rare)
nrow(gnomad_nonsyn_rare)
```


26288



730362



```R
omim_recessive <- fread("/srv/data/TGP-notebooks/recessive_genes.bed", sep="\t", col.names=c("chr", "start", "end", "symbol"))
omim_recessive <- unique(omim_recessive, by = "symbol")
nrow(omim_recessive)
```


2512


## Download PanelApp data


```R

#install.packages("rvest")
library(rvest)
panels <- read_html("https://panelapp.genomicsengland.co.uk/panels/")
panels2 <- panels %>%html_nodes("h4")
allPanelsList <- (lapply(panels2[-length(panels2)], function(panel){
    panelURL <- paste0("https://panelapp.genomicsengland.co.uk",panel%>% html_nodes("a") %>%  html_attr("href") ,"download/01234/")
    print(panelURL)
    fread(panelURL)
}))
allPanelsDT <- rbindlist(allPanelsList)
idx <- grep ("Expert Review Green", allPanelsDT$"Sources(; separated)" )
allPanelsDT2 <- allPanelsDT[idx,]
allPanelsDT3 <- allPanelsDT2[grep("BIALLELIC", allPanelsDT2$Model_Of_Inheritance),]
```

    [1] "https://panelapp.genomicsengland.co.uk/panels/399/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/933/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/929/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/930/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/934/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/931/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/932/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/774/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/540/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/245/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/391/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/511/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/269/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/502/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/263/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/510/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/34/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/134/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/258/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/477/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/139/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/260/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/657/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/282/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/557/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/543/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/38/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/208/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/545/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/166/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/90/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/55/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/13/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/234/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/842/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/843/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/749/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/230/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/214/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/286/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/109/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/491/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/147/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/745/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/847/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/259/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/221/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/544/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/30/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/197/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/81/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/210/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/64/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/244/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/517/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/507/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/508/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/145/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/25/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/512/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/308/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/31/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/207/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/232/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/225/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/250/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/658/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/111/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/168/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/5/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/560/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/283/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/487/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/519/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/159/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/484/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/251/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/293/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/26/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/652/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/47/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/9/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/235/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/209/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/265/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/192/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/553/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/136/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/53/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/648/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/271/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/562/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/119/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/554/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/91/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/314/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/158/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/50/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/305/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/152/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/110/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/7/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/23/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/63/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/6/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/772/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/480/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/312/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/522/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/394/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/11/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/212/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/318/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/200/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/290/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/552/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/167/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/478/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/144/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/33/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/61/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/132/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/412/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/402/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/201/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/227/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/254/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/249/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/528/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/720/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/473/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/59/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/407/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/99/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/115/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/126/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/20/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/466/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/488/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/157/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/123/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/85/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/846/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/165/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/567/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/568/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/78/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/310/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/179/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/267/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/236/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/49/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/481/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/92/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/650/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/482/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/490/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/555/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/277/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/467/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/176/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/246/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/175/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/503/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/171/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/143/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/524/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/97/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/649/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/504/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/525/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/523/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/521/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/42/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/174/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/285/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/514/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/315/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/515/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/131/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/239/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/173/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/248/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/213/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/549/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/530/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/238/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/384/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/185/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/546/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/527/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/76/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/218/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/529/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/96/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/133/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/83/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/112/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/534/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/535/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/536/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/537/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/538/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/533/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/532/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/472/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/18/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/564/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/75/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/87/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/36/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/211/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/19/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/558/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/385/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/149/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/474/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/183/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/255/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/724/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/736/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/465/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/526/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/219/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/513/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/296/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/253/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/189/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/128/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/294/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/722/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/186/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/196/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/943/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/117/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/486/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/479/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/79/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/288/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/215/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/556/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/386/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/86/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/39/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/541/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/24/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/60/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/114/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/94/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/559/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/483/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/311/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/105/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/653/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/539/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/178/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/398/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/65/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/155/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/566/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/506/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/17/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/106/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/193/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/531/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/247/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/518/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/565/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/150/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/728/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/48/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/8/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/154/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/725/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/902/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/903/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/292/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/550/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/307/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/600/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/66/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/217/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/734/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/98/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/130/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/262/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/228/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/162/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/62/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/921/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/224/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/199/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/726/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/309/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/229/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/542/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/3/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/180/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/509/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/841/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/45/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/551/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/82/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/700/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/1/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/122/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/945/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/516/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/421/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/548/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/243/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/195/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/302/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/216/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/156/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/678/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/273/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/101/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/563/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/222/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/928/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/579/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/476/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/496/download/01234/"
    [1] "https://panelapp.genomicsengland.co.uk/panels/77/download/01234/"



```R
gnomad_nonsyn_rare[,key:=paste0(Location, ":", REF_ALLELE, ">", Allele)]
```

## Calculate cumulative frequencies


```R
calculateCumFreq <- function(tgp_clinvar, gnomad_clinvar){
    
    tgp_clinvar[, c("cum_AF", "nr_vars"):=list(sum(as.numeric(AF)), .N), by=Symbol]
    tgp_clinvar_unique <- tgp_clinvar[-which(duplicated(tgp_clinvar[,"Symbol",with=F])),]
    tgp_clinvar_unique <- tgp_clinvar_unique[order(cum_AF, decreasing =TRUE), c("Symbol", "cum_AF",  "nr_vars"),with=F]
    gnomad_clinvar[which(gnomad_clinvar$gnomAD_AF_NFE != 0), c("cum_gnomAD_AF_NFE","gnomAD_nr_vars_NFE")  := list(sum(as.numeric(gnomAD_AF_NFE)), .N), by=SYMBOL]
    gnomad_clinvar_unique <- gnomad_clinvar[-which(duplicated(gnomad_clinvar[,"SYMBOL",with=F])),]
    gnomad_clinvar_unique <- gnomad_clinvar_unique[order(cum_gnomAD_AF_NFE, decreasing =TRUE), c("SYMBOL", "cum_gnomAD_AF_NFE", "gnomAD_nr_vars_NFE"),with=F]
    
    
    #gnomad_clinvar[, c("cum_gnomAD_AF","gnomAD_nr_vars") := list(sum(as.numeric(gnomAD_AF)), .N), by=SYMBOL]
    #gnomad_clinvar[, cum_gnomAD_AF2 := cum_gnomAD_AF^2]
    #gnomad_clinvar[, cum_gnomAD_AF2_NFE := cum_gnomAD_AF_NFE^2]
    #tgp_clinvar[, cum_AF2 := cum_AF^2]

    list(tgp_clinvar_unique,gnomad_clinvar_unique )
}


tgp_clinvar <-  tgp_vep_nonsyn_rare[which(CLNSIG %in% c("Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic"))]
gnomad_clinvar  <- gnomad_nonsyn_rare[which(clinvar.vcf.gz_CLNSIG %in% c("Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic"))]
gnomad_clinvar <- gnomad_clinvar[which(gnomad_clinvar$gnomAD_AF_NFE != 0),]
rr <- calculateCumFreq (tgp_clinvar, gnomad_clinvar)
tgp_clinvar_unique <- rr[[1]]
gnomad_clinvar_unique <- rr[[2]]                      


tgp_clinvar_unique$cum_AC <- round(tgp_clinvar_unique$cum_AF*max(as.numeric(tgp_vep_nonsyn_rare$AN)))
tgp_clinvar_unique$cum_AN <- max(as.numeric(tgp_vep_nonsyn_rare$AN))#round(max(as.numeric(tgp_vep_nonsyn_rare$AN))*tgp_clinvar_unique$cum_AC)
tgp_clinvar_unique$cum_gnomad_AF_NFE <- gnomad_clinvar_unique$cum_gnomAD_AF_NFE[match(tgp_clinvar_unique$Symbol,gnomad_clinvar_unique$SYMBOL )]
tgp_clinvar_unique$gnomad_nr_vars_NFE <- gnomad_clinvar_unique$gnomAD_nr_vars_NFE[match(tgp_clinvar_unique$Symbol,gnomad_clinvar_unique$SYMBOL )]
nfe_cnt <- 2*32299
tgp_clinvar_unique$cum_gnomad_AC_NFE <- round(tgp_clinvar_unique$cum_gnomad_AF_NFE*nfe_cnt)
tgp_clinvar_unique$cum_gnomad_AN_NFE <- nfe_cnt
tgp_clinvar_unique$cum_AF_fold_change <- log2(tgp_clinvar_unique$cum_AF/tgp_clinvar_unique$cum_gnomad_AF_NFE)

                     
suppressWarnings(pvals <- sapply(1:nrow(tgp_clinvar_unique), function(i){
        #print(i)
        row <- tgp_clinvar_unique[i,]
        if(!is.finite(row$cum_gnomad_AC_NFE)){return (NA)}
        M <- as.table(rbind(c(row$cum_AC, row$cum_AN), c(row$cum_gnomad_AC_NFE, row$cum_gnomad_AN_NFE)))
        dimnames(M) <- list(pop = c("PL", "NFE"),  value = c("AC", "AN"))
      (Xsq <- chisq.test(M))  # Prints test summary
    Xsq$p.value
}))
tgp_clinvar_unique$pvals <- pvals*nrow(tgp_clinvar_unique) #  Bonferonni correction
```


```R
summary(tgp_clinvar_unique$cum_AF_fold_change)
```


       Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
    -3.1363 -0.1086  0.8486  0.9852  2.0095  5.2297       9 



```R

```


```R

```


```R

```


```R
tgp_clinvar_unique[which(tgp_clinvar_unique$pvals < 0.01),]
```


<table>
<thead><tr><th scope=col>Symbol</th><th scope=col>cum_AF</th><th scope=col>nr_vars</th><th scope=col>cum_AC</th><th scope=col>cum_AN</th><th scope=col>cum_gnomad_AF_NFE</th><th scope=col>gnomad_nr_vars_NFE</th><th scope=col>cum_gnomad_AC_NFE</th><th scope=col>cum_gnomad_AN_NFE</th><th scope=col>cum_AF_fold_change</th><th scope=col>pvals</th></tr></thead>
<tbody>
	<tr><td>C2                                                        </td><td>0.023329800                                               </td><td>1                                                         </td><td>44                                                        </td><td>1886                                                      </td><td>0.0075581600                                              </td><td> 1                                                        </td><td>488                                                       </td><td>64598                                                     </td><td>1.626067                                                  </td><td>0.00000000010413876209229218366564776703444704548939814259</td></tr>
	<tr><td>TGM5                                                      </td><td>0.016436923                                               </td><td>2                                                         </td><td>31                                                        </td><td>1886                                                      </td><td>0.0040888476                                              </td><td> 5                                                        </td><td>264                                                       </td><td>64598                                                     </td><td>2.007174                                                  </td><td>0.00000000000640893529886687046770499652539629744686655588</td></tr>
	<tr><td>C8B                                                       </td><td>0.008483563                                               </td><td>2                                                         </td><td>16                                                        </td><td>1886                                                      </td><td>0.0024491506                                              </td><td> 3                                                        </td><td>158                                                       </td><td>64598                                                     </td><td>1.792389                                                  </td><td>0.00072884167742414231790787360054650889651384204626083374</td></tr>
	<tr><td>CBS                                                       </td><td>0.008055240                                               </td><td>1                                                         </td><td>15                                                        </td><td>1886                                                      </td><td>0.0002146807                                              </td><td> 5                                                        </td><td> 14                                                       </td><td>64598                                                     </td><td>5.229663                                                  </td><td>0.00000000000000000000000000000000000000000000000009082388</td></tr>
	<tr><td>OTOG                                                      </td><td>0.007953343                                               </td><td>2                                                         </td><td>15                                                        </td><td>1886                                                      </td><td>0.0015798706                                              </td><td>14                                                        </td><td>102                                                       </td><td>64598                                                     </td><td>2.331755                                                  </td><td>0.00000026090673382056531751232455150490974205013117170893</td></tr>
	<tr><td>GP6                                                       </td><td>0.007953340                                               </td><td>1                                                         </td><td>15                                                        </td><td>1886                                                      </td><td>0.0019975296                                              </td><td> 4                                                        </td><td>129                                                       </td><td>64598                                                     </td><td>1.993344                                                  </td><td>0.00008962331566893800270136305474721893915557302534580231</td></tr>
	<tr><td>MASP1                                                     </td><td>0.006362670                                               </td><td>1                                                         </td><td>12                                                        </td><td>1886                                                      </td><td>0.0005730213                                              </td><td> 3                                                        </td><td> 37                                                       </td><td>64598                                                     </td><td>3.472972                                                  </td><td>0.00000000000000197856810733732156561274713487032410770742</td></tr>
	<tr><td>C19orf12                                                  </td><td>0.004783850                                               </td><td>2                                                         </td><td> 9                                                        </td><td>1886                                                      </td><td>0.0002478660                                              </td><td> 4                                                        </td><td> 16                                                       </td><td>64598                                                     </td><td>4.270540                                                  </td><td>0.00000000000000000360632921425608510510921365516159829915</td></tr>
	<tr><td>SORD                                                      </td><td>0.004772010                                               </td><td>2                                                         </td><td> 9                                                        </td><td>1886                                                      </td><td>0.0005419470                                              </td><td> 1                                                        </td><td> 35                                                       </td><td>64598                                                     </td><td>3.138373                                                  </td><td>0.00000002358205094023312235399078215514201417590811615810</td></tr>
	<tr><td>ALOX12B                                                   </td><td>0.004772000                                               </td><td>2                                                         </td><td> 9                                                        </td><td>1886                                                      </td><td>0.0009913307                                              </td><td>30                                                        </td><td> 64                                                       </td><td>64598                                                     </td><td>2.267156                                                  </td><td>0.00286763538781958492868140453424530278425663709640502930</td></tr>
	<tr><td>FKBP14                                                    </td><td>0.004772000                                               </td><td>1                                                         </td><td> 9                                                        </td><td>1886                                                      </td><td>0.0008674370                                              </td><td> 1                                                        </td><td> 56                                                       </td><td>64598                                                     </td><td>2.459763                                                  </td><td>0.00032772621626388887634481750232851027249125763773918152</td></tr>
	<tr><td>PIGV                                                      </td><td>0.004241786                                               </td><td>3                                                         </td><td> 8                                                        </td><td>1886                                                      </td><td>0.0004646057                                              </td><td> 6                                                        </td><td> 30                                                       </td><td>64598                                                     </td><td>3.190593                                                  </td><td>0.00000017767759029256624497213913876220159338004123128485</td></tr>
	<tr><td>ABCB11                                                    </td><td>0.003711566                                               </td><td>4                                                         </td><td> 7                                                        </td><td>1886                                                      </td><td>0.0004028471                                              </td><td>12                                                        </td><td> 26                                                       </td><td>64598                                                     </td><td>3.203724                                                  </td><td>0.00000269795976183464327036024867079166966732373111881316</td></tr>
	<tr><td>XPA                                                       </td><td>0.003711563                                               </td><td>3                                                         </td><td> 7                                                        </td><td>1886                                                      </td><td>0.0003717554                                              </td><td> 8                                                        </td><td> 24                                                       </td><td>64598                                                     </td><td>3.319601                                                  </td><td>0.00000060004871020621404435605762694994069761378341354430</td></tr>
	<tr><td>ENAM                                                      </td><td>0.003711560                                               </td><td>2                                                         </td><td> 7                                                        </td><td>1886                                                      </td><td>0.0006815260                                              </td><td> 2                                                        </td><td> 44                                                       </td><td>64598                                                     </td><td>2.445185                                                  </td><td>0.00980380524443690698488396151333290617913007736206054688</td></tr>
	<tr><td>PROP1                                                     </td><td>0.003711560                                               </td><td>2                                                         </td><td> 7                                                        </td><td>1886                                                      </td><td>0.0006198380                                              </td><td> 5                                                        </td><td> 40                                                       </td><td>64598                                                     </td><td>2.582063                                                  </td><td>0.00274243891321816887426954245654542319243773818016052246</td></tr>
	<tr><td>HPGD                                                      </td><td>0.003181333                                               </td><td>2                                                         </td><td> 6                                                        </td><td>1886                                                      </td><td>0.0002788313                                              </td><td> 2                                                        </td><td> 18                                                       </td><td>64598                                                     </td><td>3.512167                                                  </td><td>0.00000154206306322945463381025085064512936128267028834671</td></tr>
	<tr><td>SLC26A3                                                   </td><td>0.003181333                                               </td><td>2                                                         </td><td> 6                                                        </td><td>1886                                                      </td><td>0.0003717904                                              </td><td>11                                                        </td><td> 24                                                       </td><td>64598                                                     </td><td>3.097070                                                  </td><td>0.00015502483707905190583008214755267317741527222096920013</td></tr>
	<tr><td>SLC19A1                                                   </td><td>0.002653930                                               </td><td>1                                                         </td><td> 5                                                        </td><td>1886                                                      </td><td>0.0003098180                                              </td><td> 1                                                        </td><td> 20                                                       </td><td>64598                                                     </td><td>3.098637                                                  </td><td>0.00238272711099553249208815586257514951284974813461303711</td></tr>
	<tr><td>AP4B1                                                     </td><td>0.002651110                                               </td><td>1                                                         </td><td> 5                                                        </td><td>1886                                                      </td><td>0.0003097648                                              </td><td> 6                                                        </td><td> 20                                                       </td><td>64598                                                     </td><td>3.097351                                                  </td><td>0.00238272711099553249208815586257514951284974813461303711</td></tr>
	<tr><td>STAG3                                                     </td><td>0.002120890                                               </td><td>1                                                         </td><td> 4                                                        </td><td>1886                                                      </td><td>0.0000929219                                              </td><td> 5                                                        </td><td>  6                                                       </td><td>64598                                                     </td><td>4.512507                                                  </td><td>0.00000043749047735085456468692830141187233294886027579196</td></tr>
	<tr><td>SBF2                                                      </td><td>0.001590670                                               </td><td>1                                                         </td><td> 3                                                        </td><td>1886                                                      </td><td>0.0000619819                                              </td><td> 3                                                        </td><td>  4                                                       </td><td>64598                                                     </td><td>4.681644                                                  </td><td>0.00007693537835650435155031984235662889659579377621412277</td></tr>
</tbody>
</table>




```R
fwrite(tgp_clinvar_unique[tgp_clinvar_unique$pvals < 0.01,], "Supp_Table_tgp_clinvar_unique.xls", row.names=F)
```


```R
#tgp_cadd <-  tgp_vep_nonsyn_rare[which(CADD_phred > 20)]
#gnomad_cadd  <- gnomad_nonsyn_rare[which(CADD_phred > 20)]
#rr1 <- calculateCumFreq (tgp_cadd, gnomad_cadd)
#tgp_cadd_unique <- rr1[[1]]
#gnomad_cadd_unique <- rr1[[2]] 

#tgp_clinvar_unique$cum_gnomad_AF2 <- gnomad_clinvar_unique$cum_gnomAD_AF2[match(tgp_clinvar_unique$Symbol,gnomad_clinvar_unique$SYMBOL )]
#tgp_clinvar_unique$cum_gnomad_AF2_NFE <- gnomad_clinvar_unique$cum_gnomAD_AF2_NFE[match(tgp_clinvar_unique$Symbol,gnomad_clinvar_unique$SYMBOL )]
#tgp_clinvar_unique$cum_gnomad_AF <- gnomad_clinvar_unique$cum_gnomAD_AF[match(tgp_clinvar_unique$Symbol,gnomad_clinvar_unique$SYMBOL )]
#tgp_clinvar_unique$gnomad_nr_vars <- gnomad_clinvar_unique$gnomAD_nr_vars[match(tgp_clinvar_unique$Symbol,gnomad_clinvar_unique$SYMBOL )]

#tgp_clinvar_unique$cum_gnomad_AN_NFE[which(!is.finite(tgp_clinvar_unique$cum_gnomad_AN_NFE))] <- 0 


#tgp_cadd_unique$cum_gnomad_AF2 <- gnomad_cadd_unique$cum_gnomAD_AF2[match(tgp_cadd_unique$Symbol,gnomad_cadd_unique$SYMBOL )]
#tgp_cadd_unique$cum_gnomad_AF2_NFE <- gnomad_cadd_unique$cum_gnomAD_AF2_NFE[match(tgp_cadd_unique$Symbol,gnomad_cadd_unique$SYMBOL )]

```


```R
tcu <- tgp_clinvar_unique[which(abs(tgp_clinvar_unique$cum_AF_fold_change) >=1  & tgp_clinvar_unique$pvals < 0.01  & tgp_clinvar_unique$nr_vars >=3),]
nrow(tcu)
#tcu[order(tcu$cum_AF_fold_change, decreasing=F),]

 #SLC26A2, ADSL, CNGB3, PCDH15, FANCD2, GPR179, PIGV and ABCB11.
tcu <- tgp_clinvar_unique[tgp_clinvar_unique$pvals < 0.01,]
#tcu
```


3



```R

```


```R

```

### Cum frequency differences
* Difference of cumulative frequencies between Polish and GnomAD (first column) and Polish and NFE  (second columns)
* Cumulative frequencies have been calculated for pathogenic and likely pathogenic clinvar variants (first row) and for predicted deleterious (CADD >20) variants (second row)



```R
library(ggrepel)
```


```R
set.seed(42)
p <- ggplot(tcu, aes(cum_AF, cum_gnomad_AF_NFE)) +
    geom_point() +
    # geom_text(aes(label = Symbol),hjust = 0, nudge_x = 0.0001)+ 
    #geom_abline() + #xlim(0, 0.008) + ylim(0, 0.003) +
    ggtitle("Cumulative allele frequencies of selected ClinVar variants in AR genes")+
    ylab("Per-gene cumulative allele frequencies in GnomAD NFE population")+
    xlab("Per-gene cumulative allele frequencies in Polish population")+
    scale_x_log10() + scale_y_log10()

p<- p +   geom_label_repel(
  aes(
    fill = "red", 
    label = Symbol),
  fontface = 'bold.italic', 
    color = 'white',
  size = 3,
     segment.color = "black",
    max.overlaps = 10,
    min.segment.length = 0,
  box.padding = unit(0.1, "lines"),
 point.padding = unit(0.1, "lines")
)

p +  theme(legend.position = "none") + theme(plot.title = element_text(size=13))
```


    
![png](output_44_0.png)
    



```R
allPanelsDT3$cum_AF <- tgp_clinvar_unique$cum_AF[match(allPanelsDT3$"Gene Symbol",tgp_clinvar_unique$Symbol )]
#allPanelsDT3$cum_AC <- tgp_clinvar_unique$cum_AC[match(allPanelsDT3$"Gene Symbol",tgp_clinvar_unique$Symbol )]
#allPanelsDT3$cum_AN <- tgp_clinvar_unique$cum_AN[match(allPanelsDT3$"Gene Symbol",tgp_clinvar_unique$Symbol )]
allPanelsDT3$nr_vars <- tgp_clinvar_unique$nr_vars[match(allPanelsDT3$"Gene Symbol",tgp_clinvar_unique$Symbol )]

allPanelsDT3$cum_gnomad_AF_NFE <- tgp_clinvar_unique$cum_gnomad_AF_NFE[match(allPanelsDT3$"Gene Symbol",tgp_clinvar_unique$Symbol )]
#allPanelsDT3$cum_gnomad_AC_NFE <- tgp_clinvar_unique$cum_gnomad_AC_NFE[match(allPanelsDT3$"Gene Symbol",tgp_clinvar_unique$Symbol )]
#allPanelsDT3$cum_gnomad_AN_NFE <- tgp_clinvar_unique$cum_gnomad_AN_NFE[match(allPanelsDT3$"Gene Symbol",tgp_clinvar_unique$Symbol )]
allPanelsDT3$gnomad_nr_vars_NFE <- tgp_clinvar_unique$gnomad_nr_vars_NFE[match(allPanelsDT3$"Gene Symbol",tgp_clinvar_unique$Symbol )]

allPanelsDT3[,cum_AF_sum:=sum(cum_AF, na.rm=T), by=Level4]
allPanelsDT3[,cum_AC_sum:=round(allPanelsDT3$cum_AF_sum*max(as.numeric(tgp_vep_nonsyn_rare$AN)))]
allPanelsDT3[,cum_AN_sum:=max(as.numeric(tgp_vep_nonsyn_rare$AN))]#round(max(as.numeric(tgp_vep_nonsyn_rare$AN))*allPanelsDT3$cum_AC_sum)]
allPanelsDT3[,nr_vars_sum:=sum(nr_vars, na.rm=T), by=Level4]

allPanelsDT3[,cum_gnomad_AF_NFE_sum:=sum(cum_gnomad_AF_NFE, na.rm=T), by=Level4]
allPanelsDT3[,cum_gnomad_AC_NFE_sum:=round(allPanelsDT3$cum_gnomad_AF_NFE_sum*nfe_cnt)]
allPanelsDT3[,cum_gnomad_AN_NFE_sum:=nfe_cnt]#round(nfe_cnt*allPanelsDT3$cum_gnomad_AC_NFE_sum)]
allPanelsDT3[,gnomad_nr_vars_NFE_sum:=sum(gnomad_nr_vars_NFE, na.rm=T), by=Level4]

allPanelsDT3_unique <- allPanelsDT3[-which(duplicated(Level4)),]
allPanelsDT3_unique$Symbol <- allPanelsDT3_unique$"Gene Symbol"
allPanelsDT3_unique$cum_AF_fold_change <- log2(allPanelsDT3_unique$cum_AF_sum/allPanelsDT3_unique$cum_gnomad_AF_NFE_sum)

suppressWarnings(pvals <- sapply(1:nrow(allPanelsDT3_unique), function(i){
       # print(i)
        row <- allPanelsDT3_unique[i,]
        if(row$cum_AC_sum ==0   ){return (NA)}
        M <- as.table(rbind(c(row$cum_AC_sum, row$cum_AN_sum), c(row$cum_gnomad_AC_NFE_sum, row$cum_gnomad_AN_NFE_sum)))
        dimnames(M) <- list(pop = c("PL", "NFE"),  value = c("AC", "AN"))
   # print(M)
        
      (Xsq <- chisq.test(M))  # Prints test summary
    Xsq$p.value
}))
allPanelsDT3_unique$pvals <- pvals*nrow(allPanelsDT3_unique) #  Bonferonni correction
```


```R
length(intersect(omim_recessive$symbol, (unique(allPanelsDT3$Gene))))
```


1823



```R

```


```R

```


```R
SuppTable_recessive_panels <- allPanelsDT3_unique[which(allPanelsDT3_unique$pvals < 0.01),]
```


```R
SuppTable_recessive_panels2 <- SuppTable_recessive_panels[,c("Level4", "cum_AF_sum", "cum_gnomad_AF_NFE_sum", "cum_AF_fold_change", "cum_AC_sum", "cum_AN_sum", "cum_gnomad_AC_NFE_sum", "cum_gnomad_AN_NFE_sum","nr_vars_sum" ,"gnomad_nr_vars_NFE_sum" , "pvals")]
```


```R
nrow(SuppTable_recessive_panels)
```


30



```R
fwrite (SuppTable_recessive_panels2, file="SuppTable_recessive_panels.xls", row.names=F)
```


```R
selPanels <- SuppTable_recessive_panels2#[which(abs(cum_AF_fold_change) > 1 & pvals < 0.01 & nr_vars_sum>=3 ),]
#selPanels[order(cum_AF_fold_change, decreasing=T), c("Level4", "cum_AF_sum", "cum_gnomad_AF_NFE_sum", "cum_AF_fold_change","nr_vars_sum" ,"gnomad_nr_vars_NFE_sum" , "pvals"),with=F]
```


```R
nrow(selPanels)
```


30



```R
set.seed(42)
p <- ggplot(selPanels, aes(cum_AF_sum, cum_gnomad_AF_NFE_sum)) +
    geom_point() +
    # geom_text(aes(label = Symbol),hjust = 0, nudge_x = 0.0001)+ 
    #geom_abline() + 
    #xlim(0,0.018) + ylim(0,0.01)+
    ggtitle("Cumulative allele frequencies of selected ClinVar variants in AR gene panels")+
    ylab("Per-panel cumulative allele frequencies in GnomAD NFE population")+
    xlab("Per-panel cumulative allele frequencies in Polish population")+
  scale_y_log10() + scale_x_log10()

p<- p +   geom_label_repel(
  aes(
    fill = "red", 
    label = Level4),
  fontface = 'bold.italic', 
    color = 'white',
  size = 3,
     segment.color = "black",
    #max.overlaps = 0,
    min.segment.length = 0,
  box.padding = unit(0.1, "lines"),
 point.padding = unit(0.1, "lines")
)
p +  theme(legend.position = "none") + theme(plot.title = element_text(size=11))
```

    Warning message:
    “Transformation introduced infinite values in continuous y-axis”


    
![png](output_55_1.png)
    


# OLD


```R
p1 <- ggplot(tcu, aes(cum_AF, cum_gnomad_AF)) +
    geom_point() +
     geom_text(aes(label = Symbol),hjust = 0, nudge_x = 0.00001)+ 
    geom_abline()#+
      #xlim(0, 0.0003) 
p2 <- ggplot(tcu, aes(cum_AF, cum_gnomad_AF_NFE)) +
    geom_point() +
    geom_text(aes(label = Symbol),hjust = 0, nudge_x = 0.00001) + 
    geom_abline()#+
    #xlim(0, 0.0003) 

p3 <- ggplot(tgp_cadd_unique, aes(cum_AF2, cum_gnomad_AF2)) +
    geom_point() +
     geom_text(aes(label = Symbol),hjust = 0, nudge_x = 0.0001)+ 
    geom_abline()+ xlim(0, 0.01) + ylim(0, 0.01) 

p4 <- ggplot(tgp_cadd_unique, aes(cum_AF2, cum_gnomad_AF2_NFE)) +
    geom_point() +
    geom_text(aes(label = Symbol),hjust = 0, nudge_x = 0.0001) + 
    geom_abline()+ xlim(0, 0.01) + ylim(0, 0.01) 

    
grid.arrange(p1, p2,p3,p4, nrow = 2)
```

    Warning message:
    “Removed 11 rows containing missing values (geom_point).”Warning message:
    “Removed 11 rows containing missing values (geom_text).”Warning message:
    “Removed 5 rows containing missing values (geom_point).”Warning message:
    “Removed 5 rows containing missing values (geom_text).”


    
![png](output_57_1.png)
    



```R
allPanelsDT3$clinvar_cum_AF2 <- tgp_clinvar_unique$cum_AF2[match(allPanelsDT3$"Gene Symbol",tgp_clinvar_unique$Symbol )]
allPanelsDT3$clinvar_cum_gnomad_AF2 <- tgp_clinvar_unique$cum_gnomad_AF2[match(allPanelsDT3$"Gene Symbol",tgp_clinvar_unique$Symbol )]
allPanelsDT3$clinvar_cum_gnomad_AF2_NFE <- tgp_clinvar_unique$cum_gnomad_AF2_NFE[match(allPanelsDT3$"Gene Symbol",tgp_clinvar_unique$Symbol )]

allPanelsDT3$cadd_cum_AF2 <- tgp_cadd_unique$cum_AF2[match(allPanelsDT3$"Gene Symbol",tgp_cadd_unique$Symbol )]
allPanelsDT3$cadd_cum_gnomad_AF2 <- tgp_cadd_unique$cum_gnomad_AF2[match(allPanelsDT3$"Gene Symbol",tgp_cadd_unique$Symbol )]
allPanelsDT3$cadd_cum_gnomad_AF2_NFE <- tgp_cadd_unique$cum_gnomad_AF2_NFE[match(allPanelsDT3$"Gene Symbol",tgp_cadd_unique$Symbol )]

allPanelsDT3[,clinvar_cum_AF2_sum:=sum(clinvar_cum_AF2, na.rm=T), by=Level4]
allPanelsDT3[,clinvar_cum_gnomad_AF2_sum:=sum(clinvar_cum_gnomad_AF2, na.rm=T), by=Level4]
allPanelsDT3[,clinvar_cum_gnomad_AF2_NFE_sum:=sum(clinvar_cum_gnomad_AF2_NFE, na.rm=T), by=Level4]

allPanelsDT3[,cadd_cum_AF2_sum:=sum(cadd_cum_AF2, na.rm=T), by=Level4]
allPanelsDT3[,cadd_cum_gnomad_AF2_sum:=sum(cadd_cum_gnomad_AF2, na.rm=T), by=Level4]
allPanelsDT3[,cadd_cum_gnomad_AF2_NFE_sum:=sum(cadd_cum_gnomad_AF2_NFE, na.rm=T), by=Level4]
allPanelsDT3_unique <- allPanelsDT3[-which(duplicated(Level4)),]
allPanelsDT3_unique$Symbol <- allPanelsDT3_unique$"Gene Symbol"
```


```R

```


```R
p1 <- ggplot(allPanelsDT3_unique, aes(clinvar_cum_AF2_sum, clinvar_cum_gnomad_AF2_sum)) +
    geom_point() +
     geom_text(aes(label = Level4),hjust = 0, nudge_x = 0.00001)+ 
    geom_abline()
      #xlim(0, 0.0003) 
p2 <- ggplot(allPanelsDT3_unique, aes(clinvar_cum_AF2_sum, clinvar_cum_gnomad_AF2_NFE_sum)) +
    geom_point() +
    geom_text(aes(label = Level4),hjust = 0, nudge_x = 0.00001) + 
    geom_abline()
    #xlim(0, 0.0003) 

p3 <- ggplot(allPanelsDT3_unique, aes(cadd_cum_AF2_sum, cadd_cum_gnomad_AF2_sum)) +
    geom_point() +
     geom_text(aes(label = Level4),hjust = 0, nudge_x = 0.0001)+ 
    geom_abline()+ xlim(0, 0.01) + ylim(0, 0.01) 

p4 <- ggplot(allPanelsDT3_unique, aes(cadd_cum_AF2_sum, cadd_cum_gnomad_AF2_NFE_sum)) +
    geom_point() +
    geom_text(aes(label = Level4),hjust = 0, nudge_x = 0.0001) + 
    geom_abline()+ xlim(0, 0.01) + ylim(0, 0.01) 

    
grid.arrange(p1, p2,p3,p4, nrow = 2)
```

    Warning message:
    “Removed 90 rows containing missing values (geom_point).”Warning message:
    “Removed 90 rows containing missing values (geom_text).”Warning message:
    “Removed 47 rows containing missing values (geom_point).”Warning message:
    “Removed 47 rows containing missing values (geom_text).”


    
![png](output_60_1.png)
    


### Estimated disease prelevance in gene panels (from geneApp)



```R
idx_list <- split(1:nrow(allPanelsDT3_unique), 1:10)

allPanelsDT3_unique <- allPanelsDT3_unique[order(clinvar_cum_gnomad_AF2_sum, decreasing=T)]
lapply(1:25, function(i){
    start <- ((i-1)*10+1)
    idx <- c(start:min(nrow(allPanelsDT3_unique),start+10))
    
melt_panel_clinvar <- melt(allPanelsDT3_unique[idx],  measure.vars = c("clinvar_cum_AF2_sum", "clinvar_cum_gnomad_AF2_sum", "clinvar_cum_gnomad_AF2_NFE_sum"), value.name = c("frequency"))
#melt_panel_cadd <- melt(allPanelsDT3_unique[idx],  measure.vars = c("cadd_cum_AF2_sum", "cadd_cum_gnomad_AF2_sum", "cadd_cum_gnomad_AF2_NFE_sum"), value.name = c("frequency"))

#library(repr)
#options(repr.plot.width=10, repr.plot.height=15)

p1 <- ggplot(data=melt_panel_clinvar, aes(x=Level4, y=frequency, fill=variable)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))

p1
})

#p2 <-     ggplot(data=melt_panel_cadd, aes(x=Level4, y=frequency, fill=variable)) +
#geom_bar(stat="identity", position=position_dodge()) +
#theme(axis.text.x = element_text(angle = 45, hjust=1))


```

    Warning message in split.default(1:nrow(allPanelsDT3_unique), 1:10):
    “data length is not a multiple of split variable”


    
![png](output_62_1.png)
    



    
![png](output_62_2.png)
    



    
![png](output_62_3.png)
    



    
![png](output_62_4.png)
    



    
![png](output_62_5.png)
    



    
![png](output_62_6.png)
    



    
![png](output_62_7.png)
    



    
![png](output_62_8.png)
    



    
![png](output_62_9.png)
    



    
![png](output_62_10.png)
    



    
![png](output_62_11.png)
    



    
![png](output_62_12.png)
    



    
![png](output_62_13.png)
    



    
![png](output_62_14.png)
    



    
![png](output_62_15.png)
    



    
![png](output_62_16.png)
    



    
![png](output_62_17.png)
    



    
![png](output_62_18.png)
    



    
![png](output_62_19.png)
    



    
![png](output_62_20.png)
    



    
![png](output_62_21.png)
    



    
![png](output_62_22.png)
    



    
![png](output_62_23.png)
    



    [[1]]
    
    [[2]]
    
    [[3]]
    
    [[4]]
    
    [[5]]
    
    [[6]]
    
    [[7]]
    
    [[8]]
    
    [[9]]
    
    [[10]]
    
    [[11]]
    
    [[12]]
    
    [[13]]
    
    [[14]]
    
    [[15]]
    
    [[16]]
    
    [[17]]
    
    [[18]]
    
    [[19]]
    
    [[20]]
    
    [[21]]
    
    [[22]]
    
    [[23]]
    
    [[24]]
    
    [[25]]




    
![png](output_62_25.png)
    



    
![png](output_62_26.png)
    



```R
idx_list <- split(1:nrow(allPanelsDT3_unique), 1:10)

allPanelsDT3_unique <- allPanelsDT3_unique[order(cadd_cum_gnomad_AF2_sum, decreasing=T)]
lapply(1:25, function(i){
    start <- ((i-1)*10+1)
    idx <- c(start:min(nrow(allPanelsDT3_unique),start+10))
    
#melt_panel_clinvar <- melt(allPanelsDT3_unique[idx],  measure.vars = c("clinvar_cum_AF2_sum", "clinvar_cum_gnomad_AF2_sum", "clinvar_cum_gnomad_AF2_NFE_sum"), value.name = c("frequency"))
melt_panel_cadd <- melt(allPanelsDT3_unique[idx],  measure.vars = c("cadd_cum_AF2_sum", "cadd_cum_gnomad_AF2_sum", "cadd_cum_gnomad_AF2_NFE_sum"), value.name = c("frequency"))

#library(repr)
#options(repr.plot.width=10, repr.plot.height=15)

p1 <- ggplot(data=melt_panel_cadd, aes(x=Level4, y=frequency, fill=variable)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))

p1
})



```

    Warning message in split.default(1:nrow(allPanelsDT3_unique), 1:10):
    “data length is not a multiple of split variable”


    
![png](output_63_1.png)
    



    
![png](output_63_2.png)
    



    
![png](output_63_3.png)
    



    
![png](output_63_4.png)
    



    
![png](output_63_5.png)
    



    
![png](output_63_6.png)
    



    
![png](output_63_7.png)
    



    
![png](output_63_8.png)
    



    
![png](output_63_9.png)
    



    
![png](output_63_10.png)
    



    
![png](output_63_11.png)
    



    
![png](output_63_12.png)
    



    
![png](output_63_13.png)
    



    
![png](output_63_14.png)
    



    
![png](output_63_15.png)
    



    
![png](output_63_16.png)
    



    
![png](output_63_17.png)
    



    
![png](output_63_18.png)
    



    
![png](output_63_19.png)
    



    
![png](output_63_20.png)
    



    
![png](output_63_21.png)
    



    
![png](output_63_22.png)
    



    
![png](output_63_23.png)
    



    [[1]]
    
    [[2]]
    
    [[3]]
    
    [[4]]
    
    [[5]]
    
    [[6]]
    
    [[7]]
    
    [[8]]
    
    [[9]]
    
    [[10]]
    
    [[11]]
    
    [[12]]
    
    [[13]]
    
    [[14]]
    
    [[15]]
    
    [[16]]
    
    [[17]]
    
    [[18]]
    
    [[19]]
    
    [[20]]
    
    [[21]]
    
    [[22]]
    
    [[23]]
    
    [[24]]
    
    [[25]]




    
![png](output_63_25.png)
    



    
![png](output_63_26.png)
    


# OLD analysis:


```R
# Conseuqneces summary
csq_table <- as.data.table(table(tgp_vep$Consequence))
head(csq_table)
```


    Error in table(tgp_vep$Consequence): object 'tgp_vep' not found
    Traceback:


    1. as.data.table(table(tgp_vep$Consequence))

    2. table(tgp_vep$Consequence)



```R
# Select "main" consequence per each row 
csq_table$csq_first <- sapply(strsplit(csq_table$V1, ","), `[`, 1)
csq_first <- as.data.table(csq_table %>% mutate(csq_first = toupper(csq_first)) %>% 
    group_by(csq_first) %>% summarise(N = sum(N)))
dim(csq_first)
head(csq_first)
```


<ol class=list-inline>
	<li>24</li>
	<li>2</li>
</ol>




<table>
<thead><tr><th scope=col>csq_first</th><th scope=col>N</th></tr></thead>
<tbody>
	<tr><td>3_PRIME_UTR_VARIANT              </td><td> 341486                          </td></tr>
	<tr><td>5_PRIME_UTR_VARIANT              </td><td>  57721                          </td></tr>
	<tr><td>CODING_SEQUENCE_VARIANT          </td><td>      9                          </td></tr>
	<tr><td>DOWNSTREAM_GENE_VARIANT          </td><td>2131350                          </td></tr>
	<tr><td>FRAMESHIFT_VARIANT               </td><td>   4416                          </td></tr>
	<tr><td>INCOMPLETE_TERMINAL_CODON_VARIANT</td><td>      4                          </td></tr>
</tbody>
</table>




```R
# Plot consequences
p1 <- ggplot(csq_first, aes(reorder(csq_first, -N), N)) + 
    geom_col() +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    scale_y_log10() +
    labs(x = "Consequence", y = "Frequency (log10)")
p2 <- ggplot(csq_first, aes(reorder(csq_first, -N), N)) + 
    geom_col() +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(x = "Consequence", y = "Frequency")
    #theme(plot.title = "Liczba wariantów dla poszczególnych głownych konsekwencji", axis.title.x = element_blank(), axis.title.y = element_blank())


grid.arrange(p2, p1, nrow = 1, top = textGrob("Total no. of variants obs. by VEP consequence",gp=gpar(fontsize=20)))
```


    
![png](output_67_0.png)
    



```R

```


```R
# Recessive clinvar pathogenic/likely pathogenic variants
clinvar_recessive <- tgp_vep[which(Symbol %in% omim_recessive$symbol & CLNSIG %in% c("Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic)"))]
```


```R
dim(tgp_vep[which(AC==0)])
```


<ol class=list-inline>
	<li>429163</li>
	<li>154</li>
</ol>




```R
# 10 variants with AC=0
dim(clinvar_recessive)
dim(clinvar_recessive[which(AC>0)])
head(clinvar_recessive[which(AC==0)])
```


```R
table(clinvar_recessive[,Consequence])
```


    
                                            5_prime_UTR_variant 
                                                              3 
                                        downstream_gene_variant 
                                                             12 
                                             frameshift_variant 
                                                            118 
                       frameshift_variant,splice_region_variant 
                                                              3 
                                               inframe_deletion 
                                                              6 
                                              inframe_insertion 
                                                              1 
                                                 intron_variant 
                                                             18 
                                               missense_variant 
                                                            251 
                         missense_variant,splice_region_variant 
                                                              8 
                                        splice_acceptor_variant 
                                                             19 
                                           splice_donor_variant 
                                                             40 
    splice_donor_variant,coding_sequence_variant,intron_variant 
                                                              2 
                           splice_region_variant,intron_variant 
                                                             17 
                       splice_region_variant,synonymous_variant 
                                                              3 
                                                     start_lost 
                                                              6 
                                                    stop_gained 
                                                            121 
                                 stop_gained,frameshift_variant 
                                                              2 
                              stop_gained,splice_region_variant 
                                                              2 
                                                      stop_lost 
                                                              2 
                                             synonymous_variant 
                                                              3 
                                          upstream_gene_variant 
                                                             11 



```R
clinvar_recessive[which(Consequence == "upstream_gene_variant"),ClinVar_url]
```


<ol class=list-inline>
	<li>'https://www.ncbi.nlm.nih.gov/clinvar/variation/12076/'</li>
	<li>'https://www.ncbi.nlm.nih.gov/clinvar/variation/265135/'</li>
	<li>'https://www.ncbi.nlm.nih.gov/clinvar/variation/626999/'</li>
	<li>'https://www.ncbi.nlm.nih.gov/clinvar/variation/116/'</li>
	<li>'https://www.ncbi.nlm.nih.gov/clinvar/variation/188834/'</li>
	<li>'https://www.ncbi.nlm.nih.gov/clinvar/variation/573257/'</li>
	<li>'https://www.ncbi.nlm.nih.gov/clinvar/variation/162203/'</li>
	<li>'https://www.ncbi.nlm.nih.gov/clinvar/variation/223044/'</li>
	<li>'https://www.ncbi.nlm.nih.gov/clinvar/variation/16663/'</li>
	<li>'https://www.ncbi.nlm.nih.gov/clinvar/variation/223019/'</li>
	<li>'https://www.ncbi.nlm.nih.gov/clinvar/variation/31077/'</li>
</ol>




```R
dim(clinvar_recessive[which(AF > gnomAD_AF_NFE), c("key","gnomAD_AF","gnomAD_AF_NFE", "AF")])
head(clinvar_recessive[which(AF > gnomAD_AF_NFE), c("key","gnomAD_AF","gnomAD_AF_NFE", "AF")])
```


```R
clinvar_gene_table <- as.data.table(sort(table(clinvar_recessive$Symbol), decreasing = TRUE))
```


```R
# Genes with the highest number of Clinvar pathogenic varianst
clinvar_gene_table_2 <- clinvar_gene_table[which(N > 2)]
ggplot(clinvar_gene_table[which(N > 2)], aes(reorder(V1, -N), N)) + 
    geom_col() +
    scale_x_discrete(guide = guide_axis(angle = 60)) +
    labs(title = "Genes with most P/LP ClinVar variants in TGP cohort", x = "Symbol", y = "Count")
```


    
![png](output_76_0.png)
    



```R
clinvar_recessive$AF <- as.numeric(clinvar_recessive$AF)
```


```R
melt_clinvar <- melt(clinvar_recessive, id.vars = c("key2"), measure.vars = c("gnomAD_AF", "gnomAD_AF_NFE", "AF"), variable.name = c("database"), value.name = c("frequency"))
```


```R
dim(melt_clinvar)
#head(melt_clinvar)
```


<ol class=list-inline>
	<li>1944</li>
	<li>3</li>
</ol>




```R
ggplot(melt_clinvar, aes(database,frequency, fill=frequency)) +
    geom_boxplot() +
    geom_line(aes(group=key2)) +
    geom_point(aes(fill=frequency,group=database),size=2,shape=21, position = position_dodge(0.2)) +
    theme(legend.position = "none") +
    scale_y_log10() +
    labs(title = "ClinVar P/LP variant freq. across gnomAD database and TGP cohort", x = "Population", y = "Frequency (log10)")
```

    Warning message:
    “Transformation introduced infinite values in continuous y-axis”Warning message:
    “Transformation introduced infinite values in continuous y-axis”Warning message:
    “Transformation introduced infinite values in continuous y-axis”Warning message:
    “Removed 211 rows containing non-finite values (stat_boxplot).”


    
![png](output_80_1.png)
    



```R
clinvar_recessive_overrep_TGP <- clinvar_recessive[which(AF > gnomAD_AF_NFE & AF > gnomAD_AF), c("key2","gnomAD_AF","gnomAD_AF_NFE", "AF")]
#head(clinvar_recessive_overrep_TGP)
melt_clinvar_overrep_TGP <- melt(clinvar_recessive_overrep_TGP, id.vars = c("key2"), measure.vars = c("gnomAD_AF", "gnomAD_AF_NFE", "AF"), variable.name = c("database"), value.name = c("frequency"))
#head(melt_clinvar_overrep_TGP)
```


```R
clinvar_recessive_diff <- clinvar_recessive[, diff_AF := gnomAD_AF - AF][, diff_AF_NFE := gnomAD_AF_NFE - AF]
```


```R
#head(clinvar_recessive_diff[which(diff_AF < 0 & diff_AF_NFE < 0), c("key","gnomAD_AF","gnomAD_AF_NFE", "AF", "diff_AF", "diff_AF_NFE")])
dim(clinvar_recessive_diff[which(diff_AF < 0 & diff_AF_NFE < 0)])
```


<ol class=list-inline>
	<li>551</li>
	<li>156</li>
</ol>




```R
#ggplot(clinvar_recessive_overrep_TGP, aes(database, frequency, fill = frequency)) +
#    geom_boxplot() +
#    geom_line(aes(group = key2)) +
#    geom_point(aes(fill = frequency, group = database), size = 2, shape = 21, position = position_dodge(0.2)) +
#    theme(legend.position = "none") +
#    scale_y_log10() +
#    labs(title = "ClinVar P/LP variant freq. across gnomAD database and TGP cohort", x = "Population", y = "Frequency (log10)")

#ggplot(clinvar_recessive_overrep_TGP, aes(x=database, y=frequency)) +
#    geom_boxplot() +
#    geom_jitter(colour="lightblue", height = 0) +
#    scale_y_log10() +
#    labs(title = "ClinVar P/LP variant freq. across gnomAD database and TGP cohort", x = "Population", y = "Frequency (log10)")
```


```R
ggplot(melt_clinvar, aes(x=database, y=frequency)) +
    geom_boxplot() +
    geom_jitter(colour="lightblue", height = 0) +
    scale_y_log10() +
    labs(title = "ClinVar P/LP variant freq. across gnomAD database and TGP cohort", x = "Population", y = "Frequency (log10)")
```

    Warning message:
    “Transformation introduced infinite values in continuous y-axis”Warning message:
    “Transformation introduced infinite values in continuous y-axis”Warning message:
    “Removed 211 rows containing non-finite values (stat_boxplot).”Warning message:
    “Removed 211 rows containing missing values (geom_point).”


    
![png](output_85_1.png)
    



```R
ggplot(melt_clinvar, aes(x=database, y=frequency)) +
    geom_violin() +
    geom_count(alpha=0.5) +
    scale_y_log10() +
    labs(title = "ClinVar P/LP variant freq. across gnomAD database and TGP cohort", x = "Population", y = "Frequency (log10)")
```

    Warning message:
    “Transformation introduced infinite values in continuous y-axis”Warning message:
    “Transformation introduced infinite values in continuous y-axis”Warning message:
    “Removed 211 rows containing non-finite values (stat_ydensity).”Warning message:
    “Removed 211 rows containing non-finite values (stat_sum).”


    
![png](output_86_1.png)
    



```R

```


```R
tgp_vep_gene_table <- as.data.table(sort(table(tgp_vep_nonsyn_rare$Symbol), decreasing = TRUE))
ggplot(tgp_vep_gene_table[2:31], aes(reorder(V1, -N), N)) + 
    geom_col() +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    labs(title = "Total no. of variants per gene in TGP cohort", x = "Symbol", y = "Count")
```


    
![png](output_88_0.png)
    



```R
tgp_vep_gene_table_singleton <- as.data.table(sort(table(tgp_vep_nonsyn_rare[which(AF < 0.01 & AC == 1)]$Symbol), decreasing = TRUE))
```


```R
ggplot(tgp_vep_gene_table_singleton[2:31], aes(reorder(V1, -N), N)) + 
    geom_col() +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    labs(title = "Total no. of singleton variants per gene in TGP cohort", x = "Symbol", y = "Count")
```


    
![png](output_90_0.png)
    



```R
tgp_clinvar[which(Symbol == "OTOG" & CLNSIG %in% c("Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic")), c("AAchange", "AF", "AC", "gnomAD_AF", "gnomAD_AF_NFE")]
```


<table>
<thead><tr><th scope=col>AAchange</th><th scope=col>AF</th><th scope=col>AC</th><th scope=col>gnomAD_AF</th><th scope=col>gnomAD_AF_NFE</th></tr></thead>
<tbody>
	<tr><td>p.Gln834Ter        </td><td>0.00742312         </td><td>14                 </td><td>0.0003906800       </td><td>0.0006813250       </td></tr>
	<tr><td>p.Ala1838ProfsTer31</td><td>0.000530223        </td><td>1                  </td><td>0.0000349055       </td><td>0.0000774593       </td></tr>
</tbody>
</table>




```R

```
