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

    
    Attaching package: ‘dplyr’
    
    The following object is masked from ‘package:glue’:
    
        collapse
    
    The following objects are masked from ‘package:data.table’:
    
        between, first, last
    
    The following objects are masked from ‘package:stats’:
    
        filter, lag
    
    The following objects are masked from ‘package:base’:
    
        intersect, setdiff, setequal, union
    
    Loading required package: ggplot2
    
    Attaching package: ‘plotly’
    
    The following object is masked from ‘package:ggplot2’:
    
        last_plot
    
    The following object is masked from ‘package:stats’:
    
        filter
    
    The following object is masked from ‘package:graphics’:
    
        layout
    
    
    Attaching package: ‘gridExtra’
    
    The following object is masked from ‘package:dplyr’:
    
        combine
    


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

## Start from here


```R
#load("~/tgp_vep_nonsyn_rare_2021_06_27.RData")
load("~/gnomad_nonsyn_rare_2021_06_27.RData" )
load("~/tgp_recessive_2021_06_27.RData" )
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
tgp_clinvar_unique$cum_AN <- round(max(as.numeric(tgp_vep_nonsyn_rare$AN))*tgp_clinvar_unique$cum_AC)
tgp_clinvar_unique$cum_gnomad_AF_NFE <- gnomad_clinvar_unique$cum_gnomAD_AF_NFE[match(tgp_clinvar_unique$Symbol,gnomad_clinvar_unique$SYMBOL )]
tgp_clinvar_unique$gnomad_nr_vars_NFE <- gnomad_clinvar_unique$gnomAD_nr_vars_NFE[match(tgp_clinvar_unique$Symbol,gnomad_clinvar_unique$SYMBOL )]
nfe_cnt <- 32299
tgp_clinvar_unique$cum_gnomad_AC_NFE <- round(tgp_clinvar_unique$cum_gnomad_AF_NFE*nfe_cnt)
tgp_clinvar_unique$cum_gnomad_AN_NFE <- round(nfe_cnt*tgp_clinvar_unique$cum_gnomad_AC_NFE)
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
length(which(tgp_clinvar_unique$pvals < 0.01))
```


152



```R
tgp_clinvar_unique[which(tgp_clinvar_unique$pvals < 0.01),]
```


<table>
<thead><tr><th scope=col>Symbol</th><th scope=col>cum_AF</th><th scope=col>nr_vars</th><th scope=col>cum_AC</th><th scope=col>cum_AN</th><th scope=col>cum_gnomad_AF_NFE</th><th scope=col>gnomad_nr_vars_NFE</th><th scope=col>cum_gnomad_AC_NFE</th><th scope=col>cum_gnomad_AN_NFE</th><th scope=col>cum_AF_fold_change</th><th scope=col>pvals</th></tr></thead>
<tbody>
	<tr><td>C2                                                                                                                             </td><td>0.023329800                                                                                                                    </td><td> 1                                                                                                                             </td><td>44                                                                                                                             </td><td>82984                                                                                                                          </td><td>0.0075581600                                                                                                                   </td><td>  1                                                                                                                            </td><td>244                                                                                                                            </td><td> 7880956                                                                                                                       </td><td> 1.62606697                                                                                                                    </td><td>0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000002141842</td></tr>
	<tr><td>TGM5                                                                                                                           </td><td>0.016436923                                                                                                                    </td><td> 2                                                                                                                             </td><td>31                                                                                                                             </td><td>58466                                                                                                                          </td><td>0.0040888476                                                                                                                   </td><td>  5                                                                                                                            </td><td>132                                                                                                                            </td><td> 4263468                                                                                                                       </td><td> 2.00717405                                                                                                                    </td><td>0.00000000000000000000000000000000000000000000000000000000000000000000000000000025439502116676947905558192120088882771738319583</td></tr>
	<tr><td>MPO                                                                                                                            </td><td>0.015379840                                                                                                                    </td><td> 4                                                                                                                             </td><td>29                                                                                                                             </td><td>54694                                                                                                                          </td><td>0.0129199716                                                                                                                   </td><td>  8                                                                                                                            </td><td>417                                                                                                                            </td><td>13468683                                                                                                                       </td><td> 0.25143760                                                                                                                    </td><td>0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000014890020510988255311577648500487216368704</td></tr>
	<tr><td>VWF                                                                                                                            </td><td>0.014316599                                                                                                                    </td><td> 4                                                                                                                             </td><td>27                                                                                                                             </td><td>50922                                                                                                                          </td><td>0.0076959417                                                                                                                   </td><td> 31                                                                                                                            </td><td>249                                                                                                                            </td><td> 8042451                                                                                                                       </td><td> 0.89551904                                                                                                                    </td><td>0.00000000000000000000000000000000000000000000000000000000000000000000000000014733761846323446039134618701815961688642572228603</td></tr>
	<tr><td>DHCR7                                                                                                                          </td><td>0.014316013                                                                                                                    </td><td> 3                                                                                                                             </td><td>27                                                                                                                             </td><td>50922                                                                                                                          </td><td>0.0101287216                                                                                                                   </td><td> 43                                                                                                                            </td><td>327                                                                                                                            </td><td>10561773                                                                                                                       </td><td> 0.49917766                                                                                                                    </td><td>0.00000000000000000000000000000000000000000000000000000000000000000000000000000206895866775148326815995266795601855217938897207</td></tr>
	<tr><td>PAH                                                                                                                            </td><td>0.012725350                                                                                                                    </td><td>11                                                                                                                             </td><td>24                                                                                                                             </td><td>45264                                                                                                                          </td><td>0.0114455812                                                                                                                   </td><td>106                                                                                                                            </td><td>370                                                                                                                            </td><td>11950630                                                                                                                       </td><td> 0.15291461                                                                                                                    </td><td>0.00000000000000000000000000000000000000000000000000000000000000000000020065391942528339883036651236457771239559768707859813592</td></tr>
	<tr><td>CFTR                                                                                                                           </td><td>0.012195123                                                                                                                    </td><td> 3                                                                                                                             </td><td>23                                                                                                                             </td><td>43378                                                                                                                          </td><td>0.0213872883                                                                                                                   </td><td> 70                                                                                                                            </td><td>691                                                                                                                            </td><td>22318609                                                                                                                       </td><td>-0.81044927                                                                                                                    </td><td>0.00000000000000000000000000000000000000000000000000000000000000000000207439595759666801097498432177758635577968584804924075649</td></tr>
	<tr><td>CYP21A2                                                                                                                        </td><td>0.011961077                                                                                                                    </td><td> 4                                                                                                                             </td><td>23                                                                                                                             </td><td>43378                                                                                                                          </td><td>0.0132150295                                                                                                                   </td><td> 11                                                                                                                            </td><td>427                                                                                                                            </td><td>13791673                                                                                                                       </td><td>-0.14383235                                                                                                                    </td><td>0.00000000000000000000000000000000000000000000000000000000000000000004919313776847481362557850066538303273259798526216536574726</td></tr>
	<tr><td>ALDOB                                                                                                                          </td><td>0.009544009                                                                                                                    </td><td> 5                                                                                                                             </td><td>18                                                                                                                             </td><td>33948                                                                                                                          </td><td>0.0063657069                                                                                                                   </td><td> 15                                                                                                                            </td><td>206                                                                                                                            </td><td> 6653594                                                                                                                       </td><td> 0.58427467                                                                                                                    </td><td>0.00000000000000000000000000000000000000000000000001043208886356037628120153618215290329491972353139079259108033322571200206567</td></tr>
	<tr><td>C8B                                                                                                                            </td><td>0.008483563                                                                                                                    </td><td> 2                                                                                                                             </td><td>16                                                                                                                             </td><td>30176                                                                                                                          </td><td>0.0024491506                                                                                                                   </td><td>  3                                                                                                                            </td><td> 79                                                                                                                            </td><td> 2551621                                                                                                                       </td><td> 1.79238882                                                                                                                    </td><td>0.00000000000000000000000000000000000000030190970732010409191878782323982133265267159696029473575387564200096464798500505878468</td></tr>
	<tr><td>PADI3                                                                                                                          </td><td>0.008483560                                                                                                                    </td><td> 2                                                                                                                             </td><td>16                                                                                                                             </td><td>30176                                                                                                                          </td><td>0.0183699566                                                                                                                   </td><td>  7                                                                                                                            </td><td>593                                                                                                                            </td><td>19153307                                                                                                                       </td><td>-1.11460652                                                                                                                    </td><td>0.00000000000000000000000000000000000000000000002582872706667123548112641636144582442489458083411487531437091714695152684180788</td></tr>
	<tr><td>CBS                                                                                                                            </td><td>0.008055240                                                                                                                    </td><td> 1                                                                                                                             </td><td>15                                                                                                                             </td><td>28290                                                                                                                          </td><td>0.0002146807                                                                                                                   </td><td>  5                                                                                                                            </td><td>  7                                                                                                                            </td><td>  226093                                                                                                                       </td><td> 5.22966317                                                                                                                    </td><td>0.00000000000014181807686266376649545292366082095133561838395941379076248267665505409240722656250000000000000000000000000000000</td></tr>
	<tr><td>ABCA4                                                                                                                          </td><td>0.007953355                                                                                                                    </td><td> 9                                                                                                                             </td><td>15                                                                                                                             </td><td>28290                                                                                                                          </td><td>0.0075767334                                                                                                                   </td><td>103                                                                                                                            </td><td>245                                                                                                                            </td><td> 7913255                                                                                                                       </td><td> 0.06998758                                                                                                                    </td><td>0.00000000000000000000000000000000000000000123136858575660530611047043797798215236602678248066283551896406385490789330074044508</td></tr>
	<tr><td>OTOG                                                                                                                           </td><td>0.007953343                                                                                                                    </td><td> 2                                                                                                                             </td><td>15                                                                                                                             </td><td>28290                                                                                                                          </td><td>0.0015798706                                                                                                                   </td><td> 14                                                                                                                            </td><td> 51                                                                                                                            </td><td> 1647249                                                                                                                       </td><td> 2.33175499                                                                                                                    </td><td>0.00000000000000000000000000000000009247482169364185709852098431280980573442425018862891374099274291437656642734262509079846623</td></tr>
	<tr><td>GP6                                                                                                                            </td><td>0.007953340                                                                                                                    </td><td> 1                                                                                                                             </td><td>15                                                                                                                             </td><td>28290                                                                                                                          </td><td>0.0019975296                                                                                                                   </td><td>  4                                                                                                                            </td><td> 65                                                                                                                            </td><td> 2099435                                                                                                                       </td><td> 1.99334396                                                                                                                    </td><td>0.00000000000000000000000000000000000131256023911555491314883250268676812569468816093400098975935640244728170406428481682633711</td></tr>
	<tr><td>COL18A1                                                                                                                        </td><td>0.006892893                                                                                                                    </td><td> 2                                                                                                                             </td><td>13                                                                                                                             </td><td>24518                                                                                                                          </td><td>0.0026180830                                                                                                                   </td><td> 10                                                                                                                            </td><td> 85                                                                                                                            </td><td> 2745415                                                                                                                       </td><td> 1.39659879                                                                                                                    </td><td>0.00000000000000000000000000000000205124477400919983939114654966761739373117117585197873557531819653040903538326128341807812194</td></tr>
	<tr><td>TYR                                                                                                                            </td><td>0.006362680                                                                                                                    </td><td> 4                                                                                                                             </td><td>12                                                                                                                             </td><td>22632                                                                                                                          </td><td>0.0047262731                                                                                                                   </td><td> 36                                                                                                                            </td><td>153                                                                                                                            </td><td> 4941747                                                                                                                       </td><td> 0.42893157                                                                                                                    </td><td>0.00000000000000000000000000000000984048456587864289926377713534198373474963716092286870553970670896986451925582197754714675053</td></tr>
	<tr><td>MASP1                                                                                                                          </td><td>0.006362670                                                                                                                    </td><td> 1                                                                                                                             </td><td>12                                                                                                                             </td><td>22632                                                                                                                          </td><td>0.0005730213                                                                                                                   </td><td>  3                                                                                                                            </td><td> 19                                                                                                                            </td><td>  613681                                                                                                                       </td><td> 3.47297163                                                                                                                    </td><td>0.00000000000000000000315569309601283555491664750633367127293233956292679768775499490884861586437182268127799034118652343750000</td></tr>
	<tr><td>ATP7B                                                                                                                          </td><td>0.006362669                                                                                                                    </td><td> 4                                                                                                                             </td><td>12                                                                                                                             </td><td>22632                                                                                                                          </td><td>0.0047390235                                                                                                                   </td><td> 63                                                                                                                            </td><td>153                                                                                                                            </td><td> 4941747                                                                                                                       </td><td> 0.42504226                                                                                                                    </td><td>0.00000000000000000000000000000000984048456587864289926377713534198373474963716092286870553970670896986451925582197754714675053</td></tr>
	<tr><td>SLC26A2                                                                                                                        </td><td>0.005832456                                                                                                                    </td><td> 4                                                                                                                             </td><td>11                                                                                                                             </td><td>20746                                                                                                                          </td><td>0.0025259415                                                                                                                   </td><td>  7                                                                                                                            </td><td> 82                                                                                                                            </td><td> 2648518                                                                                                                       </td><td> 1.20728229                                                                                                                    </td><td>0.00000000000000000000000000037844978621945383950495161487382262103939493342117593229437586039478531884684697761400684612453915</td></tr>
	<tr><td>ABCC6                                                                                                                          </td><td>0.005302232                                                                                                                    </td><td> 5                                                                                                                             </td><td>10                                                                                                                             </td><td>18860                                                                                                                          </td><td>0.0050648990                                                                                                                   </td><td> 61                                                                                                                            </td><td>164                                                                                                                            </td><td> 5297036                                                                                                                       </td><td> 0.06606630                                                                                                                    </td><td>0.00000000000000000000000000474105667706053154519429738613331919473605907939852579623179709473226569119153062104032869683578610</td></tr>
	<tr><td>C19orf12                                                                                                                       </td><td>0.004783850                                                                                                                    </td><td> 2                                                                                                                             </td><td> 9                                                                                                                             </td><td>16974                                                                                                                          </td><td>0.0002478660                                                                                                                   </td><td>  4                                                                                                                            </td><td>  8                                                                                                                            </td><td>  258392                                                                                                                       </td><td> 4.27053986                                                                                                                    </td><td>0.00000000002695514571709420661634299426381682074360668366352911107242107391357421875000000000000000000000000000000000000000000</td></tr>
	<tr><td>MPZL2                                                                                                                          </td><td>0.004772010                                                                                                                    </td><td> 2                                                                                                                             </td><td> 9                                                                                                                             </td><td>16974                                                                                                                          </td><td>0.0011305156                                                                                                                   </td><td>  2                                                                                                                            </td><td> 37                                                                                                                            </td><td> 1195063                                                                                                                       </td><td> 2.07761617                                                                                                                    </td><td>0.00000000000000000003067594699341852056652627617517025448219165511425289170552299644612048723502084612846374511718750000000000</td></tr>
	<tr><td>SORD                                                                                                                           </td><td>0.004772010                                                                                                                    </td><td> 2                                                                                                                             </td><td> 9                                                                                                                             </td><td>16974                                                                                                                          </td><td>0.0005419470                                                                                                                   </td><td>  1                                                                                                                            </td><td> 18                                                                                                                            </td><td>  581382                                                                                                                       </td><td> 3.13837339                                                                                                                    </td><td>0.00000000000000014765513883215730425751710933861530779746113178213517969084023206960409879684448242187500000000000000000000000</td></tr>
	<tr><td>SLC26A4                                                                                                                        </td><td>0.004772009                                                                                                                    </td><td> 5                                                                                                                             </td><td> 9                                                                                                                             </td><td>16974                                                                                                                          </td><td>0.0048638286                                                                                                                   </td><td> 43                                                                                                                            </td><td>157                                                                                                                            </td><td> 5070943                                                                                                                       </td><td>-0.02749563                                                                                                                    </td><td>0.00000000000000000000000513427017453163915776052897965676412754645787061744843735722329189363755119757115608081221580505371094</td></tr>
	<tr><td>ALOX12B                                                                                                                        </td><td>0.004772000                                                                                                                    </td><td> 2                                                                                                                             </td><td> 9                                                                                                                             </td><td>16974                                                                                                                          </td><td>0.0009913307                                                                                                                   </td><td> 30                                                                                                                            </td><td> 32                                                                                                                            </td><td> 1033568                                                                                                                       </td><td> 2.26715573                                                                                                                    </td><td>0.00000000000000000013318365755227897165665039680807174055253251147175455387436371346154828643193468451499938964843750000000000</td></tr>
	<tr><td>M1AP                                                                                                                           </td><td>0.004772000                                                                                                                    </td><td> 1                                                                                                                             </td><td> 9                                                                                                                             </td><td>16974                                                                                                                          </td><td>0.0039952900                                                                                                                   </td><td>  1                                                                                                                            </td><td>129                                                                                                                            </td><td> 4166571                                                                                                                       </td><td> 0.25629382                                                                                                                    </td><td>0.00000000000000000000001009256305322406859412926978281310262389736037388445185823740424518185632152267316996585577726364135742</td></tr>
	<tr><td>FKBP14                                                                                                                         </td><td>0.004772000                                                                                                                    </td><td> 1                                                                                                                             </td><td> 9                                                                                                                             </td><td>16974                                                                                                                          </td><td>0.0008674370                                                                                                                   </td><td>  1                                                                                                                            </td><td> 28                                                                                                                            </td><td>  904372                                                                                                                       </td><td> 2.45976316                                                                                                                    </td><td>0.00000000000000000057392277955507753845204492677021954366036928513517315596226264773349612369202077388763427734375000000000000</td></tr>
	<tr><td>GNRHR                                                                                                                          </td><td>0.004241790                                                                                                                    </td><td> 2                                                                                                                             </td><td> 8                                                                                                                             </td><td>15088                                                                                                                          </td><td>0.0058397254                                                                                                                   </td><td>  8                                                                                                                            </td><td>189                                                                                                                            </td><td> 6104511                                                                                                                       </td><td>-0.46122733                                                                                                                    </td><td>0.00000000000000000000335029747022857191378076861023988466049833777696650834269096255746411827658448601141571998596191406250000</td></tr>
	<tr><td>PIGV                                                                                                                           </td><td>0.004241786                                                                                                                    </td><td> 3                                                                                                                             </td><td> 8                                                                                                                             </td><td>15088                                                                                                                          </td><td>0.0004646057                                                                                                                   </td><td>  6                                                                                                                            </td><td> 15                                                                                                                            </td><td>  484485                                                                                                                       </td><td> 3.19059308                                                                                                                    </td><td>0.00000000000005320855579008787686253948991302215313270364965836378701169451232999563217163085937500000000000000000000000000000</td></tr>
	<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><td>CABP2        </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0010532555 </td><td> 2           </td><td> 34          </td><td>1098166      </td><td> 0.59477912  </td><td>0.00004917893</td></tr>
	<tr><td>RNASEH2B     </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0023698735 </td><td> 9           </td><td> 77          </td><td>2487023      </td><td>-0.57517549  </td><td>0.00002415008</td></tr>
	<tr><td>EIF2B2       </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0007435570 </td><td> 5           </td><td> 24          </td><td> 775176      </td><td> 1.09711932  </td><td>0.00008032737</td></tr>
	<tr><td>IVD          </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0017653392 </td><td>13           </td><td> 57          </td><td>1841043      </td><td>-0.15031085  </td><td>0.00002961165</td></tr>
	<tr><td>B9D1         </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0002477842 </td><td> 3           </td><td>  8          </td><td> 258392      </td><td> 2.68247846  </td><td>0.00113765104</td></tr>
	<tr><td>RBBP8        </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0005110876 </td><td> 3           </td><td> 17          </td><td> 549083      </td><td> 1.63799207  </td><td>0.00015175756</td></tr>
	<tr><td>NARS1        </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0015493642 </td><td> 2           </td><td> 50          </td><td>1614950      </td><td> 0.03795826  </td><td>0.00003298030</td></tr>
	<tr><td>PIGN         </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0009451138 </td><td>18           </td><td> 31          </td><td>1001269      </td><td> 0.75107461  </td><td>0.00005527431</td></tr>
	<tr><td>FKRP         </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0021597293 </td><td>10           </td><td> 70          </td><td>2260930      </td><td>-0.44121593  </td><td>0.00002560865</td></tr>
	<tr><td>SCP2         </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0001703472 </td><td> 2           </td><td>  6          </td><td> 193794      </td><td> 3.22308443  </td><td>0.00309042695</td></tr>
	<tr><td>LEPR         </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0005419300 </td><td> 1           </td><td> 18          </td><td> 581382      </td><td> 1.55345615  </td><td>0.00013501234</td></tr>
	<tr><td>NCF4         </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0003408031 </td><td> 4           </td><td> 11          </td><td> 355289      </td><td> 2.22262420  </td><td>0.00043531750</td></tr>
	<tr><td>TRMU         </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0004493085 </td><td> 9           </td><td> 15          </td><td> 484485      </td><td> 1.82385631  </td><td>0.00019938153</td></tr>
	<tr><td>DNAJB2       </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0002633050 </td><td> 1           </td><td>  9          </td><td> 290691      </td><td> 2.59482774  </td><td>0.00078274817</td></tr>
	<tr><td>GLB1         </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0011615249 </td><td>29           </td><td> 38          </td><td>1227362      </td><td> 0.45361448  </td><td>0.00004322010</td></tr>
	<tr><td>FRAS1        </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0004336632 </td><td>16           </td><td> 14          </td><td> 452186      </td><td> 1.87498764  </td><td>0.00023412525</td></tr>
	<tr><td>BDP1         </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0004802030 </td><td> 1           </td><td> 16          </td><td> 516784      </td><td> 1.72791824  </td><td>0.00017269849</td></tr>
	<tr><td>CYP7B1       </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0014254749 </td><td>12           </td><td> 46          </td><td>1485754      </td><td> 0.15819193  </td><td>0.00003556323</td></tr>
	<tr><td>SURF1        </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0018428156 </td><td>11           </td><td> 60          </td><td>1937940      </td><td>-0.21227715  </td><td>0.00002848357</td></tr>
	<tr><td>ERCC6L2      </td><td>0.001590670  </td><td>1            </td><td>3            </td><td>5658         </td><td>0.0005576443 </td><td> 2           </td><td> 18          </td><td> 581382      </td><td> 1.51221749  </td><td>0.00013501234</td></tr>
	<tr><td>ATM          </td><td>0.001590669  </td><td>3            </td><td>3            </td><td>5658         </td><td>0.0017039015 </td><td>63           </td><td> 55          </td><td>1776445      </td><td>-0.09920828  </td><td>0.00003045638</td></tr>
	<tr><td>GNPTAB       </td><td>0.001590669  </td><td>3            </td><td>3            </td><td>5658         </td><td>0.0010997070 </td><td>19           </td><td> 36          </td><td>1162764      </td><td> 0.53251447  </td><td>0.00004595072</td></tr>
	<tr><td>DDX11        </td><td>0.001590669  </td><td>3            </td><td>3            </td><td>5658         </td><td>0.0003872345 </td><td> 7           </td><td> 13          </td><td> 419887      </td><td> 2.03835426  </td><td>0.00028052129</td></tr>
	<tr><td>BRCA2        </td><td>0.001590669  </td><td>3            </td><td>3            </td><td>5658         </td><td>0.0014874492 </td><td>60           </td><td> 48          </td><td>1550352      </td><td> 0.09679326  </td><td>0.00003419678</td></tr>
	<tr><td>GALC         </td><td>0.001590669  </td><td>3            </td><td>3            </td><td>5658         </td><td>0.0008378020 </td><td>23           </td><td> 27          </td><td> 872073      </td><td> 0.92495243  </td><td>0.00006698791</td></tr>
	<tr><td>POLG         </td><td>0.001590669  </td><td>3            </td><td>3            </td><td>5658         </td><td>0.0040879575 </td><td>27           </td><td>132          </td><td>4263468      </td><td>-1.36174654  </td><td>0.00001882540</td></tr>
	<tr><td>GPR179       </td><td>0.001590669  </td><td>3            </td><td>3            </td><td>5658         </td><td>0.0002478108 </td><td> 5           </td><td>  8          </td><td> 258392      </td><td> 2.68232269  </td><td>0.00113765104</td></tr>
	<tr><td>GCDH         </td><td>0.001590669  </td><td>3            </td><td>3            </td><td>5658         </td><td>0.0024368682 </td><td>40           </td><td> 79          </td><td>2551621      </td><td>-0.61539457  </td><td>0.00002379256</td></tr>
	<tr><td>ARSB         </td><td>0.001590669  </td><td>3            </td><td>3            </td><td>5658         </td><td>0.0009140548 </td><td>15           </td><td> 30          </td><td> 968970      </td><td> 0.79928109  </td><td>0.00005774196</td></tr>
	<tr><td>DNAH11       </td><td>0.001590669  </td><td>3            </td><td>3            </td><td>5658         </td><td>0.0013166912 </td><td>37           </td><td> 43          </td><td>1388857      </td><td> 0.27271662  </td><td>0.00003795792</td></tr>
</tbody>
</table>




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
tcu[order(tcu$cum_AF_fold_change, decreasing=F),]
```


<table>
<thead><tr><th scope=col>Symbol</th><th scope=col>cum_AF</th><th scope=col>nr_vars</th><th scope=col>cum_AC</th><th scope=col>cum_AN</th><th scope=col>cum_gnomad_AF_NFE</th><th scope=col>gnomad_nr_vars_NFE</th><th scope=col>cum_gnomad_AC_NFE</th><th scope=col>cum_gnomad_AN_NFE</th><th scope=col>cum_AF_fold_change</th><th scope=col>pvals</th></tr></thead>
<tbody>
	<tr><td>POLG                                </td><td>0.001590669                         </td><td>3                                   </td><td> 3                                  </td><td> 5658                               </td><td>0.0040879575                        </td><td>27                                  </td><td>132                                 </td><td>4263468                             </td><td>-1.361747                           </td><td>0.0000188253986940597513180875588468</td></tr>
	<tr><td>BCHE                                </td><td>0.002651116                         </td><td>3                                   </td><td> 5                                  </td><td> 9430                               </td><td>0.0061805733                        </td><td>15                                  </td><td>200                                 </td><td>6459800                             </td><td>-1.221141                           </td><td>0.0000000000069300189039680994598879</td></tr>
	<tr><td>SLC26A2                             </td><td>0.005832456                         </td><td>4                                   </td><td>11                                  </td><td>20746                               </td><td>0.0025259415                        </td><td> 7                                  </td><td> 82                                 </td><td>2648518                             </td><td> 1.207282                           </td><td>0.0000000000000000000000000003784498</td></tr>
	<tr><td>ADSL                                </td><td>0.002122016                         </td><td>3                                   </td><td> 4                                  </td><td> 7544                               </td><td>0.0008672026                        </td><td>11                                  </td><td> 28                                 </td><td> 904372                             </td><td> 1.290995                           </td><td>0.0000001282386477429749010698391883</td></tr>
	<tr><td>DDX11                               </td><td>0.001590669                         </td><td>3                                   </td><td> 3                                  </td><td> 5658                               </td><td>0.0003872345                        </td><td> 7                                  </td><td> 13                                 </td><td> 419887                             </td><td> 2.038354                           </td><td>0.0002805212883797294858533033501402</td></tr>
	<tr><td>CNGB3                               </td><td>0.003712690                         </td><td>3                                   </td><td> 7                                  </td><td>13202                               </td><td>0.0008219265                        </td><td>23                                  </td><td> 27                                 </td><td> 872073                             </td><td> 2.175384                           </td><td>0.0000000000000107425844857754702852</td></tr>
	<tr><td>PCDH15                              </td><td>0.002120892                         </td><td>4                                   </td><td> 4                                  </td><td> 7544                               </td><td>0.0004492958                        </td><td>18                                  </td><td> 15                                 </td><td> 484485                             </td><td> 2.238934                           </td><td>0.0000009835424129108588841083208112</td></tr>
	<tr><td>FANCD2                              </td><td>0.002121455                         </td><td>4                                   </td><td> 4                                  </td><td> 7544                               </td><td>0.0003871877                        </td><td> 7                                  </td><td> 13                                 </td><td> 419887                             </td><td> 2.453949                           </td><td>0.0000017760033602144646465160722840</td></tr>
	<tr><td>GPR179                              </td><td>0.001590669                         </td><td>3                                   </td><td> 3                                  </td><td> 5658                               </td><td>0.0002478108                        </td><td> 5                                  </td><td>  8                                 </td><td> 258392                             </td><td> 2.682323                           </td><td>0.0011376510433942635999032821914057</td></tr>
	<tr><td>PIGV                                </td><td>0.004241786                         </td><td>3                                   </td><td> 8                                  </td><td>15088                               </td><td>0.0004646057                        </td><td> 6                                  </td><td> 15                                 </td><td> 484485                             </td><td> 3.190593                           </td><td>0.0000000000000532085557900878768625</td></tr>
	<tr><td>ABCB11                              </td><td>0.003711566                         </td><td>4                                   </td><td> 7                                  </td><td>13202                               </td><td>0.0004028471                        </td><td>12                                  </td><td> 13                                 </td><td> 419887                             </td><td> 3.203724                           </td><td>0.0000000000087145363968780358734649</td></tr>
	<tr><td>XPA                                 </td><td>0.003711563                         </td><td>3                                   </td><td> 7                                  </td><td>13202                               </td><td>0.0003717554                        </td><td> 8                                  </td><td> 12                                 </td><td> 387588                             </td><td> 3.319601                           </td><td>0.0000000000205469765423096345832060</td></tr>
</tbody>
</table>



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
    geom_abline() + xlim(0, 0.008) + ylim(0, 0.008) +
    ggtitle("Cumulative allele frequencies of selected ClinVar variants in AR genes")+
    ylab("Per-gene cumulative allele frequencies in GnomAD NFE population")+
    xlab("Per-gene cumulative allele frequencies in Polish population")

p<- p +   geom_label_repel(
  aes(
    fill = "red", 
    label = Symbol),
  fontface = 'bold.italic', 
    color = 'white',
  size = 3,
     segment.color = "black",
    #max.overlaps = 0,
    min.segment.length = 0,
  box.padding = unit(0.6, "lines"),
 point.padding = unit(0.6, "lines")
)
p +  theme(legend.position = "none") + theme(plot.title = element_text(size=13))
```


    
![png](output_37_0.png)
    



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
allPanelsDT3[,cum_AN_sum:=round(max(as.numeric(tgp_vep_nonsyn_rare$AN))*allPanelsDT3$cum_AC_sum)]
allPanelsDT3[,nr_vars_sum:=sum(nr_vars, na.rm=T), by=Level4]

allPanelsDT3[,cum_gnomad_AF_NFE_sum:=sum(cum_gnomad_AF_NFE, na.rm=T), by=Level4]
allPanelsDT3[,cum_gnomad_AC_NFE_sum:=round(allPanelsDT3$cum_gnomad_AF_NFE_sum*nfe_cnt)]
allPanelsDT3[,cum_gnomad_AN_NFE_sum:=round(nfe_cnt*allPanelsDT3$cum_gnomad_AC_NFE_sum)]
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
length(which(allPanelsDT3_unique$pvals < 0.01))
```


176



```R
SuppTable_recessive_panels <- allPanelsDT3_unique[which(allPanelsDT3_unique$pvals < 0.01),]
```


```R
SuppTable_recessive_panels[,c("Level4", "cum_AF_sum", "cum_gnomad_AF_NFE_sum", "cum_AF_fold_change","nr_vars_sum" ,"gnomad_nr_vars_NFE_sum" , "pvals")]
```


<table>
<thead><tr><th scope=col>Level4</th><th scope=col>cum_AF_sum</th><th scope=col>cum_gnomad_AF_NFE_sum</th><th scope=col>cum_AF_fold_change</th><th scope=col>nr_vars_sum</th><th scope=col>gnomad_nr_vars_NFE_sum</th><th scope=col>pvals</th></tr></thead>
<tbody>
	<tr><td>Adult onset movement disorder                                                                                                                                                                                                         </td><td>0.014316020                                                                                                                                                                                                                           </td><td>0.0132885165                                                                                                                                                                                                                          </td><td> 0.10745041                                                                                                                                                                                                                           </td><td>14                                                                                                                                                                                                                                    </td><td>189                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000000000000000000037555308403016628724156794984099641021960792348092790295909258901894060667008021429009481303228470799242150211568350819052068037806808849708582257598</td></tr>
	<tr><td>Adult solid tumours cancer susceptibility                                                                                                                                                                                             </td><td>0.015380422                                                                                                                                                                                                                           </td><td>0.0086585732                                                                                                                                                                                                                          </td><td> 0.82889387                                                                                                                                                                                                                           </td><td>19                                                                                                                                                                                                                                    </td><td>207                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000000000000000000000036781854179258396522588198464989495844526699396569354138559617160590736330131306051655601821482717778324898036783843654312108694075657438026727096</td></tr>
	<tr><td>Adult solid tumours for rare disease                                                                                                                                                                                                  </td><td>0.001590669                                                                                                                                                                                                                           </td><td>0.0017039015                                                                                                                                                                                                                          </td><td>-0.09920828                                                                                                                                                                                                                           </td><td> 3                                                                                                                                                                                                                                    </td><td> 63                                                                                                                                                                                                                                   </td><td>0.000016173837218752095122890552469918645783764077350497245788574218750000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Albinism or congenital nystagmus                                                                                                                                                                                                      </td><td>0.002651678                                                                                                                                                                                                                           </td><td>0.0023085035                                                                                                                                                                                                                          </td><td> 0.19994768                                                                                                                                                                                                                           </td><td> 5                                                                                                                                                                                                                                    </td><td> 28                                                                                                                                                                                                                                   </td><td>0.000000000012104799394219024061535360466116724312626851123297910817200317978858947753906250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Amelogenesis imperfecta                                                                                                                                                                                                               </td><td>0.003711561                                                                                                                                                                                                                           </td><td>0.0049098637                                                                                                                                                                                                                          </td><td>-0.40365690                                                                                                                                                                                                                           </td><td> 7                                                                                                                                                                                                                                    </td><td> 63                                                                                                                                                                                                                                   </td><td>0.000000000000000002883752867618456435439732445467761846608440243081364351418471869692439213395118713378906250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Anophthalmia or microphthalmia                                                                                                                                                                                                        </td><td>0.003181339                                                                                                                                                                                                                           </td><td>0.0022306231                                                                                                                                                                                                                          </td><td> 0.51218734                                                                                                                                                                                                                           </td><td> 4                                                                                                                                                                                                                                    </td><td> 38                                                                                                                                                                                                                                   </td><td>0.000000000000016555995411505194787940781865345999072106292562822638103625649819150567054748535156250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Arthrogryposis                                                                                                                                                                                                                        </td><td>0.021749362                                                                                                                                                                                                                           </td><td>0.0154366620                                                                                                                                                                                                                          </td><td> 0.49461226                                                                                                                                                                                                                           </td><td>30                                                                                                                                                                                                                                    </td><td>254                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000197570975369675683457631971719557921212514553976379705349352658934885010047627459810518557702018583407675265</td></tr>
	<tr><td>Ataxia and cerebellar anomalies - narrow panel                                                                                                                                                                                        </td><td>0.027042527                                                                                                                                                                                                                           </td><td>0.0206442645                                                                                                                                                                                                                          </td><td> 0.38948895                                                                                                                                                                                                                           </td><td>30                                                                                                                                                                                                                                    </td><td>339                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000009804690728283747512154757338582234553769431544910066551316104668809040962408</td></tr>
	<tr><td>Autosomal recessive congenital ichthyosis                                                                                                                                                                                             </td><td>0.012195121                                                                                                                                                                                                                           </td><td>0.0085041907                                                                                                                                                                                                                          </td><td> 0.52005822                                                                                                                                                                                                                           </td><td>10                                                                                                                                                                                                                                    </td><td> 92                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000002060678844001047663391653636805566884873937499473301875782715767676013816741819173941499690092599294825782884138888151775079373408226326970661212130825491352466761</td></tr>
	<tr><td>Autosomal recessive primary hypertrophic osteoarthropathy                                                                                                                                                                             </td><td>0.003181333                                                                                                                                                                                                                           </td><td>0.0002788313                                                                                                                                                                                                                          </td><td> 3.51216697                                                                                                                                                                                                                           </td><td> 2                                                                                                                                                                                                                                    </td><td>  2                                                                                                                                                                                                                                   </td><td>0.000000004765254926611347604275401897014688423759309898741776123642921447753906250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Bardet Biedl syndrome                                                                                                                                                                                                                 </td><td>0.003711565                                                                                                                                                                                                                           </td><td>0.0021845153                                                                                                                                                                                                                          </td><td> 0.76471442                                                                                                                                                                                                                           </td><td> 6                                                                                                                                                                                                                                    </td><td> 44                                                                                                                                                                                                                                   </td><td>0.000000000000000026128240194603115132413952981288155052041673819082819030246156444263760931789875030517578125000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Bleeding and platelet disorders                                                                                                                                                                                                       </td><td>0.029694749                                                                                                                                                                                                                           </td><td>0.0156659331                                                                                                                                                                                                                          </td><td> 0.92257714                                                                                                                                                                                                                           </td><td>16                                                                                                                                                                                                                                    </td><td>104                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000722286440903274013101672783784742697669313299295042998636985759009</td></tr>
	<tr><td>Brain cancer pertinent cancer susceptibility                                                                                                                                                                                          </td><td>0.001590669                                                                                                                                                                                                                           </td><td>0.0017039015                                                                                                                                                                                                                          </td><td>-0.09920828                                                                                                                                                                                                                           </td><td> 3                                                                                                                                                                                                                                    </td><td> 63                                                                                                                                                                                                                                   </td><td>0.000016173837218752095122890552469918645783764077350497245788574218750000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Brain channelopathy                                                                                                                                                                                                                   </td><td>0.006362669                                                                                                                                                                                                                           </td><td>0.0047390235                                                                                                                                                                                                                          </td><td> 0.42504226                                                                                                                                                                                                                           </td><td> 4                                                                                                                                                                                                                                    </td><td> 63                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000005225781953614354727987798244791176227089202151944897112451339790214380594897994021866272085219407017575576901435852050781250000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>CAKUT                                                                                                                                                                                                                                 </td><td>0.018557802                                                                                                                                                                                                                           </td><td>0.0133037085                                                                                                                                                                                                                          </td><td> 0.48019738                                                                                                                                                                                                                           </td><td> 8                                                                                                                                                                                                                                    </td><td> 68                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000334867419758338423988527002609465902459439472346790597473244530345124744469327496195079679262532911611062619685379928384385608</td></tr>
	<tr><td>Cardiomyopathies - including childhood onset                                                                                                                                                                                          </td><td>0.022270501                                                                                                                                                                                                                           </td><td>0.0215378402                                                                                                                                                                                                                          </td><td> 0.04826043                                                                                                                                                                                                                           </td><td>25                                                                                                                                                                                                                                    </td><td>246                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000469595228482859193737929449762615374306307411874075289664520201623422673236159986849240500470046309131</td></tr>
	<tr><td>Cataracts                                                                                                                                                                                                                             </td><td>0.034998669                                                                                                                                                                                                                           </td><td>0.0233262338                                                                                                                                                                                                                          </td><td> 0.58534667                                                                                                                                                                                                                           </td><td>22                                                                                                                                                                                                                                    </td><td>192                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000116721907454689719287010025148528</td></tr>
	<tr><td>Cerebellar hypoplasia                                                                                                                                                                                                                 </td><td>0.004772579                                                                                                                                                                                                                           </td><td>0.0043922801                                                                                                                                                                                                                          </td><td> 0.11979902                                                                                                                                                                                                                           </td><td> 5                                                                                                                                                                                                                                    </td><td> 45                                                                                                                                                                                                                                   </td><td>0.000000000000000000000003795916264804964085781877310157099399442507344722576253227922524549620142408912215614691376686096191406250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Cerebral malformations                                                                                                                                                                                                                </td><td>0.018558924                                                                                                                                                                                                                           </td><td>0.0164598456                                                                                                                                                                                                                          </td><td> 0.17316227                                                                                                                                                                                                                           </td><td> 9                                                                                                                                                                                                                                    </td><td>122                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000010099949266997450373626822555299265149604770296436513622391216973412681295355352302059880122622400056359111684097882047679401</td></tr>
	<tr><td>Childhood onset dystonia or chorea or related movement disorder                                                                                                                                                                       </td><td>0.032344168                                                                                                                                                                                                                           </td><td>0.0298752877                                                                                                                                                                                                                          </td><td> 0.11455300                                                                                                                                                                                                                           </td><td>35                                                                                                                                                                                                                                    </td><td>400                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000011204187883477117529597489724674195452439151</td></tr>
	<tr><td>Childhood solid tumours cancer susceptibility                                                                                                                                                                                         </td><td>0.026554361                                                                                                                                                                                                                           </td><td>0.0151516815                                                                                                                                                                                                                          </td><td> 0.80947090                                                                                                                                                                                                                           </td><td>33                                                                                                                                                                                                                                    </td><td>409                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000037298817057614986941266011292723907616888087781141656148713589917610759174203785888</td></tr>
	<tr><td>Cholestasis                                                                                                                                                                                                                           </td><td>0.025981498                                                                                                                                                                                                                           </td><td>0.0181832150                                                                                                                                                                                                                          </td><td> 0.51487731                                                                                                                                                                                                                           </td><td>24                                                                                                                                                                                                                                    </td><td>166                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000118846054614141391352540036340501016466545471251299802198557989989534236307242830226</td></tr>
	<tr><td>Clefting                                                                                                                                                                                                                              </td><td>0.040300914                                                                                                                                                                                                                           </td><td>0.0207870621                                                                                                                                                                                                                          </td><td> 0.95512669                                                                                                                                                                                                                           </td><td>23                                                                                                                                                                                                                                    </td><td>152                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000009390199</td></tr>
	<tr><td>ClinGen Gene Validity Curations                                                                                                                                                                                                       </td><td>0.009014923                                                                                                                                                                                                                           </td><td>0.0102306951                                                                                                                                                                                                                          </td><td>-0.18251709                                                                                                                                                                                                                           </td><td>14                                                                                                                                                                                                                                    </td><td>203                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000132390399378602836069984394194444173746435185961947406893396751808249921614008187998791193618473723559371579877127105755678761547788901964395336108282208442687988281250000000000000</td></tr>
	<tr><td>Confirmed Fanconi anaemia or Bloom syndrome                                                                                                                                                                                           </td><td>0.013786366                                                                                                                                                                                                                           </td><td>0.0074970420                                                                                                                                                                                                                          </td><td> 0.87884883                                                                                                                                                                                                                           </td><td>15                                                                                                                                                                                                                                    </td><td>218                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000000000000064081954705274443209048627764490722379606042828850523901504913182042661512537190692352794121364656827980801049063005302125523997317709292844106503850354233</td></tr>
	<tr><td>Congenital adrenal hypoplasia                                                                                                                                                                                                         </td><td>0.001590673                                                                                                                                                                                                                           </td><td>0.0013318314                                                                                                                                                                                                                          </td><td> 0.25622583                                                                                                                                                                                                                           </td><td> 2                                                                                                                                                                                                                                    </td><td> 16                                                                                                                                                                                                                                   </td><td>0.000020157525990592170225031848351804342200921382755041122436523437500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Congenital disorders of glycosylation                                                                                                                                                                                                 </td><td>0.021741977                                                                                                                                                                                                                           </td><td>0.0164135756                                                                                                                                                                                                                          </td><td> 0.40559358                                                                                                                                                                                                                           </td><td>21                                                                                                                                                                                                                                    </td><td>129                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000056406852778347270208210971646828920927140807619419633822771017768296254054851813274067139515721161490409389</td></tr>
	<tr><td>Congenital hyperinsulinism                                                                                                                                                                                                            </td><td>0.003711563                                                                                                                                                                                                                           </td><td>0.0070149095                                                                                                                                                                                                                          </td><td>-0.91839763                                                                                                                                                                                                                           </td><td> 3                                                                                                                                                                                                                                    </td><td> 27                                                                                                                                                                                                                                   </td><td>0.000000000000000001634953432640934863293760897077792265436133125465652509572533901405222422908991575241088867187500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Congenital hypothyroidism                                                                                                                                                                                                             </td><td>0.011669978                                                                                                                                                                                                                           </td><td>0.0105178604                                                                                                                                                                                                                          </td><td> 0.14996059                                                                                                                                                                                                                           </td><td>11                                                                                                                                                                                                                                    </td><td> 74                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000000000000000000000000000000000000000000167329452467634644134209242923408726440349079036161083639992462498066115881703674312555117166994563792127128111929486534334443120644213198189669228127170479050039864</td></tr>
	<tr><td>Congenital muscular dystrophy                                                                                                                                                                                                         </td><td>0.005304486                                                                                                                                                                                                                           </td><td>0.0038509392                                                                                                                                                                                                                          </td><td> 0.46200262                                                                                                                                                                                                                           </td><td> 8                                                                                                                                                                                                                                    </td><td> 70                                                                                                                                                                                                                                   </td><td>0.000000000000000000000000008225717000888006478159351027724376482097435658069942454501160845576430893227870555506342498119920492172241210937500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td><td>...</td></tr>
	<tr><td>Rare multisystem ciliopathy disorders                                                                                                                                                                                                                                                                                            </td><td>0.039767289                                                                                                                                                                                                                                                                                                                      </td><td>0.0334573980                                                                                                                                                                                                                                                                                                                     </td><td> 0.2492570                                                                                                                                                                                                                                                                                                                       </td><td> 40                                                                                                                                                                                                                                                                                                                              </td><td> 464                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000086625612176492735451195438283762886156005890687307254926243670293727610200124522554653117553</td></tr>
	<tr><td>Rare multisystem ciliopathy Super panel                                                                                                                                                                                                                                                                                          </td><td>0.106576510                                                                                                                                                                                                                                                                                                                      </td><td>0.0945028298                                                                                                                                                                                                                                                                                                                     </td><td> 0.1734601                                                                                                                                                                                                                                                                                                                       </td><td> 99                                                                                                                                                                                                                                                                                                                              </td><td>1116                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Renal ciliopathies                                                                                                                                                                                                                                                                                                               </td><td>0.036585947                                                                                                                                                                                                                                                                                                                      </td><td>0.0312489437                                                                                                                                                                                                                                                                                                                     </td><td> 0.2274822                                                                                                                                                                                                                                                                                                                       </td><td> 35                                                                                                                                                                                                                                                                                                                              </td><td> 389                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000177525681504626013514480209151242995487793956044673495894799377834115570642906591816067101675714662418666148940</td></tr>
	<tr><td>Renal superpanel - broad                                                                                                                                                                                                                                                                                                         </td><td>0.172333812                                                                                                                                                                                                                                                                                                                      </td><td>0.1237773353                                                                                                                                                                                                                                                                                                                     </td><td> 0.4774586                                                                                                                                                                                                                                                                                                                       </td><td>177                                                                                                                                                                                                                                                                                                                              </td><td>1715                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Renal superpanel - narrow                                                                                                                                                                                                                                                                                                        </td><td>0.098098618                                                                                                                                                                                                                                                                                                                      </td><td>0.0705301415                                                                                                                                                                                                                                                                                                                     </td><td> 0.4759929                                                                                                                                                                                                                                                                                                                       </td><td> 97                                                                                                                                                                                                                                                                                                                              </td><td> 901                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Renal tubulopathies                                                                                                                                                                                                                                                                                                              </td><td>0.008487533                                                                                                                                                                                                                                                                                                                      </td><td>0.0059804221                                                                                                                                                                                                                                                                                                                     </td><td> 0.5050980                                                                                                                                                                                                                                                                                                                       </td><td>  9                                                                                                                                                                                                                                                                                                                              </td><td>  76                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000043324491031348012649376239905780254474847107326238033603606868289214533622521395899915697152233495824337015714096225937890238810723531059920787811279296875000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Respiratory ciliopathies including non-CF bronchiectasis                                                                                                                                                                                                                                                                         </td><td>0.016436911                                                                                                                                                                                                                                                                                                                      </td><td>0.0246401653                                                                                                                                                                                                                                                                                                                     </td><td>-0.5840727                                                                                                                                                                                                                                                                                                                       </td><td> 10                                                                                                                                                                                                                                                                                                                              </td><td> 129                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000984809685693337251810195194101009499875370316647467259477121034045320973432091583501094871467337636013067912265608311910645128777418016600849882764328418896350983621424107147400210238973376198524187374240817394188576747348928</td></tr>
	<tr><td>Retinal disorders                                                                                                                                                                                                                                                                                                                </td><td>0.054615794                                                                                                                                                                                                                                                                                                                      </td><td>0.0447445225                                                                                                                                                                                                                                                                                                                     </td><td> 0.2876071                                                                                                                                                                                                                                                                                                                       </td><td> 71                                                                                                                                                                                                                                                                                                                              </td><td> 827                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001391319</td></tr>
	<tr><td>Rhabdomyolysis and metabolic muscle disorders                                                                                                                                                                                                                                                                                    </td><td>0.018558371                                                                                                                                                                                                                                                                                                                      </td><td>0.0244241076                                                                                                                                                                                                                                                                                                                     </td><td>-0.3962358                                                                                                                                                                                                                                                                                                                       </td><td> 21                                                                                                                                                                                                                                                                                                                              </td><td> 197                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000695621126262151012862471854941657605657055612974205231418904837509752124267452169227990284198085149707541647019964608991188915095507459847865547839875999282625864603834297740499792778439932659634519479774156371662</td></tr>
	<tr><td>Sarcoma cancer susceptibility                                                                                                                                                                                                                                                                                                    </td><td>0.001629922                                                                                                                                                                                                                                                                                                                      </td><td>0.0013490120                                                                                                                                                                                                                                                                                                                     </td><td> 0.2728997                                                                                                                                                                                                                                                                                                                       </td><td>  3                                                                                                                                                                                                                                                                                                                              </td><td>  32                                                                                                                                                                                                                                                                                                                             </td><td>0.0000197061427672924177048260907518795193027472123503684997558593750000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Sarcoma susceptibility                                                                                                                                                                                                                                                                                                           </td><td>0.001629922                                                                                                                                                                                                                                                                                                                      </td><td>0.0013490120                                                                                                                                                                                                                                                                                                                     </td><td> 0.2728997                                                                                                                                                                                                                                                                                                                       </td><td>  3                                                                                                                                                                                                                                                                                                                              </td><td>  32                                                                                                                                                                                                                                                                                                                             </td><td>0.0000197061427672924177048260907518795193027472123503684997558593750000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Severe early-onset obesity                                                                                                                                                                                                                                                                                                       </td><td>0.004241785                                                                                                                                                                                                                                                                                                                      </td><td>0.0033936906                                                                                                                                                                                                                                                                                                                     </td><td> 0.3218165                                                                                                                                                                                                                                                                                                                       </td><td>  6                                                                                                                                                                                                                                                                                                                              </td><td>  94                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000076970713691856802244469132097713447218288554741584038087869046318445498400251381099224090576171875000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Severe microcephaly                                                                                                                                                                                                                                                                                                              </td><td>0.036594405                                                                                                                                                                                                                                                                                                                      </td><td>0.0238246689                                                                                                                                                                                                                                                                                                                     </td><td> 0.6191669                                                                                                                                                                                                                                                                                                                       </td><td> 30                                                                                                                                                                                                                                                                                                                              </td><td> 331                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000002185347719970120189398624580144407441993205829259923382410710735469178685989587143242592397121232802166070411549371</td></tr>
	<tr><td>Severe Paediatric Disorders                                                                                                                                                                                                                                                                                                      </td><td>0.506858543                                                                                                                                                                                                                                                                                                                      </td><td>0.3307293722                                                                                                                                                                                                                                                                                                                     </td><td> 0.6159320                                                                                                                                                                                                                                                                                                                       </td><td>476                                                                                                                                                                                                                                                                                                                              </td><td>3512                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Skeletal ciliopathies                                                                                                                                                                                                                                                                                                            </td><td>0.024390251                                                                                                                                                                                                                                                                                                                      </td><td>0.0192126648                                                                                                                                                                                                                                                                                                                     </td><td> 0.3442470                                                                                                                                                                                                                                                                                                                       </td><td> 13                                                                                                                                                                                                                                                                                                                              </td><td> 123                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000102908896991959155337257009004505722182725461855860341022304437311603438257919335329890830474557282317122069330175002299223047133453833922115321687283581982685609404203753915492711278</td></tr>
	<tr><td>Skeletal dysplasia                                                                                                                                                                                                                                                                                                               </td><td>0.081173537                                                                                                                                                                                                                                                                                                                      </td><td>0.0528609157                                                                                                                                                                                                                                                                                                                     </td><td> 0.6188081                                                                                                                                                                                                                                                                                                                       </td><td> 77                                                                                                                                                                                                                                                                                                                              </td><td> 585                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Structural basal ganglia disorders                                                                                                                                                                                                                                                                                               </td><td>0.014316016                                                                                                                                                                                                                                                                                                                      </td><td>0.0127047505                                                                                                                                                                                                                                                                                                                     </td><td> 0.1722620                                                                                                                                                                                                                                                                                                                       </td><td> 15                                                                                                                                                                                                                                                                                                                              </td><td> 171                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000625045279894423731739461066983464342442183055740720470379270327455563326246879256566026981445281613849796013682328771468673824921996314243081040085092246810534843515074429335250154248755132593015973463401735443767393007874488830566406250000</td></tr>
	<tr><td>Structural eye disease                                                                                                                                                                                                                                                                                                           </td><td>0.030331404                                                                                                                                                                                                                                                                                                                      </td><td>0.0180514988                                                                                                                                                                                                                                                                                                                     </td><td> 0.7486937                                                                                                                                                                                                                                                                                                                       </td><td> 24                                                                                                                                                                                                                                                                                                                              </td><td> 230                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000085623940501646025235002681143334082258623932122535847573634425254157732118096683721847832084721431909815752234252585151367363010496808105639555304191482</td></tr>
	<tr><td>Thoracic aortic aneurysm and dissection                                                                                                                                                                                                                                                                                          </td><td>0.013357463                                                                                                                                                                                                                                                                                                                      </td><td>0.0010976062                                                                                                                                                                                                                                                                                                                     </td><td> 3.6052136                                                                                                                                                                                                                                                                                                                       </td><td>  3                                                                                                                                                                                                                                                                                                                              </td><td>   7                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000014982785765464259047821898554778986222514126794275740189323375967471438120381023933933959871647516770398625014342744096729198588491271948441863059997558593750000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Thoracic dystrophies                                                                                                                                                                                                                                                                                                             </td><td>0.003181901                                                                                                                                                                                                                                                                                                                      </td><td>0.0009950244                                                                                                                                                                                                                                                                                                                     </td><td> 1.6770851                                                                                                                                                                                                                                                                                                                       </td><td>  6                                                                                                                                                                                                                                                                                                                              </td><td>  25                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000003835202872046562123580721629929505284302936052309718206743127666413784027099609375000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Tubulointerstitial kidney disease                                                                                                                                                                                                                                                                                                </td><td>0.003711565                                                                                                                                                                                                                                                                                                                      </td><td>0.0021380095                                                                                                                                                                                                                                                                                                                     </td><td> 0.7957594                                                                                                                                                                                                                                                                                                                       </td><td>  6                                                                                                                                                                                                                                                                                                                              </td><td>  39                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000291500181873259336862922550524555490593429039056835166965342409639561083167791366577148437500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Tumour predisposition - childhood onset                                                                                                                                                                                                                                                                                          </td><td>0.026024138                                                                                                                                                                                                                                                                                                                      </td><td>0.0149348766                                                                                                                                                                                                                                                                                                                     </td><td> 0.8011651                                                                                                                                                                                                                                                                                                                       </td><td> 32                                                                                                                                                                                                                                                                                                                              </td><td> 407                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000308042883075400485560607033114271873588062766993619800061846805271362314743216556110582208693803412639423456618962618385864100437452154201446614533157080637679852436623164237068</td></tr>
	<tr><td>Undiagnosed metabolic disorders                                                                                                                                                                                                                                                                                                  </td><td>0.219637974                                                                                                                                                                                                                                                                                                                      </td><td>0.1575700524                                                                                                                                                                                                                                                                                                                     </td><td> 0.4791341                                                                                                                                                                                                                                                                                                                       </td><td>219                                                                                                                                                                                                                                                                                                                              </td><td>1685                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Unexplained kidney failure in young people                                                                                                                                                                                                                                                                                       </td><td>0.021740854                                                                                                                                                                                                                                                                                                                      </td><td>0.0134298829                                                                                                                                                                                                                                                                                                                     </td><td> 0.6949619                                                                                                                                                                                                                                                                                                                       </td><td> 27                                                                                                                                                                                                                                                                                                                              </td><td> 270                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000046548333261152926607714585489309812172330550868693039569527940035935811415899186007173420603021910886426574441982097042764463861692085370857473874336286744174302343704290979046567179123795347088920096</td></tr>
	<tr><td>Unexplained paediatric onset end-stage renal disease                                                                                                                                                                                                                                                                             </td><td>0.052494340                                                                                                                                                                                                                                                                                                                      </td><td>0.0398173109                                                                                                                                                                                                                                                                                                                     </td><td> 0.3987661                                                                                                                                                                                                                                                                                                                       </td><td> 53                                                                                                                                                                                                                                                                                                                              </td><td> 544                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000113192730409935993850</td></tr>
	<tr><td>Vascular skin disorders                                                                                                                                                                                                                                                                                                          </td><td>0.003181342                                                                                                                                                                                                                                                                                                                      </td><td>0.0023078478                                                                                                                                                                                                                                                                                                                     </td><td> 0.4630874                                                                                                                                                                                                                                                                                                                       </td><td>  5                                                                                                                                                                                                                                                                                                                              </td><td>  74                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000148239707717852143350469277834994704558658866477838245145903783850371837615966796875000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>White matter disorders - adult onset                                                                                                                                                                                                                                                                                             </td><td>0.044009086                                                                                                                                                                                                                                                                                                                      </td><td>0.0325754208                                                                                                                                                                                                                                                                                                                     </td><td> 0.4340176                                                                                                                                                                                                                                                                                                                       </td><td> 45                                                                                                                                                                                                                                                                                                                              </td><td> 393                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000003122102034962548857220073037273400508718320743188420381040376115769872</td></tr>
	<tr><td>White matter disorders and cerebral calcification - narrow panel                                                                                                                                                                                                                                                                 </td><td>0.055147172                                                                                                                                                                                                                                                                                                                      </td><td>0.0334932075                                                                                                                                                                                                                                                                                                                     </td><td> 0.7194184                                                                                                                                                                                                                                                                                                                       </td><td> 58                                                                                                                                                                                                                                                                                                                              </td><td> 414                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000017776550325</td></tr>
	<tr><td>White matter disorders - childhood onset                                                                                                                                                                                                                                                                                         </td><td>0.509803012                                                                                                                                                                                                                                                                                                                      </td><td>0.3489477862                                                                                                                                                                                                                                                                                                                     </td><td> 0.5469287                                                                                                                                                                                                                                                                                                                       </td><td>526                                                                                                                                                                                                                                                                                                                              </td><td>4104                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Xeroderma pigmentosum, Trichothiodystrophy or Cockayne syndrome                                                                                                                                                                                                                                                                  </td><td>0.009017185                                                                                                                                                                                                                                                                                                                      </td><td>0.0042135425                                                                                                                                                                                                                                                                                                                     </td><td> 1.0976434                                                                                                                                                                                                                                                                                                                       </td><td> 10                                                                                                                                                                                                                                                                                                                              </td><td>  78                                                                                                                                                                                                                                                                                                                             </td><td>0.0000000000000000000000000000000000000000000002598585381139261369133051724813083606193893194943431392339441920181798567731362020471161078540714462817944744836974428625619992772044497542083263397216796875000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
</tbody>
</table>




```R
selPanels <- allPanelsDT3_unique[which(abs(cum_AF_fold_change) > 1 & pvals < 0.01 & nr_vars_sum>=3 ),]
selPanels[order(cum_AF_fold_change, decreasing=T), c("Level4", "cum_AF_sum", "cum_gnomad_AF_NFE_sum", "cum_AF_fold_change","nr_vars_sum" ,"gnomad_nr_vars_NFE_sum" , "pvals"),with=F]
```


<table>
<thead><tr><th scope=col>Level4</th><th scope=col>cum_AF_sum</th><th scope=col>cum_gnomad_AF_NFE_sum</th><th scope=col>cum_AF_fold_change</th><th scope=col>nr_vars_sum</th><th scope=col>gnomad_nr_vars_NFE_sum</th><th scope=col>pvals</th></tr></thead>
<tbody>
	<tr><td>Optic neuropathy                                                                                                                                                                                                                                                                                               </td><td>0.005314073                                                                                                                                                                                                                                                                                                    </td><td>0.0002788392                                                                                                                                                                                                                                                                                                   </td><td> 4.252313                                                                                                                                                                                                                                                                                                      </td><td> 3                                                                                                                                                                                                                                                                                                             </td><td>  5                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000003104548656263312302509563807338213397408153981782419350565760396420955657958984375000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Thoracic aortic aneurysm and dissection                                                                                                                                                                                                                                                                        </td><td>0.013357463                                                                                                                                                                                                                                                                                                    </td><td>0.0010976062                                                                                                                                                                                                                                                                                                   </td><td> 3.605214                                                                                                                                                                                                                                                                                                      </td><td> 3                                                                                                                                                                                                                                                                                                             </td><td>  7                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000000000000000000000000000000000014982785765464259047821898554778986222514126794275740189323375967471438120381023933933959871647516770398625014342744096729198588491271948441863059997558593750000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Intestinal failure                                                                                                                                                                                                                                                                                             </td><td>0.003711556                                                                                                                                                                                                                                                                                                    </td><td>0.0004802474                                                                                                                                                                                                                                                                                                   </td><td> 2.950174                                                                                                                                                                                                                                                                                                      </td><td> 3                                                                                                                                                                                                                                                                                                             </td><td> 14                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000005532505564278401040702824804962014197283123873027932404511375352740287780761718750000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Ehlers Danlos syndromes                                                                                                                                                                                                                                                                                        </td><td>0.005834704                                                                                                                                                                                                                                                                                                    </td><td>0.0010560352                                                                                                                                                                                                                                                                                                   </td><td> 2.466002                                                                                                                                                                                                                                                                                                      </td><td> 3                                                                                                                                                                                                                                                                                                             </td><td>  8                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000000000000031866996039611835578004745171833990206648939174807824448359394316036830030469673147308640182018280029296875000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Peeling skin syndrome                                                                                                                                                                                                                                                                                          </td><td>0.017497373                                                                                                                                                                                                                                                                                                    </td><td>0.0043831040                                                                                                                                                                                                                                                                                                   </td><td> 1.997114                                                                                                                                                                                                                                                                                                      </td><td> 3                                                                                                                                                                                                                                                                                                             </td><td>  7                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000004076126021137917666636908333517177925249307191327034792533021781033672760546774728999552240486888845343899318007825811803071409464532746627147489469914673529784724948774549460064765325228977975976188649382021822464850</td></tr>
	<tr><td>Gastrointestinal epithelial barrier disorders                                                                                                                                                                                                                                                                  </td><td>0.003712124                                                                                                                                                                                                                                                                                                    </td><td>0.0010537867                                                                                                                                                                                                                                                                                                   </td><td> 1.816662                                                                                                                                                                                                                                                                                                      </td><td> 7                                                                                                                                                                                                                                                                                                             </td><td> 21                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000011164997174417336983211579572288610058044064235818759733831484481925144791603088378906250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Thoracic dystrophies                                                                                                                                                                                                                                                                                           </td><td>0.003181901                                                                                                                                                                                                                                                                                                    </td><td>0.0009950244                                                                                                                                                                                                                                                                                                   </td><td> 1.677085                                                                                                                                                                                                                                                                                                      </td><td> 6                                                                                                                                                                                                                                                                                                             </td><td> 25                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000003835202872046562123580721629929505284302936052309718206743127666413784027099609375000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Pain syndromes                                                                                                                                                                                                                                                                                                 </td><td>0.001591796                                                                                                                                                                                                                                                                                                    </td><td>0.0005297068                                                                                                                                                                                                                                                                                                   </td><td> 1.587390                                                                                                                                                                                                                                                                                                      </td><td> 3                                                                                                                                                                                                                                                                                                             </td><td> 11                                                                                                                                                                                                                                                                                                            </td><td>0.0000805907370580582490355789349045778635627357289195060729980468750000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Epidermolysis bullosa                                                                                                                                                                                                                                                                                          </td><td>0.018027592                                                                                                                                                                                                                                                                                                    </td><td>0.0060248865                                                                                                                                                                                                                                                                                                   </td><td> 1.581201                                                                                                                                                                                                                                                                                                      </td><td> 5                                                                                                                                                                                                                                                                                                             </td><td> 38                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000361529997506698848115888817342206856162149267153029003714489186031111094686084301130856188219709289327386931781939261499448487022206677025445134699607909762767800198760160024401676016904069027655497182908774703</td></tr>
	<tr><td>Congenital myopathy                                                                                                                                                                                                                                                                                            </td><td>0.007953346                                                                                                                                                                                                                                                                                                    </td><td>0.0027123761                                                                                                                                                                                                                                                                                                   </td><td> 1.552005                                                                                                                                                                                                                                                                                                      </td><td> 5                                                                                                                                                                                                                                                                                                             </td><td> 56                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000000000000000000000000000079160970164403256871994350685601297327963422373612170193187371197001798300339079353787914939320714737686079942591277358587831258773803710937500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Epidermolysis bullosa and congenital skin fragility                                                                                                                                                                                                                                                            </td><td>0.020148488                                                                                                                                                                                                                                                                                                    </td><td>0.0072506723                                                                                                                                                                                                                                                                                                   </td><td> 1.474485                                                                                                                                                                                                                                                                                                      </td><td> 8                                                                                                                                                                                                                                                                                                             </td><td> 63                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000173787647411240331042768646534077965901504037252141338922282270454305363517476727541454150513743635150806181686607422593273462364061866514206974120838305060983629103399646569840159819263280350648115</td></tr>
	<tr><td>Hereditary spastic paraplegia - adult onset                                                                                                                                                                                                                                                                    </td><td>0.014858651                                                                                                                                                                                                                                                                                                    </td><td>0.0053592603                                                                                                                                                                                                                                                                                                   </td><td> 1.471197                                                                                                                                                                                                                                                                                                      </td><td>12                                                                                                                                                                                                                                                                                                             </td><td>113                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000004122513541809179652805296921940257526037831225314593596294620998451097374447662416144857025802237783597889596034956546578496201096816002365247299893403944171171711358629228479119730908179081163567047951801214367151260375976562</td></tr>
	<tr><td>Hereditary spastic paraplegia - childhood onset                                                                                                                                                                                                                                                                </td><td>0.014858651                                                                                                                                                                                                                                                                                                    </td><td>0.0053592603                                                                                                                                                                                                                                                                                                   </td><td> 1.471197                                                                                                                                                                                                                                                                                                      </td><td>12                                                                                                                                                                                                                                                                                                             </td><td>113                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000004122513541809179652805296921940257526037831225314593596294620998451097374447662416144857025802237783597889596034956546578496201096816002365247299893403944171171711358629228479119730908179081163567047951801214367151260375976562</td></tr>
	<tr><td>Palmoplantar keratoderma and erythrokeratodermas                                                                                                                                                                                                                                                               </td><td>0.009544006                                                                                                                                                                                                                                                                                                    </td><td>0.0037327661                                                                                                                                                                                                                                                                                                   </td><td> 1.354350                                                                                                                                                                                                                                                                                                      </td><td> 5                                                                                                                                                                                                                                                                                                             </td><td> 44                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000000000000000000000000000000000000032391299628367962548350574154539812733773096437159856083194859285340279846669897854329864494009088708280490839459416597479260446768023484764853492379188537597656250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Primary immunodeficiency                                                                                                                                                                                                                                                                                       </td><td>0.054613528                                                                                                                                                                                                                                                                                                    </td><td>0.0233774972                                                                                                                                                                                                                                                                                                   </td><td> 1.224138                                                                                                                                                                                                                                                                                                      </td><td>33                                                                                                                                                                                                                                                                                                             </td><td>239                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001387946</td></tr>
	<tr><td>Multiple Epiphyseal Dysplasia                                                                                                                                                                                                                                                                                  </td><td>0.005832456                                                                                                                                                                                                                                                                                                    </td><td>0.0025259415                                                                                                                                                                                                                                                                                                   </td><td> 1.207282                                                                                                                                                                                                                                                                                                      </td><td> 4                                                                                                                                                                                                                                                                                                             </td><td>  7                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000000000000000002009754753370975339684698117687548264891773096514798630111531914427400589864651514471205473455484025180339813232421875000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Cytopenia - NOT Fanconi anaemia                                                                                                                                                                                                                                                                                </td><td>0.006362683                                                                                                                                                                                                                                                                                                    </td><td>0.0028964414                                                                                                                                                                                                                                                                                                   </td><td> 1.135354                                                                                                                                                                                                                                                                                                      </td><td> 5                                                                                                                                                                                                                                                                                                             </td><td> 29                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000000000000000000001577243205255980114009179672857890487025970248716977211427286121752449453367786295649582983813274950080085545778274536132812500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Hereditary spastic paraplegia                                                                                                                                                                                                                                                                                  </td><td>0.016979544                                                                                                                                                                                                                                                                                                    </td><td>0.0077601070                                                                                                                                                                                                                                                                                                   </td><td> 1.129649                                                                                                                                                                                                                                                                                                      </td><td>14                                                                                                                                                                                                                                                                                                             </td><td>123                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000036038058677251167544260529652315537514389810806223681335436529418603849981060958773693013369559715639835696808152264222046827415722359169481947672682856598826598852357563872586834527915054223313027287501799284213</td></tr>
	<tr><td>Haematuria                                                                                                                                                                                                                                                                                                     </td><td>0.003182469                                                                                                                                                                                                                                                                                                    </td><td>0.0014559813                                                                                                                                                                                                                                                                                                   </td><td> 1.128155                                                                                                                                                                                                                                                                                                      </td><td> 5                                                                                                                                                                                                                                                                                                             </td><td> 50                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000676608359289451903475375679708776349041048428500921829709113808348774909973144531250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Nephrocalcinosis or nephrolithiasis                                                                                                                                                                                                                                                                            </td><td>0.007423122                                                                                                                                                                                                                                                                                                    </td><td>0.0034071421                                                                                                                                                                                                                                                                                                   </td><td> 1.123464                                                                                                                                                                                                                                                                                                      </td><td> 7                                                                                                                                                                                                                                                                                                             </td><td> 33                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000000000000000000000000002060940286195191438170006482250017661886891367808316183970955792662431847952032468439457838762975810942279508708452340215444564819335937500000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Xeroderma pigmentosum, Trichothiodystrophy or Cockayne syndrome                                                                                                                                                                                                                                                </td><td>0.009017185                                                                                                                                                                                                                                                                                                    </td><td>0.0042135425                                                                                                                                                                                                                                                                                                   </td><td> 1.097643                                                                                                                                                                                                                                                                                                      </td><td>10                                                                                                                                                                                                                                                                                                             </td><td> 78                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000000000000000000000000000000000002598585381139261369133051724813083606193893194943431392339441920181798567731362020471161078540714462817944744836974428625619992772044497542083263397216796875000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Proteinuric renal disease                                                                                                                                                                                                                                                                                      </td><td>0.012196827                                                                                                                                                                                                                                                                                                    </td><td>0.0058237265                                                                                                                                                                                                                                                                                                   </td><td> 1.066491                                                                                                                                                                                                                                                                                                      </td><td>14                                                                                                                                                                                                                                                                                                             </td><td>117                                                                                                                                                                                                                                                                                                            </td><td>0.0000000000000000000000000000000000000000000000000000000000000004264419922479812023761066260323947148959112539977797846227851329727493272064977370240516202367208130624668535913832779079023109759657360859677970211471182749619901475313099581399001181125640869140625000000000000000000000000000000000000000</td></tr>
	<tr><td>Pyruvate dehydrogenase (PDH) deficiency                                                                                                                                                                                                                                                                        </td><td>0.001591232                                                                                                                                                                                                                                                                                                    </td><td>0.0034383292                                                                                                                                                                                                                                                                                                   </td><td>-1.111563                                                                                                                                                                                                                                                                                                      </td><td> 3                                                                                                                                                                                                                                                                                                             </td><td> 26                                                                                                                                                                                                                                                                                                            </td><td>0.0000106870690963555931642147214799543064600584330037236213684082031250000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
	<tr><td>Gastrointestinal neuromuscular disorders                                                                                                                                                                                                                                                                       </td><td>0.001590669                                                                                                                                                                                                                                                                                                    </td><td>0.0040879575                                                                                                                                                                                                                                                                                                   </td><td>-1.361747                                                                                                                                                                                                                                                                                                      </td><td> 3                                                                                                                                                                                                                                                                                                             </td><td> 27                                                                                                                                                                                                                                                                                                            </td><td>0.0000099972138675092457192790373898816369546693749725818634033203125000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000</td></tr>
</tbody>
</table>




```R
nrow(selPanels)
```


24



```R
set.seed(42)
p <- ggplot(selPanels, aes(cum_AF_sum, cum_gnomad_AF_NFE_sum)) +
    geom_point() +
    # geom_text(aes(label = Symbol),hjust = 0, nudge_x = 0.0001)+ 
    geom_abline() + 
    #xlim(0,0.018) + ylim(0,0.018)+
    ggtitle("Cumulative allele frequencies of selected ClinVar variants in AR gene panels")+
    ylab("Per-panel cumulative allele frequencies in GnomAD NFE population")+
    xlab("Per-panel cumulative allele frequencies in Polish population")

p<- p +   geom_label_repel(
  aes(
    fill = "red", 
    label = Level4),
  fontface = 'bold.italic', 
    color = 'white',
  size = 2,
     segment.color = "black",
    #max.overlaps = 0,
    min.segment.length = 0,
  box.padding = unit(0.1, "lines"),
 point.padding = unit(0.1, "lines")
)
p +  theme(legend.position = "none") + theme(plot.title = element_text(size=11))
```

    Warning message:
    “ggrepel: 4 unlabeled data points (too many overlaps). Consider increasing max.overlaps”


    
![png](output_45_1.png)
    


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


    
![png](output_47_1.png)
    



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


    
![png](output_50_1.png)
    


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


    
![png](output_52_1.png)
    



    
![png](output_52_2.png)
    



    
![png](output_52_3.png)
    



    
![png](output_52_4.png)
    



    
![png](output_52_5.png)
    



    
![png](output_52_6.png)
    



    
![png](output_52_7.png)
    



    
![png](output_52_8.png)
    



    
![png](output_52_9.png)
    



    
![png](output_52_10.png)
    



    
![png](output_52_11.png)
    



    
![png](output_52_12.png)
    



    
![png](output_52_13.png)
    



    
![png](output_52_14.png)
    



    
![png](output_52_15.png)
    



    
![png](output_52_16.png)
    



    
![png](output_52_17.png)
    



    
![png](output_52_18.png)
    



    
![png](output_52_19.png)
    



    
![png](output_52_20.png)
    



    
![png](output_52_21.png)
    



    
![png](output_52_22.png)
    



    
![png](output_52_23.png)
    



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




    
![png](output_52_25.png)
    



    
![png](output_52_26.png)
    



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


    
![png](output_53_1.png)
    



    
![png](output_53_2.png)
    



    
![png](output_53_3.png)
    



    
![png](output_53_4.png)
    



    
![png](output_53_5.png)
    



    
![png](output_53_6.png)
    



    
![png](output_53_7.png)
    



    
![png](output_53_8.png)
    



    
![png](output_53_9.png)
    



    
![png](output_53_10.png)
    



    
![png](output_53_11.png)
    



    
![png](output_53_12.png)
    



    
![png](output_53_13.png)
    



    
![png](output_53_14.png)
    



    
![png](output_53_15.png)
    



    
![png](output_53_16.png)
    



    
![png](output_53_17.png)
    



    
![png](output_53_18.png)
    



    
![png](output_53_19.png)
    



    
![png](output_53_20.png)
    



    
![png](output_53_21.png)
    



    
![png](output_53_22.png)
    



    
![png](output_53_23.png)
    



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




    
![png](output_53_25.png)
    



    
![png](output_53_26.png)
    


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


    
![png](output_57_0.png)
    



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


    
![png](output_66_0.png)
    



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


    
![png](output_70_1.png)
    



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


    
![png](output_75_1.png)
    



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


    
![png](output_76_1.png)
    



```R

```


```R
tgp_vep_gene_table <- as.data.table(sort(table(tgp_vep_nonsyn_rare$Symbol), decreasing = TRUE))
ggplot(tgp_vep_gene_table[2:31], aes(reorder(V1, -N), N)) + 
    geom_col() +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    labs(title = "Total no. of variants per gene in TGP cohort", x = "Symbol", y = "Count")
```


    
![png](output_78_0.png)
    



```R
tgp_vep_gene_table_singleton <- as.data.table(sort(table(tgp_vep_nonsyn_rare[which(AF < 0.01 & AC == 1)]$Symbol), decreasing = TRUE))
```


```R
ggplot(tgp_vep_gene_table_singleton[2:31], aes(reorder(V1, -N), N)) + 
    geom_col() +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    labs(title = "Total no. of singleton variants per gene in TGP cohort", x = "Symbol", y = "Count")
```


    
![png](output_80_0.png)
    

