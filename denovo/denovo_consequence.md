``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(knitr)
t<-read.table('denovo_consequences_per_sample.tsv', stringsAsFactors = T)
colnames(t)<-c('sample', 'conseq','impact','symbol')
t['coding'] = t$impact!='MODIFIER'
```

``` r
kable(t %>% count(impact))
```

| impact   |      n|
|:---------|------:|
| HIGH     |      8|
| LOW      |     72|
| MODERATE |    113|
| MODIFIER |  31653|

``` r
kable(t %>% count(conseq))
```

| conseq                                                                   |      n|
|:-------------------------------------------------------------------------|------:|
| 3\_prime\_UTR\_variant                                                   |    177|
| 3\_prime\_UTR\_variant,NMD\_transcript\_variant                          |      3|
| 5\_prime\_UTR\_variant                                                   |     32|
| downstream\_gene\_variant                                                |   2131|
| intergenic\_variant                                                      |  11185|
| intron\_variant                                                          |  10222|
| intron\_variant,NMD\_transcript\_variant                                 |    152|
| intron\_variant,non\_coding\_transcript\_variant                         |   4353|
| missense\_variant                                                        |    112|
| missense\_variant,splice\_region\_variant                                |      1|
| non\_coding\_transcript\_exon\_variant                                   |    251|
| regulatory\_region\_variant                                              |    874|
| splice\_acceptor\_variant                                                |      1|
| splice\_acceptor\_variant,non\_coding\_transcript\_variant               |      1|
| splice\_donor\_variant                                                   |      1|
| splice\_donor\_variant,non\_coding\_transcript\_variant                  |      1|
| splice\_region\_variant,5\_prime\_UTR\_variant                           |      1|
| splice\_region\_variant,intron\_variant                                  |     12|
| splice\_region\_variant,intron\_variant,non\_coding\_transcript\_variant |      5|
| splice\_region\_variant,synonymous\_variant                              |      1|
| stop\_gained                                                             |      4|
| synonymous\_variant                                                      |     53|
| TF\_binding\_site\_variant                                               |     74|
| upstream\_gene\_variant                                                  |   2199|

Number of rare coding de-novo per-sample

``` r
t[t$coding,] %>% count(sample, .drop=F) %>% summarize(mean=mean(n), min=min(n), max=max(n))
```

    ##       mean min max
    ## 1 2.097826   0   8

``` r
hist((t[t$coding,] %>% count(sample, .drop=F))$n, 
     main = '# rare de-novo exonic variants per sample',
     xlab = '# de-novo', ylab='# samples')
```

![](denovo_consequence_files/figure-markdown_github/coding.per.sample.hist-1.png)

``` r
kable(t %>% count(sample,impact, .drop=F) 
        %>% group_by(impact) 
        %>% summarize(mean=mean(n), min=min(n), max=max(n)))
```

| impact   |         mean|  min|  max|
|:---------|------------:|----:|----:|
| HIGH     |    0.0869565|    0|    1|
| LOW      |    0.7826087|    0|    3|
| MODERATE |    1.2282609|    0|    6|
| MODIFIER |  344.0543478|  280|  456|
