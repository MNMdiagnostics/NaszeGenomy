Variants in disease causing genes ewline Results for 943 unrelated
individuals
================

### Samples count

    ## Warning: `funs()` was deprecated in dplyr 0.8.0.
    ## Please use a list of either functions or lambdas: 
    ## 
    ##   # Simple named list: 
    ##   list(mean = mean, median = median)
    ## 
    ##   # Auto named with `tibble::lst()`: 
    ##   tibble::lst(mean, median)
    ## 
    ##   # Using lambdas
    ##   list(~ mean(., trim = .2), ~ median(., na.rm = TRUE))

| Variant    | min\_count | mean\_count | max\_count |
|:-----------|-----------:|------------:|-----------:|
| Indels     |     751498 |      768905 |     781426 |
| Singletons |        436 |       16473 |      82329 |
| SNP        |    3637424 |     3715552 |    3776871 |

## Cummulative allele frequency

<!-- ```{r af_hist,echo=FALSE} -->
<!-- af <- read.table('../input/multisample_20210519.dv.bcfnorm.filtered.ACgt0.AF_list.tsv') -->
<!-- colnames(af) <- c('AF','id', 'allele_frequency', 'SNP', 'number_of_transitions', 'number_of transversions', 'indel', 'repeat-consistent','repeat-inconsistent', 'not_applicable') -->
<!-- af$allele_frequency[1] <- 0.000530223 -->
<!-- af_plot <- af %>% select(allele_frequency,SNP,indel)  -->
<!-- af_plot$SNP <- cumsum(af_plot$SNP)/1e+6 -->
<!-- af_plot$indel <- cumsum(af_plot$indel)/1e+6 -->
<!-- type.colors <- c(SNP = "#27384A", indel ="#48C095") -->
<!-- af_plot %>% -->
<!--   pivot_longer(-allele_frequency,names_to = 'Variant class', -->
<!--                values_to = 'Cummulative number of variants (millions)') %>% -->
<!--   ggplot(aes(allele_frequency,`Cummulative number of variants (millions)`, -->
<!--              col=`Variant class`)) + -->
<!--   geom_line() + -->
<!--   xlab('Allele frequency') + -->
<!--   scale_color_manual(values = type.colors) + -->
<!--   scale_y_continuous(breaks = seq(0,30,2)) + -->
<!--   theme_classic() -->
<!-- ``` -->

![](variants_af_files/figure-gfm/af_hist_pct-1.jpeg)<!-- -->

    ## `summarise()` has grouped output by 'svtype'. You can override using the `.groups` argument.
    ## `summarise()` has grouped output by 'svtype'. You can override using the `.groups` argument.

![](variants_af_files/figure-gfm/sv.af.hist-1.jpeg)<!-- -->

## ACMG

![](variants_af_files/figure-gfm/ACMG-1.jpeg)<!-- -->

### ClinVar variants

![](variants_af_files/figure-gfm/clinvar-1.jpeg)<!-- -->

| stars |   n |
|:------|----:|
| 1     | 325 |
| 2     | 293 |
| 3     |  22 |
| 4     |   1 |
| NA    | 122 |

### Putative variants

### % IMPACT variants

![](variants_af_files/figure-gfm/unnamed-chunk-3-1.jpeg)<!-- -->

<!-- ### 7. Variants per functional category -->
<!-- ```{r, echo=FALSE} -->
<!-- type.colors <- c(Exonic = "#27384A", Intronic ="#48C095", Noncoding = "#B6B6B6") -->
<!-- # consequence <- read.table('../input/diseases/consequence_ready.tsv',header=T) -->
<!-- # -->
<!-- # consequence$group <- factor(consequence$group) -->
<!-- # consequence$group <- ordered(consequence$group, levels = c("<0.1%", "0.1-0.5%", ">0.5%")) -->
<!-- # -->
<!-- # write.table(consequence,'consequence_summary.tsv',quote = F,col.names = T,sep='\t',row.names = F) -->
<!-- # -->
<!-- # cons_list <- consequence %>% group_by(Consequence) %>% -->
<!-- #   summarise(n = sum(n)) -->
<!-- #   write.table(cons_list,'consequence_list.tsv', quote = F,col.names = T,sep='\t',row.names = F) -->
<!-- # -->
<!-- #  consequence %>% -->
<!-- #   ggplot(aes(fill=type,y=n,x=group)) + -->
<!-- #   geom_bar(position="fill", stat="identity") + -->
<!-- #   theme_classic() + -->
<!-- #   scale_fill_manual(values = type.colors) + -->
<!-- #   xlab('Allele frequencies') + -->
<!-- #   ylab('% of variants') -->
<!-- # -->
<!-- # glimpse(consequence) -->
<!-- ``` -->

## Number of variants per impact

| VARIANT\_CLASS | AF       | HIGH |   LOW | MODERATE | MODIFIER |
|:---------------|:---------|-----:|------:|---------:|---------:|
| deletion       | &gt;0.5% |  500 |  1090 |      680 |  1697005 |
| insertion      | &gt;0.5% |  327 |  1280 |      630 |  1934758 |
| SNV            | &gt;0.5% | 1446 | 45903 |    39096 | 15294308 |
| deletion       | 0.1-0.5% |  813 |   582 |      859 |   885811 |
| insertion      | 0.1-0.5% |  407 |   648 |      573 |  1010671 |
| SNV            | 0.1-0.5% | 1723 | 32131 |    41239 |  9509163 |
| deletion       | &lt;0.1% | 2662 |   978 |     1731 |  1423061 |
| insertion      | &lt;0.1% | 1252 |   800 |     1015 |  1102882 |
| SNV            | &lt;0.1% | 5109 | 69962 |   103413 | 19741686 |

### Variants per consequence

    ## `summarise()` has grouped output by 'Coding_var_category'. You can override using the `.groups` argument.

![](variants_af_files/figure-gfm/unnamed-chunk-4-1.jpeg)<!-- -->

    ## `summarise()` has grouped output by 'noncoding.var_category'. You can override using the `.groups` argument.

![](variants_af_files/figure-gfm/unnamed-chunk-4-2.jpeg)<!-- -->

## NBS: chr8\_89971213\_ATTTGT\_A

![](variants_af_files/figure-gfm/NBS-1.jpeg)<!-- -->

## Cystic fybrosis: chr7\_117559590\_ATCT\_A

![](variants_af_files/figure-gfm/Mucoviscidosis-1.jpeg)<!-- -->

## CFTR deletions

![](variants_af_files/figure-gfm/CFTR-1.jpeg)<!-- -->

| Uploaded\_variation       |    PL\_AF | gnomAD\_AF | gnomAD\_AFR\_AF | gnomAD\_AMR\_AF | gnomAD\_ASJ\_AF | gnomAD\_EAS\_AF | gnomAD\_FIN\_AF | gnomAD\_NFE\_AF | gnomAD\_OTH\_AF | gnomAD\_SAS\_AF |
|:--------------------------|----------:|-----------:|----------------:|----------------:|----------------:|----------------:|----------------:|----------------:|----------------:|----------------:|
| chr7\_117536514\_AGATT\_A | 0.1935310 |  0.2696000 |       0.2159000 |        0.384300 |        0.170500 |       0.4153000 |        0.300100 |       0.1875000 |        0.245500 |       0.3452000 |
| chr7\_117548606\_ATG\_A   | 0.2500000 |  0.2222000 |       0.4383000 |        0.144200 |        0.159400 |       0.0236600 |        0.260300 |       0.2703000 |        0.222900 |       0.1627000 |
| chr7\_117548606\_ATGTG\_A | 0.0005519 |  0.0013450 |       0.0056910 |        0.001099 |        0.000312 |       0.0000694 |        0.001146 |       0.0014410 |        0.001383 |       0.0003595 |
| chr7\_117548628\_GTT\_G   | 0.0282818 |  0.0223100 |       0.0669800 |        0.011230 |        0.017300 |       0.0014670 |        0.028040 |       0.0283400 |        0.022420 |       0.0063610 |
| chr7\_117548834\_TG\_T    | 0.0005308 |  0.0000085 |       0.0000000 |        0.000000 |        0.000000 |       0.0000000 |        0.000000 |       0.0000191 |        0.000000 |       0.0000000 |
| chr7\_117559590\_ATCT\_A  | 0.0095440 |  0.0070680 |       0.0029540 |        0.003673 |        0.005557 |       0.0000000 |        0.002222 |       0.0122700 |        0.006856 |       0.0019930 |
| chr7\_117666871\_CT\_C    | 0.0005302 |  0.0009436 |       0.0000617 |        0.000116 |        0.019400 |       0.0000000 |        0.000000 |       0.0002149 |        0.001802 |       0.0000000 |
