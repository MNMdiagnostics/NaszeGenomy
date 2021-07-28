Haplogroups
================

1.  Classification quality

<!-- -->

    ## Warning: Removed 2 rows containing non-finite values (stat_ydensity).

    ## Warning: Removed 2 rows containing non-finite values (stat_boxplot).

![](haplogroups_files/figure-gfm/quality_violin-1.tiff)<!-- -->

``` r
haplo_input <- haplo_input %>% filter(Quality > 0.80)

haplo_summary <- haplo_input %>%
                 group_by(clad) %>% 
                 summarise(n=n()) %>% arrange(n)

haplo_summary$n_perc <- (haplo_summary$n/sum(haplo_summary$n)) * 100
```

2.  % individuals per haplogroup

![](haplogroups_files/figure-gfm/percent_haplogroup-1.tiff)<!-- -->

4.  Subclades

<!-- -->

    ## `summarise()` has grouped output by 'subclad'. You can override using the `.groups` argument.

    ## `summarise()` has grouped output by 'clad'. You can override using the `.groups` argument.

    ## `summarise()` has grouped output by 'clad', 'subclad'. You can override using the `.groups` argument.

![](haplogroups_files/figure-gfm/percent_subclads-1.tiff)<!-- -->

<!-- ```{r subclads,echo=FALSE, eval=F} -->
<!-- # n <- paste('(',h_summary$n,')',sep = '') -->
<!-- # labels <- paste(round(h_summary$n_perc,2),'%',sep=' ',n) -->
<!-- #  -->
<!-- h_summary <- haplo_input  -->
<!-- h_summary$subclad <- substring(h_summary$subclad,1,3) -->
<!-- h_summary$subclad <- gsub('[a-zA-Z]$','',h_summary$subclad) -->
<!-- h_summary <- h_summary %>% group_by(Haplogroup,subclad) %>% summarise(n=n()) -->
<!-- h_summary$subclad <- ifelse(h_summary$n < 7, 'other',h_summary$subclad) -->
<!-- h_summary <- h_summary %>% group_by(Haplogroup,subclad) %>% summarise(n=sum(n)) -->
<!-- h_summary$n_perc <- (h_summary$n/sum(h_summary$n)) * 100 -->
<!-- top <- haplo_summary %>% top_n(n = 5,n) %>% select(Haplogroup) -->
<!-- top_h <- h_summary  %>% filter(Haplogroup %in% top$Haplogroup)  -->
<!-- top_h[1,2] <- 'H' -->
<!-- library(RColorBrewer) -->
<!-- nb.cols <- 25 -->
<!-- mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols) -->
<!-- top_h %>% -->
<!--   ggplot(aes(y=reorder(Haplogroup, -n),x=n,fill=subclad)) + -->
<!--   geom_bar(position="fill", stat="identity") + -->
<!--   scale_fill_manual(values = mycolors) + -->
<!--   theme_classic() + -->
<!--   ylab('Haplogroup') + -->
<!--   xlab('# individuals')  -->
<!-- ``` -->
