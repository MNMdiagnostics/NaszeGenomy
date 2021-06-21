Runs of homozygosity on 943 individuals
================

| stat                |  min |  median |      mean |        max |
|:--------------------|-----:|--------:|----------:|-----------:|
| Length              | 83.0 | 54619.0 | 124596.53 | 63912719.0 |
| Number\_of\_markers |  2.0 |   100.0 |    166.18 |    45797.0 |
| Quality             |  0.6 |    51.1 |     49.86 |       98.7 |

ROHs quality histogram

![](roh_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Results

### All results below are ROHs filter for autosomes with Quality &gt; 25%

1.  Average sum of ROHs per genome

![](roh_files/figure-gfm/total_roh-1.png)<!-- -->

|     | mean\_SROH\_Mb  |
|:----|:----------------|
|     | Min. : 5.821    |
|     | 1st Qu.: 28.883 |
|     | Median : 33.855 |
|     | Mean : 33.207   |
|     | 3rd Qu.: 37.978 |
|     | Max. :131.589   |

| Range    | mean\_length |
|:---------|-------------:|
| 2Mb-5Mb  |     2.892791 |
| &gt;10Mb |    16.239240 |

2.  Sum of ROHs per chromosome and individual

<!-- -->

    ## `summarise()` has grouped output by 'sample_id', 'Chromosome'. You can override using the `.groups` argument.

![](roh_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

    ## `summarise()` has grouped output by 'sample_id'. You can override using the `.groups` argument.

    ## Adding missing grouping variables: `sample_id`

|     | sample\_id       | mean\_SROH\_Mb |
|:----|:-----------------|:---------------|
|     | Length:6575      | Min. : 2.000   |
|     | Class :character | 1st Qu.: 2.358 |
|     | Mode :character  | Median : 2.972 |
|     | NA               | Mean : 4.763   |
|     | NA               | 3rd Qu.: 3.900 |
|     | NA               | Max. :38.193   |

3.  Sum of ROH length per range

<!-- -->

    ## `summarise()` has grouped output by 'Range'. You can override using the `.groups` argument.

![](roh_files/figure-gfm/SROH-1.png)<!-- -->

|     | 2Mb-5Mb        | &gt;10Mb       |
|:----|:---------------|:---------------|
|     | Min. : 5.815   | Min. : 10.00   |
|     | 1st Qu.:15.537 | 1st Qu.: 15.89 |
|     | Median :18.727 | Median : 16.76 |
|     | Mean :19.026   | Mean : 17.77   |
|     | 3rd Qu.:21.975 | 3rd Qu.: 17.17 |
|     | Max. :42.199   | Max. :104.20   |

4.  Relationship between number of ROHs and total length of genome
    covered by them

![](roh_files/figure-gfm/cummulative-1.png)<!-- -->

5.  ROH genome coverage

![](roh_files/figure-gfm/genome_coverage-1.png)<!-- -->

![](roh_files/figure-gfm/genome_coverage_heatmap-1.png)<!-- -->

6.  % of ROhs per category in sample

![](roh_files/figure-gfm/roh_count-1.png)<!-- -->

7.  ROH length sum

![](roh_files/figure-gfm/roh_sum-1.png)<!-- -->

<!-- 7. Number of ROHs per sample -->
<!-- ```{r roh_count_per_sample, echo=FALSE}  -->
<!-- suppressMessages( -->
<!-- roh_avg <- roh %>% group_by(sample_id, Range) %>% -->
<!--   summarise(avg_length=mean(Length), n=n()) %>% -->
<!--   arrange(Range,n) -->
<!-- ) -->
<!-- roh_avg %>% -->
<!--   ggplot(aes(x=1:nrow(roh_avg),y=n)) + -->
<!--   geom_point(col='#48C095',alpha=0.5) + -->
<!--   facet_wrap(~Range, nrow = 3,scales = 'free') +  -->
<!--   ylab('Number of ROHs') + -->
<!--   xlab('Individual') + -->
<!--   theme_classic() + -->
<!--   theme(axis.text.x = element_blank()) -->
<!-- ``` -->
<!-- 8. Average ROHs per sample -->
<!-- ```{r average_roh_per_sample, echo=FALSE}  -->
<!-- roh_avg <- roh %>% group_by(sample_id, Range) %>% -->
<!--   summarise(avg_length=mean(Length), n=n()) %>% -->
<!--   arrange(Range,avg_length) -->
<!-- roh_avg %>% -->
<!--   ggplot(aes(x=1:nrow(roh_avg),y=avg_length)) + -->
<!--   geom_point(col='#48C095',alpha=0.5) + -->
<!--   facet_wrap(~Range, nrow = 3,scales = 'free') +  -->
<!--   ylab('Mean ROH length') + -->
<!--   xlab('Individual') + -->
<!--   theme_classic() + -->
<!--   theme(axis.text.x = element_blank()) -->
<!-- ``` -->

8.  Cummulative sum

<!-- -->

    ## `summarise()` has grouped output by 'sample_id', 'Range'. You can override using the `.groups` argument.

    ## `geom_smooth()` using method = 'loess' and formula 'y ~ x'

![](roh_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

9.  Cosanguinity in population

![](roh_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
