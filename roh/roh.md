Runs of homozygosity on 1082 individuals
================

| stat                |  min |  median |      mean |        max |
|:--------------------|-----:|--------:|----------:|-----------:|
| Length              | 83.0 | 54628.0 | 124447.48 | 63912719.0 |
| Number\_of\_markers |  2.0 |   100.0 |    166.01 |    45797.0 |
| Quality             |  0.6 |    51.1 |     49.86 |       98.7 |

ROHs quality histogram

![](roh_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Results

### All results below are ROHs filter for autosomes with Quality &gt; 25% andNumber of markers &gt;= 50

Number of ROHs with length in specific ranges

![](roh_files/figure-gfm/roh-1.png)<!-- -->

% of ROhs per category in sample

    ## `summarise()` has grouped output by 'Range'. You can override using the `.groups` argument.

![](roh_files/figure-gfm/roh_count-1.png)<!-- -->

ROH length sum

    ## `summarise()` has grouped output by 'Range'. You can override using the `.groups` argument.

![](roh_files/figure-gfm/roh_sum-1.png)<!-- -->

Sum of ROH length per range

    ## `summarise()` has grouped output by 'Range'. You can override using the `.groups` argument.

![](roh_files/figure-gfm/SROH-1.png)<!-- -->

Relationship between number of ROHs and total length of genome covered
by them

![](roh_files/figure-gfm/cummulative-1.png)<!-- -->

Number of ROHs per sample

    ## `summarise()` has grouped output by 'sample_id'. You can override using the `.groups` argument.

![](roh_files/figure-gfm/roh_count_per_sample-1.png)<!-- -->

Average ROHs per sample

    ## `summarise()` has grouped output by 'sample_id'. You can override using the `.groups` argument.

![](roh_files/figure-gfm/average_roh_per_sample-1.png)<!-- -->

Cosanguinity in population

![](roh_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
