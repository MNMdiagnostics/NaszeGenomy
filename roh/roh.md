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

### All results below are ROHs filter for autosomes with Quality &gt; 25% and Number of markers &gt;= 50

1.  Average sum of ROHs

![](roh_files/figure-gfm/total_roh-1.png)<!-- -->

|     | mean\_SROH\_Mb  |
|:----|:----------------|
|     | Min. : 5.821    |
|     | 1st Qu.: 29.593 |
|     | Median : 34.647 |
|     | Mean : 34.556   |
|     | 3rd Qu.: 39.401 |
|     | Max. :144.543   |

2.  Sum of ROH length per range

<!-- -->

    ## `summarise()` has grouped output by 'Range'. You can override using the `.groups` argument.

![](roh_files/figure-gfm/SROH-1.png)<!-- -->

|     | &gt;10Mb       | 2-5Mb          |
|:----|:---------------|:---------------|
|     | Min. : 10.00   | Min. : 5.815   |
|     | 1st Qu.: 15.89 | 1st Qu.:15.537 |
|     | Median : 16.76 | Median :18.727 |
|     | Mean : 17.77   | Mean :19.026   |
|     | 3rd Qu.: 17.17 | 3rd Qu.:21.975 |
|     | Max. :104.20   | Max. :42.199   |

3.  Relationship between number of ROHs and total length of genome
    covered by them

![](roh_files/figure-gfm/cummulative-1.png)<!-- -->

4.  ROH genome coverage with runs &gt; 2Mb

![](roh_files/figure-gfm/genome_coverage-1.png)<!-- -->

5.  % of ROhs per category in sample

![](roh_files/figure-gfm/roh_count-1.png)<!-- -->

6.  ROH length sum

![](roh_files/figure-gfm/roh_sum-1.png)<!-- -->

7.  Number of ROHs per sample

![](roh_files/figure-gfm/roh_count_per_sample-1.png)<!-- -->

8.  Average ROHs per sample

<!-- -->

    ## `summarise()` has grouped output by 'sample_id'. You can override using the `.groups` argument.

![](roh_files/figure-gfm/average_roh_per_sample-1.png)<!-- -->

9.  Cosanguinity in population

![](roh_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
