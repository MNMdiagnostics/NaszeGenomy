Runs of homozygosity on 943 individuals
================

## Runs of homozygozity (ROH)

#### ROHs quality histogram

``` r
# Quality filtering: Quality>25 & markers > 50
roh <- roh %>% 
  filter(Quality > 25 & Number_of_markers >50)
```

#### ROH stats after quality filtering

| stat                |   min |  median |      mean |        max |
|:--------------------|------:|--------:|----------:|-----------:|
| Length              | 609.0 | 72459.0 | 123018.00 | 29363916.0 |
| Number\_of\_markers |  51.0 |   132.0 |    193.66 |    33173.0 |
| Quality             |  25.1 |    60.4 |     59.40 |       98.7 |

Average ROH size per genome (additionally filtered: &gt;100kb)

| stat      |       average |
|:----------|--------------:|
| count     |       1124.41 |
| maxLen    |   14239847.30 |
| meanLen   |     252746.68 |
| medianLen |     164071.12 |
| minLen    |     100073.33 |
| sumLen    | 284174602\.38 |

## Results

1.  Cummulative sum

![](roh_files/figure-gfm/cummulative_sum-1.jpeg)<!-- -->

Fig. Relation between size of runs of homozygozity and cumulative length
of RoHs in an average individual.

![](roh_files/figure-gfm/lengths_per_chromosome-1.jpeg)<!-- -->

Fig. Distribution of mean length of runs of homozygozity per chromosome

![](roh_files/figure-gfm/numbers_per_chromosome-1.jpeg)<!-- -->

Fig. Distribution of counts of runs of homozygozity per chromosome

2.  Coverage heatmap

![](roh_files/figure-gfm/genome_coverage_heatmap-1.jpeg)<!-- -->

Fig. Chromosomal location of run of homozygozity hotspots - regions
frequently covered by runs of homozygozity in this cohort.

### All results below are ROHs filter for autosomes with Quality &gt; 25%, Number of markers &gt; 50 and ROH length &gt; 2 Mb

2.  Average sum of ROHs per genome

![](roh_files/figure-gfm/total_roh-1.jpeg)<!-- -->

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

3.  Sum of ROH length per range

![](roh_files/figure-gfm/SROH-1.jpeg)<!-- -->

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

![](roh_files/figure-gfm/SROH_corr-1.jpeg)<!-- -->
