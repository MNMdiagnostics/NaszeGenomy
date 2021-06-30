
## 1. Read QC

#### 1.1 Depth statistics


| stat                  |   min | median |  mean |   max |
|:----------------------|------:|-------:|------:|------:|
| average\_depth        | 29.09 |  35.75 | 35.72 | 45.75 |
| percentage\_above\_10 | 91.40 |  91.98 | 91.92 | 92.34 |
| percentage\_above\_20 | 85.61 |  89.27 | 89.37 | 91.27 |
| percentage\_above\_30 | 53.66 |  79.59 | 78.41 | 90.00 |

![Mean coverage](https://github.com/MNMdiagnostics/NaszeGenomy/raw/main/qc/qc_files/figure-gfm/average_depth-1.png) 


## 2. Variants


## 3. Population analysis

An average Fst statistics over each sliding window built from 1,000 neighbouring SNPs is presented on Figure 3.1. This results prove that samples from our cohort are 
the most similar to the European population and the least to the African and East Asian populations.

![Cluster Plot](https://github.com/MNMdiagnostics/NaszeGenomy/blob/main/ClusterAnalysis/FST222.jpeg) 
**Figure 3.1. The average Fst statistics calculated for the sliding windows containing 1,000 SNP.**

On the Figure 3.2. we can observe only one huge cluster obtained using Density-Based Spatial Method (935 patients marked on black) and a few outliers (8 patients marked on red). 
The average Fst between cluster and outliers is equal to 0.012. 

**The results are presented on the following plot:**
![Cluster Plot](https://github.com/MNMdiagnostics/NaszeGenomy/blob/main/ClusterAnalysis/clusterPCA.jpeg) 
**Figure 3.2. Cluster analysis based on Density-Based Spatial method.**

