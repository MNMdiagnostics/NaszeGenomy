# Ancestry analysis #


The idea here is as follows:
1. obtain the 1000G reference files (in GRCh38, only SNPs, and autosomal chromosomes), normalize and left-align them with the reference genome
2. from the 1000G file, only use the variants also found in our specific cohort, and further prune them to MAF over 1% and in linkage equilibrium
3. from our cohort variant file, only select those pruned variants from step 2
4. merge the 1000G file and our cohort's file, with only those pruned variants
5. perform PCA on 1000G using this variant set, and project our cohort's genotype on the resulting PCs
6. train a random forest with 6 principal components on the 1000G dataset
7. predict the ancestry in your cohort using that trained random forest and your cohort's projection on 1000G's PCs
8. I choosed also only European populations and plot again our cohort across other Eyropean populations (CEU, FIN, GBR, IBS, TSI). 

**Fife individuals (out from 1 076) do not belong to the EUR population. Four individuals are from AMR and one from EAS population**
| patient ID | Population |
| :---: | :---: |
| 204_20349_20 | AMR | 
| 214_21210_20 | AMR |
| 358_23976_20 | AMR |
| 907_8428_21 | AMR |
|902_5417_21 | EAS |

**The results are presented on the following plot:**
![Cluster Plot](https://github.com/MNMdiagnostics/NaszeGenomy/blob/main/ClusterAnalysis/ancestryPCAFULL.jpeg) 

**Version for Paula D.:**
![Cluster Plot](https://github.com/MNMdiagnostics/NaszeGenomy/blob/main/ClusterAnalysis/ancestryPCAFULLpopV2.jpeg) 
Na podstawie wspólnych markerów genetycznych dla populacji polskiej oraz danych dostepnych w ramach projektu 1000G zostały wyznaczone współczynniki charakteryzujące genetyczną strukturę sześciu populacji (5 populacji dostępnych w ramach projektu 1000G oraz populacja polska). Następnie dane dla populacji polskiej zostały naniesione na dane z projektu 1000G. Okazało się, że struktura populacji polskiej jest jednorodna i pokrywa się kontynentalnie z populacją europejską. Najbliżej populacji polskiej jest subpopulacja europejska CEU do której należy 1 063 (99% populacji polskiej) pacjentów.

I again trained a random forest with 6 principal components on the 1000G dataset and predicted the ancestry in our cohort, but now we checked if people from our came from populations: CEU, FIN, GBR, IBS, TSI or OTHER (AFR, AMR, SAS, EAS). I got the following results:

| Population | Number of people <br /> from our cohort |
| :---: | :---: |
| CEU | 1 063 | 
| FIN | 0 |
| GBR | 0 |
| IBS | 0 |
| TSI | 2 |
| OTHER | 11|

# Fst analysis #
Fst statistics computed using the method introduced in Weir BS, Cockerham CC (1984) Estimating F-statistics for the analysis of population structure. The results are presented between each population available in the 1000G data. In the brackets is the mean of the weighted Fst statistics.

**The results are presented on the following plot:**
![Cluster Plot](https://github.com/MNMdiagnostics/NaszeGenomy/blob/main/ClusterAnalysis/FST2v2.jpeg)

Fst statistics computed as an average of sliding window of 10 000 SNPs for all populations.

**The results are presented on the following plot:**
![Cluster Plot](https://github.com/MNMdiagnostics/NaszeGenomy/blob/main/ClusterAnalysis/FST222v2.jpeg)

Comparison of mean weighted Fst statistics within the European subpopulations:

| Population | CEU | FIN | GBR | IBS | TSI |
| :---: | :---: | :---: | :---: | :---: | :---: |
| wFst | 0.002 | 0.009 | 0.007 | 0.008 | 0.009 | 

# Cluster analysis #
Based on *3 624* genomes (*2 548* from *1000* Human Genomes and *1 076* from our project) we estimated the first *20* principal component for each of the individual (as presented in description of ancestry analysis). Then we took into account only patients from our cohort we detected possible outliers in the data. Now we can observe one huge cluster with few possible outliers:
![Cluster Plot](https://github.com/MNMdiagnostics/NaszeGenomy/blob/main/ClusterAnalysis/Int.jpeg) 

So we try to find outliers (marked in red) using two different methods:
1. Unsupervised Outlier Detection (OneClassSVM Python function with the following parameters: kernel='sigmoid',gamma='scale',nu=0.05;
2. Unsupervised Outlier Detection using Local Outlier Factor (LOF) (LocalOutlierFactor Python function with tthe following parameters: n_neighbors=20, algorithm='brute'.

**The results are presented on the following plot:**
![Cluster Plot](https://github.com/MNMdiagnostics/NaszeGenomy/blob/main/ClusterAnalysis/UOD.jpeg) 
![Cluster Plot](https://github.com/MNMdiagnostics/NaszeGenomy/blob/main/ClusterAnalysis/LOF.jpeg) 

Also I did DBSCAN clustering - Density-Based Spatial Clustering of Applications with Noise. Finds core samples of high density and expands clusters from them. This is very good method for data which contains clusters of similar density. It proves that we have only one cluster and few outliers.

**The results are presented on the following plot:**
![Cluster Plot](https://github.com/MNMdiagnostics/NaszeGenomy/blob/main/ClusterAnalysis/clusterPCA.jpeg) 

**Mean weighted Fst statistics between people from cluster (black dots) and outliers (8 red dots) in above plot is equal to 0.003.**

