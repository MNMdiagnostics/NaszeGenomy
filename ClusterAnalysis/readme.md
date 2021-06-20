# Ancestry analysis #


The idea here is as follows:
1. obtain the 1000G reference files (in GRCh38, only SNPs, and autosomal chromosomes), normalize and left-align them with the reference genome
2. from the 1000G file, only use the variants also found in our specific cohort, and further prune them to MAF over 1% and in linkage equilibrium
3. from our cohort variant file, only select those pruned variants from step 2
4. merge the 1000G file and our cohort's file, with only those pruned variants
5. perform PCA on 1000G using this variant set, and project our cohort's genotype on the resulting PCs
6. train a random forest with 6 principal components on the 1000G dataset
7. predict the ancestry in your cohort using that trained random forest and your cohort's projection on 1000G's PCs

**The results are presented on the following plot:**
![Cluster Plot](https://github.com/MNMdiagnostics/NaszeGenomy/blob/main/ClusterAnalysis/ancestry.jpeg)

**One individual (ID=204_20349_20) is from AMR not EUR population**

# Fst analysis #
Fst statistics computed using the method introduced in Weir BS, Cockerham CC (1984) Estimating F-statistics for the analysis of population structure. The results are presented between each population available in the 1000G data. In the brackets is the mean of the weighted Fst statistics.

**The results are presented on the following plot:**
![Cluster Plot](https://github.com/MNMdiagnostics/NaszeGenomy/blob/main/ClusterAnalysis/FST2.jpeg)

# Cluster analysis #
Based on *3 491* genomes (*2 548* from *1000* Human Genomes and *943* from our project) we estimated the first *20* principal component for each of the individual (as presented in description of ancestry analysis). Then we took into account only patients from our cohort we detected possible outliers in the data. Now we can observe one huge cluster with few possible outliers:
![Cluster Plot](https://github.com/MNMdiagnostics/NaszeGenomy/blob/main/ClusterAnalysis/Int.jpeg) 

**The results are presented on the following plot:**
![Cluster Plot](https://github.com/MNMdiagnostics/NaszeGenomy/blob/main/ClusterAnalysis/ClusterPlot.jpeg) 
