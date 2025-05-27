# Genetic Analyses

This branch contains the scripts and documentation for the genetic analyses performed as part of the study:  

---

## Overview of analysis steps

1. **Cluster analysis**
   - [identify_related_individuals.py](identify_related_individuals.py) : Identify groups of third degree relatives (refinedIBD sharing >= 0.125 and/or genealogical kinship >= 0.0625)
   - [compute_PC-AiR.R](compute_PC-AiR.R) : Perform a principal component analysis on the cleaned genotyped data, projecting the related subset onto the unrelated sample
   - [ancestry_clustering_PCA_UMAP_DBSCAN.R](ancestry_clustering_PCA_UMAP_DBSCAN.R) : Identify the meaningfull principal components, perform UMAP and identify clusters using DBSCAN
   - [ibd-ne.sh](ibd-ne.sh) : Infer effective population size in each cluster using IBD segments
  
3. ****


---

## Data availability

The data from the Eastern Quebec schizophrenia and bipolar disorder kindred study is available upon request from MM. The data from the CARTaGENE cohort is publicly available via an independent data access committee by the [CARTaGENE cohort](https://cartagene.qc.ca/en/researchers/access-request.html). 

---
