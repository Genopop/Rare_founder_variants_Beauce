# Genealogical Analyses

This branch contains the scripts and documentation for the genealogical analyses performed as part of the study:  

---

## Overview of analysis steps

1. **Data preparation**
   - [1st_step_select_regional_group.py](1st_step_select_regional_group.py) : From the whole population pedigree dataset, select probands from the regions of interest
   - [2nd_step_completeness_filtering.py](2nd_step_completeness_filtering.py) : Filter regional subsets for completeness
   - [3rd_step_remove_related.py](3rd_step_remove_related.py) : Keep only one individual out of each sibling group
   - [4th_step_keep_sedentary_probands.py](4th_step_keep_sedentary_probands.py) : Only keep probands married in the same region as their parents and grandparents
   - [pedigree_description.py](pedigree_description.py) : Calculate the number of probands, ancestors and founders and the depth of the clean pedigrees
  
2. **Descriptive analyses**
   - [kinship_per_decades.py](kinship_per_decades.py) : Filter out individuals related to the first degree in each decade, then calculate the mean kinship coefficients
   - [kinship_confidence_intervals.py](kinship_confidence_intervals.py) : Using the output kinship matrix from the previous step, calculate the confidence interval. Each decade is computed separately for time efficiency
   - [inbreeding_per_decades.py](inbreeding_per_decades.py) : Calculate the mean inbreeding coefficients for each decade as well as the confidence intervals
   - [ancestors_diversity_ratio.py](ancestors_diversity_ratio.py) : Calulate the number of time each ancestor appear in the region
   - [ADR_confidence_intervals.py](ADR_confidence_intervals.py) : Calculate the mean ancestors' diversity ratio across all bootstrapped iterations and compute confidence intervals

---

## Required Libraries

The following packages are required:

### Python

- `pandas`
- `numpy`
- [`GeneaKit`]([https://github.com/GPhMorin/geneakit](https://github.com/Genopop/geneakit))

---

## Note

Genealogical data is not shared publicly due to privacy constraints, but is available upon request via an independent data access committee by [BALSAC](https://balsac.uqac.ca/acces-donnees/). The data necessary to reproduce Figures 2 and 3 is available in the data folder. It comprises the mean kinship, mean inbreeding and mean ADR with the 95% confidence intervals at each decade for the three regional groups (Beauce, SLSJ and Montreal).
