

# Supplementary Data: Exploring evolution to uncover insights into protein mutational stability

## Description

Here are the **Supplementary Data** for the publication [1] (2024) <https://www.biorxiv.org/content/10.1101/2024.05.28.596203v1>.

**Authors**: Pauline Hermans, Matsvei Tsishyn, Martin Schwersensky, Marianne Rooman, Fabrizio Pucci

## Content

### (1) Summary

A per-protein summary for datasets `D` (derived from Mega [2]) and `L` (derived from S4038 [3]) are given in `./summary_D.csv` and `./summary_L.csv` respectively. 
For each protein/protein domain, it includes information such as sequence length, depth of the multiple sequence alignment used, and Pearson and Spearman correlation with ΔΔG of the evolutionary scores and ΔΔG predictions.

### (2) Datasets

For each mutation of the datasets `D` and `L`, experimental ΔΔG values, RSA of the mutated residue and evolutionary scores (`LOR` and `WLOR`) are in given in `./datasets/`.

### (3) Sequences

Sequences of proteins from dataset `D` and `L` are located in `./fasta/`.

### (4) Structures

3D structures of proteins from dataset `D` and `L` are located in `./pdb/`.

### (5) Multiple sequence alignments

The multiple sequence alignments used in the publication and their relative sequences' weights are located in `./msa/`.

### (6) Estimation of $p$-values for pairwise correlation comparisons

File `./p_values.ods` contains all pairwise p-values discussed in Supplementary Section 7 which compares performances measured in Sections 3.1, 3.3 and 3.5.

## Conventions and units

- All energy values are given in kcal/mol.
- We use the convention that destabilizing mutations have positive ΔΔG values.

## References

  [1] Hermans P., Tsishyn M., Schwersensky M, Rooman, M. & Pucci, F. (2024). Exploring evolution to uncover insights into protein mutational stability.  
  [2] Tsuboyama, Kotaro, et al. "Mega-scale experimental analysis of protein folding stability in biology and design." Nature 620.7973 (2023): 434-444.  
  [3] Zheng, Feifan, et al. "Assessing computational tools for predicting protein stability changes upon missense mutations using a new dataset." Protein Science 33.1 (2024): e4861.  
