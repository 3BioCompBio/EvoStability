

# Supplementary Data: Exploring evolution to enhance mutational stability prediction

## Description

Here are the **Supplementary Data** for the publication [1] (2024) <https://www.biorxiv.org/content/10.1101/2024.05.28.596203v1>.

**Authors**: Pauline Hermans, Matsvei Tsishyn, Martin Schwersensky, Marianne Rooman, Fabrizio Pucci

## Content

### (1) Summary

A per protein summary for the data set `D` is given in `./summary.csv`. For each protein/protein domain, it includes information such as sequence length, depth of the multiple sequence alignment used, and Pearson and Spearman correlation with ΔΔG of the evolutionary scores and ΔΔG predictions.
For each protein, the relative data file, sequence files and structure file are referenced by its `protein_id` property.

### (2) Datasets

For each mutation of the dataset `D`, experimental ΔΔG values, RSA of the mutated residue, evolutionary scores (`LOR` and `WLOR`) are in given in `./datasets/`.

### (3) Sequences

Sequences of proteins from the dataset `D` are located in `./fasta/`.

### (4) Structures

3D structures of proteins from the dataset `D` are located in `./pdb/`.

### (5) Multiple sequence alignments

The multiple sequence alignments used in the publication and their relative sequences' weights are located in `./msa/`.

## Conventions and units

- All energy values are given in kcal/mol.
- We use the convention that destabilizing mutations have positive ΔΔG values.

## References

  [1] Hermans P., Tsishyn M., Schwersensky M, Rooman, M. & Pucci, F. (2024). Exploring evolution to enhance mutational stability prediction.
