# Polygenic-Risk-Score-Analysis-and-Fine-mapping  
# Polygenic Risk Score Analysis and Fine-Mapping of Height-associated Variants Across European and African Populations

This project was completed as **Project 4 for the GCD 8814: Computational Genomics** course. All analyses, visualizations, and interpretations are original work by Rachana Pandey.

## 📘 Project Overview

Human height is a complex polygenic trait shaped by numerous genetic variants with small effects. Leveraging the GIANT consortium summary statistics (Yengo et al., 2022) and 1000 Genomes genotypes, this study investigates the genetic architecture of height between European (EUR) and African (AFR) ancestry populations.

The key components of the project include:

- Genome-wide visualization of GWAS signals (Manhattan and QQ plots)
- Comparison of height-related effect sizes across ancestries
- Construction and analysis of polygenic risk scores (PRS)
- Fine-mapping of top loci to detect credible causal variants

## 🧪 Major Findings

- **GWAS Visualization:** Manhattan and QQ plots revealed strong polygenic signals for height, with significant inflation at the tail end of p-values.
- **Ancestry Differences:** While most SNP effect sizes were similar across populations, a subset showed large ancestry-specific effects.
- **PRS Transferability:** PRS scores were significantly higher in EUR compared to AFR populations when using EUR-derived summary stats—highlighting the challenge of cross-ancestry PRS transferability.
- **Fine-Mapping:** Fine-mapping around the SNP rs11645785 (which had the largest effect size difference) identified a credible set of causal variants, likely due to extended LD in EUR populations.

## 🧰 Tools and Data

- **Summary Stats:** GIANT Consortium (Yengo 2022)
- **Genotypes:** 1000 Genomes Phase 3
- **Software:** `plink2.1`, `susieR`, R packages including `ggplot2`, `data.table`

## 📁 Repository Structure


## 📜 License
This project is for academic use and coursework submission only. Do not reuse without permission.

## 🙋‍♀️ Author
**Rachana Pandey**  
PhD Student, Translational Bioinformatics  
University of Minnesota Twin Cities
