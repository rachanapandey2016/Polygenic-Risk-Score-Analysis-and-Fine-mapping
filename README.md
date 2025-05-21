# Polygenic-Risk-Score-Analysis-and-Fine-mapping  
# Polygenic Risk Score Analysis and Fine-Mapping of Height-associated Variants Across European and African Populations

This project was completed as **Project 4 for the GCD 8141: Computational Genomics** course. All analyses, visualizations, and interpretations are original work by Rachana Pandey.

## 📘 Project Overview and Methodology

Human height is a complex polygenic trait shaped by numerous genetic variants with small effects. Leveraging the GIANT consortium summary statistics (Yengo et al., 2022) and 1000 Genomes genotypes, this study investigates the genetic architecture of height between European (EUR) and African (AFR) ancestry populations. Publicly available GWAS summary statistics for height from the Yengo 2022 (GIANT consortium) and Genotype data from 1000 Genomes Project Phase 3 were used for this study. For the GWAS SNPs visualization, Manhattan plots and QQ plots were generated from the summary statistics to visualize genome-wide association signals and assess p-value inflation from GIANT consortium. Effect size differences for height between European and African population from GIANT consortium was visualized using histogram. Then for the PRS calculation using the European summary statistics from Giant, independent SNPs were selected via clumping using plink2.1, with a p-value threshold of 5×10−8 and an LD threshold of r2<0.1. 1000 Genomes Project Phase 3 genotypes were used for PRS computation and LD estimation, focusing on EUR and AFR populations. SNP effect sizes from the GIANT summary statistics were used to calculate PRS for each individual population. The PRS distribution differences were visualized using density plot. Similarly, for fine mapping a top SNP (rs11645785) with largest effect size difference between European and African population was selected, and a 500kb window was fine mapped using the susieR package to compute credible sets of likely causal variants. 

The key components of the project include:

- Genome-wide visualization of height associated GWAS signals (Manhattan and QQ plots)
- Comparison of height-related effect sizes across ancestries
- Construction and analysis of polygenic risk scores (PRS) across population of different ancestry
- Fine-mapping of top loci to detect credible causal variants

## 🧪 Major Findings

- **GWAS Visualization:** Manhattan and QQ plots revealed strong polygenic signals for height, with significant inflation at the tail end of p-values.
- **Ancestry Differences:** While most SNP effect sizes were similar across populations, a subset of top 10 SNPs showed large ancestry-specific effects.
- **PRS Transferability:** PRS scores were significantly higher in EUR compared to AFR populations when using EUR-derived summary stats—highlighting the challenge of cross-ancestry PRS transferability.
- **Fine-Mapping:** Fine-mapping around the SNP rs11645785 (which had the largest effect size difference) identified a credible set of causal variants, likely due to extended LD in EUR populations.

## 🧰 Tools and Data

- **Summary Stats:** GIANT Consortium (Yengo 2022)
- **Genotypes:** 1000 Genomes Phase 3
- **Software:** `plink2.1`, `susieR`, R packages including `ggplot2`, `data.table`

## 📁 Repository Structure  
PRS-Finemapping-Height-GCD 8141  
├── report/
│ ├── Project-4-GCD-Report.pdf # Final project report
│ └── Project4_GCD8814.Rmd # Annotated R Markdown file with full analysis pipeline
│
├── results/
│ ├── manhattan_plot.png # manhattan plot for height GWAS across all population in giant consortium
│ ├── qq_plot.png # QQ plot for for height GWAS and p-value inflation across all population in giant consortium
│ ├── effect_size_diff.png # distribution of effect size difernce between european and african in giant consortium
│ ├── PRS_distribution.png # PRS density distribution across EUR and AFR
│ ├── Top_10_SNPs #top 10 SNps with largest effect size difference in european and african in giant 2022
│ ├──  PRS_distribution #PRS distribution across 1k genome european and african population
│ └── fine_mapping_rs11645785.png # Regional finemapping of top SNP rs11645785 in giant consortium

## 📜 License
This project is for academic use and coursework submission only. Do not reuse without permission.

## 🙋‍♀️ Author
**Rachana Pandey**  
PhD Student, Translational Bioinformatics  
University of Minnesota Twin Cities
