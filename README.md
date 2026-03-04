# Polygenic-Risk-Score-Analysis-and-Fine-mapping  
# Polygenic Risk Score Analysis and Fine-Mapping of Height-associated Variants Across European and African Populations  

## 📘 Project Overview  
Polygenic risk scores are increasingly used in clinical research and population health studies; however, their predictive performance varies across ancestral groups. PRS performance often decreases when applied across ancestries due to differences in genetic architecture, allele frequencies, and linkage disequilibrium patterns This project demonstrates a reproducible workflow for integrating GWAS summary statistics with population genotype data to evaluate polygenic risk score transferability and perform locus-level fine mapping across diverse ancestral populations to imporve the robustness and interpretability. 
The key analysis of this project include:
- Genome-wide visualization of height associated GWAS signals (Manhattan and QQ plots)
- Comparison of height-related effect sizes across ancestries
- Construction and analysis of polygenic risk scores (PRS) across population of different ancestry
- Evaluation of PRS transferability across ancestries
- Fine-mapping of top loci to detect credible causal variants
## Background  
Human height is a complex polygenic trait shaped by numerous genetic variants with small effects. Leveraging the GIANT consortium summary statistics (Yengo et al., 2022) and 1000 Genomes genotypes, this study investigates the genetic architecture of height between European (EUR) and African (AFR) ancestry populations. Publicly available GWAS summary statistics for height from the Yengo 2022 (GIANT consortium) and Genotype data from 1000 Genomes Project Phase 3 were used for this study. Here, height is used for demonstration, the pipeline can be adapted to any trait with available GWAS summary statistics: Type 2 Diabetes, cancers, coronary artery disease. To do so:  
1. Replace `height_sumstats.txt.gz` with your summary statistics
2. Ensure the same column formatting (CHR, BP, SNP, A1, A2, BETA/OR, SE, P)
3. Run scripts in the order shown under `scripts/`
4. Evaluate PRS performance in your validation cohorts
All analyses, visualizations, and interpretations are original work by Rachana Pandey.

 ## Methods

For the GWAS SNPs visualization, Manhattan plots and QQ plots were generated from the summary statistics to visualize genome-wide association signals and assess p-value inflation from GIANT consortium. Effect size differences for height between European and African population from GIANT consortium was visualized using histogram. Then for the PRS calculation using the European summary statistics from Giant, independent SNPs were selected via clumping using plink2.1, with a p-value threshold of 5×10−8 and an LD threshold of r2<0.1. 1000 Genomes Project Phase 3 genotypes were used for PRS computation and LD estimation, focusing on EUR and AFR populations. SNP effect sizes from the GIANT summary statistics were used to calculate PRS for each individual population. The PRS distribution differences were visualized using density plot. Similarly, for fine mapping a top SNP (rs11645785) with largest effect size difference between European and African population was selected, and a 500kb window was fine mapped using the susieR package to compute credible sets of likely causal variants. 


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
```
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
```
## References  
1. Yengo, L., Vedantam, S., Marouli, E., Sidorenko, J., Bartell, E., Sakaue, S., ... & Lee, J. Y. (2022). A saturated map of common genetic variants associated with human height. Nature, 610(7933), 704-712.
2. Wang, G., Sarkar, A., Carbonetto, P., & Stephens, M. (2020). A simple new approach to variable selection in regression, with application to genetic fine mapping. Journal of the Royal Statistical Society Series B: Statistical Methodology, 82(5), 1273-1300.
3. Chang, C. C., Chow, C. C., Tellier, L. C., Vattikuti, S., Purcell, S. M., & Lee, J. J. (2015). Second-generation PLINK: rising to the challenge of larger and richer datasets. Gigascience, 4(1), s13742-015.


## 🙋‍♀️ Author
**Rachana Pandey**  
PhD Student, Translational Bioinformatics  
University of Minnesota Twin Cities
