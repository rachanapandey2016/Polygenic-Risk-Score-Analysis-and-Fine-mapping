
#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --mem=100g
#SBATCH --tmp=100g
#SBATCH --output=plink_clumping.out

set -euo pipefail

/home/gcd8141/public/aazaidi/1kg_genotypes/plink2.1 \
  --pfile /1kg_genotypes/all_hg38 \
  --clump /eur_cleaned.tsv \
  --clump-p1 5e-8 --clump-r2 0.1 --clump-kb 250 \
  --clump-snp-field RSID --clump-field P \
  --rm-dup exclude-all \
  --out /output_folder/Yengo_height_EUR_clumped


