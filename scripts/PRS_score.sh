
#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --mem=100g
#SBATCH --tmp=100g
#SBATCH --output=PRS_score.out

set -euo pipefail

/home/gcd8141/public/aazaidi/1kg_genotypes/plink2.1 \
  --pfile /scratch.global/pande250/GCD-8141/Project-4/all_hg38 \
  --score /output_folder/GIANT_HEIGHT_clumped.effects 1 2 3 header \
  --out /output_folder/GIANT_HEIGHT_clumped.prs