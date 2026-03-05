library(dplyr)
library(susieR)
library(data.table)
library(ggplot2)

dir.create("visualization", showWarnings = FALSE, recursive = TRUE)
dir.create("output_folder", showWarnings = FALSE, recursive = TRUE)

# INPUTS (edit paths if your environment differs)
PLINK <- "/home/gcd8141/public/aazaidi/1kg_genotypes/plink2.1"
PFILE_1KG <- "/home/gcd8141/public/aazaidi/1kg_genotypes/all_hg38"   # 1000G hg38 pfile prefix
PSAM_1KG  <- "/home/gcd8141/public/aazaidi/1kg_genotypes/all_hg38.psam"

# Using our cleaned EUR GWAS file if you created it (recommended).
# Otherwise point this to the original EUR GWAS file.
EUR_GWAS <- "./eur_cleaned.tsv"  # or "./GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR"

# Output prefixes
OUT_PREFIX_PGEN <- "./output_folder/GIANT_HEIGHT_rs11645785_region"
OUT_PREFIX_LD   <- "./output_folder/GIANT_HEIGHT_locusrs11645785_maf_0.01"
EUR_IID_FILE    <- "./output_folder/1kg_hg38.eur.iid"
REGION_SNPS_FILE <- "./output_folder/rs11645785_locus.snps"

# 1) Load EUR GWAS and define fine-mapping region around rs11645785
eur <- fread(EUR_GWAS)

# Ensure expected columns exist (minimal checks)
required_cols <- c("CHR", "POS", "RSID", "P", "BETA", "SE", "EFFECT_ALLELE")
missing_cols <- setdiff(required_cols, colnames(eur))
if (length(missing_cols) > 0) {
  stop(paste("EUR GWAS is missing required columns:", paste(missing_cols, collapse = ", ")))
}

# We will use rs11645785 at chr16:47164632 (hg38) +/- 500kb
top_snp <- "rs11645785"
chr <- 16
pos <- 47164632

region_start <- pos - 500000
region_end   <- pos + 500000

# Subset the GWAS summary statistics to the locus
locus <- eur[CHR == chr & POS > region_start & POS < region_end]
locus[, P := as.numeric(P)]

# Rename RSID -> ID for consistency with PLINK pvar later
setnames(locus, "RSID", "ID")

# 2) Plot regional association signal

regional_plot <- ggplot(locus, aes(POS, -log10(P))) +
  geom_point() +
  theme_classic() +
  labs(
    x = "Position (hg38)",
    y = bquote(-Log[10]~"p-value"),
    title = paste0("Regional association signal around ", top_snp, " (chr", chr, ")")
  ) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

ggsave(
  filename = "./visualization/regional_assoc_rs11645785.pdf",
  plot = regional_plot,
  width = 10, height = 6
)


# 3) Writing SNP list in this region for PLINK extraction / LD calculation

fwrite(locus[, .(ID)], REGION_SNPS_FILE, col.names = FALSE)


# 4) Create EUR IID keep file (EUR individuals only)

system(paste0(
  "awk -F'\\t' 'NR==1 {next} $5==\"EUR\" {print $1, $1}' ",
  PSAM_1KG,
  " > ",
  EUR_IID_FILE
))


# 5) Extract region SNPs and EUR individuals into a new pgen

system(paste(
  PLINK,
  "--pfile", PFILE_1KG,
  "--extract", REGION_SNPS_FILE,
  "--keep", EUR_IID_FILE,
  "--allow-extra-chr",
  "--make-pgen",
  "--out", OUT_PREFIX_PGEN
))


# 6) Compute LD matrix (unphased vcor) with MAF filter

system(paste(
  PLINK,
  "--pfile", OUT_PREFIX_PGEN,
  "--maf 0.01",
  "--r-unphased square zs ref-based",
  "--out", OUT_PREFIX_LD
))


# 7) Read LD variant order + pvar, align alleles, compute z-scores

# LD vars order
ldvars <- fread(paste0(OUT_PREFIX_LD, ".unphased.vcor1.vars"), header = FALSE)
colnames(ldvars) <- "rsid"

# Read .pvar (skip header lines that start with '#')
pvar <- fread(cmd = paste0("grep -v '^#' ", OUT_PREFIX_PGEN, ".pvar"), header = FALSE, sep = "\t")
pvar <- pvar[, c(1:5)]
colnames(pvar) <- c("CHR", "POS", "ID", "REF", "ALT")

# Merge pvar with GWAS locus stats by ID
pvar2 <- merge(pvar, locus, by = "ID", all.x = TRUE)

# Keep SNPs with available BETA + SE
ix.keep <- which(!is.na(pvar2$BETA) & !is.na(pvar2$SE))
pvar2 <- pvar2[ix.keep, ]

# Align effects: if GWAS effect allele equals REF keep sign, if equals ALT flip sign
pvar2[EFFECT_ALLELE == REF, BETA_aligned := BETA]
pvar2[EFFECT_ALLELE == ALT,  BETA_aligned := -BETA]

# Compute z-score
pvar2[, z := BETA_aligned / SE]
pvar2$z <- as.numeric(pvar2$z)

# Merge into LD order (do not sort; keep LD order)
pvar3 <- merge(ldvars, pvar2, by.x = "rsid", by.y = "ID", all.x = TRUE, sort = FALSE)

# Remove rows with missing z (if any)
pvar3 <- pvar3[!is.na(pvar3$z), ]

# Remove duplicated rsid just in case
pvar3 <- pvar3[!duplicated(pvar3$rsid), ]


# 8) Load LD matrix and subset to match pvar3 ordering

ldmat <- fread(cmd = paste0("zstdcat ", OUT_PREFIX_LD, ".unphased.vcor1.zst"))
ldmat <- as.matrix(ldmat)

# Subset LD matrix to the SNPs retained in pvar3 (same order as ldvars/pvar3)
# ldvars corresponds to original LD order; pvar3 is merged on ldvars, so row order matches ldvars.
# We now subset ldmat rows/cols to those rsids in pvar3$rsid:
keep_idx <- match(pvar3$rsid, ldvars$rsid)
keep_idx <- keep_idx[!is.na(keep_idx)]
ldmat <- ldmat[keep_idx, keep_idx]

# Checking dimension match
if (nrow(ldmat) != nrow(pvar3)) {
  stop("LD matrix dimension does not match number of z-scores after filtering.")
}

# Sample size estimate from GWAS (median N)
if (!("N" %in% colnames(pvar3))) {
  stop("Column N not found in locus/GWAS; required to estimate sample size for susie_rss.")
}
sample_size <- median(pvar3$N, na.rm = TRUE)


# 9) Run SuSiE fine-mapping

susie_fit <- susie_rss(z = pvar3$z, R = ldmat, n = sample_size)

# Credible set indices
cs_sets <- unlist(susie_fit$sets$cs)

# Save outputs
saveRDS(susie_fit, "./output_folder/susie_fit_rs11645785.rds")

# For a simple table of SNP-level PIP:
pip <- susie_fit$pip
pip_tbl <- data.table(
  rsid = pvar3$rsid,
  POS  = pvar3$POS,
  P    = pvar3$P,
  pip  = pip
)
fwrite(pip_tbl, "./output_folder/rs11645785_pip.tsv", sep = "\t")


# 10) Plot fine-mapping results (highlight credible set SNPs)

fine_map_plot <- ggplot(pip_tbl, aes(POS, -log10(P))) +
  geom_point() +
  geom_point(data = pip_tbl[cs_sets, ], aes(POS, -log10(P)), size = 2) +
  theme_classic() +
  labs(
    x = "Position (hg38)",
    y = bquote(-Log[10]~"p-value"),
    title = paste0("SuSiE fine-mapping around ", top_snp)
  ) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

ggsave(
  filename = "./visualization/fine_mapping_rs11645785.pdf",
  plot = fine_map_plot,
  width = 10, height = 6
)