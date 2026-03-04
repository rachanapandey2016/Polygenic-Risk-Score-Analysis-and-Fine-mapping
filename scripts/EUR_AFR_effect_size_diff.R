library(dplyr)
library(data.table)
library(ggplot2)

# Read in EUR and AFR GWAS files
eur <- fread("./GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR")
afr <- fread("./GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_AA")

#It is better to check how many rows have NA and then remove such rows
eur <- eur %>% filter(!is.na(P), !is.na(BETA))
afr <- afr %>% filter(!is.na(P), !is.na(BETA))

head(eur)
head(afr)

#Removing duplicate just to ensure
# For EUR GWAS
eur<- eur %>%
  filter(!duplicated(RSID))

# For AFR GWAS
afr<- afr %>%
  filter(!duplicated(RSID))

fwrite(eur, "./eur_cleaned.tsv", sep = "\t")

#Merge european and african dataset
# Merge based on SNP ID
merged <- inner_join(eur, afr, by = "RSID", suffix = c("_EUR", "_AFR"))

# Calculate effect size difference
merged <- merged %>%
  mutate(beta_diff = BETA_EUR - BETA_AFR)

#Visualize differences- how much individual SNP associations differ between EUR and AFR
height_effect_size<- ggplot(merged, aes(x = beta_diff)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black") +
  labs(title = "Distribution of Height Effect Size Differences (EUR vs AFR)",
       x = "Effect Size Difference (Beta_EUR - Beta_AFR)",
       y = "Count") + 
   theme(
  axis.text.x = element_text(size = 12, face = "bold"),
  axis.text.y = element_text(size = 12, face = "bold"),     
  axis.title.x = element_text(size = 12, face = "bold"),    
  axis.title.y = element_text(size = 12, face = "bold"),    
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

ggsave(filename = "./visualization/Distribution of Height Effect Size Differences.pdf",
       plot = height_effect_size,
       width = 10, height = 8)
height_effect_size

#Snps with big differences
# Top SNPs with largest absolute beta differences
top_diff_snps <- merged %>%
  arrange(desc(abs(beta_diff))) %>%
  select(RSID, BETA_EUR, BETA_AFR, beta_diff) %>%
  head(10)

print(top_diff_snps)

# Bar plot of effect size differences
top_SNP<-ggplot(top_diff_snps, aes(x = reorder(RSID, abs(beta_diff)), y = beta_diff, fill = beta_diff > 0)) +
  geom_col() +
  coord_flip() +  # Flip for easier reading
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "tomato"), name = "Direction (EUR - AFR)") +
  labs(title = "Top 10 SNPs with Largest Effect Size Differences (EUR vs AFR)",
       x = "SNP (RSID)",
       y = "Effect Size Difference (Beta_EUR - Beta_AFR)") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(face = "bold")
  )
ggsave(filename = "./visualization/Top 10 SNPs with largest Effect Size Differences.pdf",
       plot = top_SNP,
       width = 10, height = 8)
top_SNP