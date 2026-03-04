library(dplyr)
library(data.table)
library(ggplot2)
library(grid)

dir.create("visualization", showWarnings = FALSE, recursive = TRUE)

#Reading the summary statistics from all file
gwas_all <- fread("./GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL")
head(gwas_all)
colnames(gwas_all)

#1. we will first visualize the overall GWAS result for all ancestry using the manhatton plot
#converting the p value into numeric incase it is not
gwas_all$P<- as.numeric(gwas_all$P)

#Manhattan plot
manhat<-ggplot(gwas_all, aes(x = POS, y = -log10(P))) +
  geom_point(alpha = 0.5, size = 0.6) +
  facet_wrap(~CHR, scales = "free_x", nrow = 2) +
  geom_hline(yintercept = -log10(5e-8), color = "red") +
  theme_classic() +
  labs(title = "Manhattan Plot for Height GWAS",
       x = "Genomic Position",
       y = expression(-log[10](p-value))) +
  theme(axis.text.x = element_blank(),
        panel.spacing = unit(0.3, "lines")) + 
   theme(
  axis.text.x = element_text(size = 12, face = "bold"),
  axis.text.y = element_text(size = 12, face = "bold"),     
  axis.title.x = element_text(size = 12, face = "bold"),    
  axis.title.y = element_text(size = 12, face = "bold"),    
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
ggsave(filename = "./visualization/Manhattan_plot_for_all_chromosomes_GWAS.pdf",
       plot = manhat,
       width = 12, height = 10)
manhat

#Lets make QQ plot as well
# Sort by p-value
gwas_all <- gwas_all[order(P)]

# Expected p-values under null
gwas_all$expected <- -log10(ppoints(nrow(gwas_all)))

# Observed p-values
gwas_all$observed <- -log10(gwas_all$P)

# QQ plot
qq_plot<-ggplot(gwas_all, aes(x = expected, y = observed)) +
  geom_point(alpha = 0.5, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  theme_classic() +
  labs(title = "QQ Plot for Height GWAS",
       x = "Expected -log10(p)",
       y = "Observed -log10(p)") + 
   theme(
  axis.text.x = element_text(size = 12, face = "bold"),
  axis.text.y = element_text(size = 12, face = "bold"),     
  axis.title.x = element_text(size = 12, face = "bold"),    
  axis.title.y = element_text(size = 12, face = "bold"),    
  plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

ggsave(filename = "./visualization/QQplot_for_all_chromosomes_GWAS.pdf",
       plot = qq_plot,
       width = 10, height = 8)
qq_plot


