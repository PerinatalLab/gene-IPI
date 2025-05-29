# QC ?? wide GWAS 
library(ggplot2)
library(dplyr)
library(data.table)


gxe_data <- fread("results/gwas/gxe_ipi_gd_gw_screen.SVLEN_UL_DG.glm.linear")

# IPI × SNP test
gxe_addxipi <- gxe_data %>%
  filter(TEST == "ADDxIPI") %>%
  mutate(
    MAF = ifelse(A1_FREQ > 0.5, 1 - A1_FREQ, A1_FREQ)
  ) %>%
  filter(
    MAF >= 0.01,
    OBS_CT >= 1000,
    abs(BETA) < 10
  )

# top 10 SNPs -- raw, for checking
top_hits <- gxe_addxipi %>%
  arrange(P) %>%
  select(`#CHROM`, POS, ID, BETA, SE, T_STAT, P) %>%
  head(10)
print(top_hits)

# to remove doubled SNP IDs
duplicated_snps <- gxe_addxipi %>% filter(duplicated(ID))
print(duplicated_snps)

# one of the most significant doubled SNP IDs 
gxe_addxipi_unique <- gxe_addxipi %>%
  group_by(ID) %>%
  slice_min(order_by = P, n = 1) %>%
  ungroup()

# top 10 SNP -- table
top_hits_clean <- gxe_addxipi_unique %>%
  arrange(P) %>%
  select(`#CHROM`, POS, ID, MAF, OBS_CT, BETA, SE, T_STAT, P) %>%
  head(10)
print(top_hits_clean)

# regional plot
top_snps <- top_hits_clean %>% slice(1)
region_range <- 250000

region1 <- gxe_addxipi %>%
  filter(
    `#CHROM` == top_snps$`#CHROM`[1],
    POS >= top_snps$POS[1] - region_range,
    POS <= top_snps$POS[1] + region_range
  )

plot(region1$POS, -log10(region1$P),
     main = paste("Regional plot around", top_snps$ID[1]),
     xlab = "Position", ylab = "-log10(P)",
     pch = 20, col = "blue")
abline(v = top_snps$POS[1], col = "red", lty = 2)

# qq plot // components
gxe_qq <- gxe_data %>%
  filter(TEST == "ADDxIPI", !is.na(P), P > 0) %>%
  arrange(P) %>%
  mutate(
    observed = -log10(P),
    expected = -log10(ppoints(n())),
    chisq = qchisq(1 - P, df = 1)
  )

# genomic control lambda estimation
lambda_gc <- median(gxe_qq$chisq, na.rm = TRUE) / qchisq(0.5, df = 1)

# qq plot 
ggplot(gxe_qq, aes(x = expected, y = observed)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = paste0("QQ Plot - GxE (ADDxIPI), λGC = ", round(lambda_gc, 3)),
    x = "Expected -log10(P)", y = "Observed -log10(P)"
  ) +
  theme_bw()

# genome-wide significant SNPs
significant_hits <- gxe_addxipi_unique %>%
  filter(P < 5e-8)
cat("Number of significant SNPs (P < 5e-8):", nrow(significant_hits), "\n")

# top 20 SNP
top_20 <- gxe_addxipi_unique %>%
  arrange(P) %>%
  slice(1:20)
print(top_20)


write.table(top_hits_clean, "top_snps_clean.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top_hits_clean$ID, "top_rsids.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(top_20, "top_20_snps_clean.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(top_20$ID, "top20_rsids.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
