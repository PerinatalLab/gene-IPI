# GxE analysis: SNP × log(IPI), second pregnancy model, (after regenie)

library(data.table)
library(dplyr)
library(ggplot2)

setwd("/mnt/scratch/agnes/gene-IPI/")

df <- fread("results/gwas/regenie/gxe_ipi_2ndpreg_SVLEN_UL_DG.regenie")

# SNP × IPI, + MAF + BETA 
gxe <- df %>%
  filter(TEST == "ADD-INT_SNPxIPI") %>%
  mutate(
    MAF = ifelse(A1FREQ > 0.5, 1 - A1FREQ, A1FREQ)
  ) %>%
  filter(
    MAF >= 0.05,
    abs(BETA) < 10
  ) %>%
  mutate(P = 10^(-LOG10P))

# excluding duplicates
gxe_unique <- gxe %>%
  group_by(ID) %>%
  slice_min(order_by = LOG10P, n = 1) %>%
  ungroup()

# top 20 SNPs 
top_20 <- gxe_unique %>%
  arrange(desc(LOG10P)) %>%
  select(CHROM, GENPOS, ID, A1FREQ, BETA, SE, CHISQ, LOG10P) %>%
  slice(1:20)

write.csv(top_20, "top_20_gxe_ipi2nd_maf05.csv", row.names = FALSE)
print(top_20)

# regional plot at the 1st SNP
top_hit <- top_20[1, ]
region_range <- 250000
region <- gxe %>%
  filter(
    CHROM == top_hit$CHROM,
    GENPOS >= top_hit$GENPOS - region_range,
    GENPOS <= top_hit$GENPOS + region_range
  )

plot(region$GENPOS, region$LOG10P,
     main = paste("Region around", top_hit$ID),
     xlab = "Position", ylab = "-log10(P)",
     pch = 20, col = "blue")
abline(v = top_hit$GENPOS, col = "red", lty = 2)

# QQ plot + lambda GC
gxe_qq <- gxe_unique %>%
  arrange(P) %>%
  mutate(
    observed = LOG10P,
    expected = -log10(ppoints(n())),
    chisq = qchisq(pmax(1 - P, 1e-300), df = 1)
  )

lambda_gc <- median(gxe_qq$chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
cat("Genomic control lambda:", round(lambda_gc, 3), "\n")

ggplot(gxe_qq, aes(x = expected, y = observed)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = paste0("QQ plot - GxE log(IPI), λGC = ", round(lambda_gc, 3)),
    x = "Expected -log10(P)", y = "Observed -log10(P)"
  ) +
  theme_bw()

#  significant SNPs
sig_5e6 <- gxe_unique %>% filter(LOG10P > -log10(5e-6))
sig_1e6 <- gxe_unique %>% filter(LOG10P > -log10(1e-6))

cat("SNPs with P < 5e-6:", nrow(sig_5e6), "\n")
cat("SNPs with P < 1e-6:", nrow(sig_1e6), "\n")

# lambda GC by MAF bins
lambda_by_maf <- gxe_unique %>%
  mutate(MAF_bin = cut(MAF, breaks = seq(0.05, 0.5, by = 0.05))) %>%
  group_by(MAF_bin) %>%
  summarise(lambda = median(qchisq(pmax(1 - 10^(-LOG10P), 1e-300), df = 1)) / qchisq(0.5, df = 1))

ggplot(lambda_by_maf, aes(x = MAF_bin, y = lambda)) +
  geom_bar(stat = "identity") +
  labs(title = "Lambda GC by MAF bin", y = expression(lambda), x = "MAF bin") +
  theme_bw()


# effect sizes
summary(gxe_unique$BETA)
sd(gxe_unique$BETA)

# statistics for P-values
summary(10^(-gxe_unique$LOG10P))
quantile(10^(-gxe_unique$LOG10P), probs = c(0.01, 0.05, 0.10))

# extremely high log10(P) values
table(gxe_unique$LOG10P > 30)
