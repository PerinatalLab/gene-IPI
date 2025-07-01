# GxE IPI, regenie, MAF 0.05, test for robustness 

library(data.table)
library(dplyr)
library(ggplot2)


#df <- fread("~/scratch/agnes/gene-IPI/results/gwas/regenie/gxe_ipi_gd_gw_SVLEN_UL_DG.regenie")
df <- fread(snakemake@input[[1]])
            
# SNP × IPI
gxe <- df %>%
  filter(TEST == "ADD-INT_SNPxIPI") %>%
  mutate(
    MAF = ifelse(A1FREQ > 0.5, 1 - A1FREQ, A1FREQ)
  ) %>%
  filter(
    MAF >= 0.05, 
    abs(BETA) < 10
  ) %>%
  mutate(
    P = 10^(-LOG10P)
  )

# handling duplicant SNPs, keeping the most significant ones
gxe_unique <- gxe %>%
  group_by(ID) %>%
  slice_min(order_by = LOG10P, n = 1) %>%
  ungroup()

# top 20 SNP
top_20 <- gxe_unique %>%  
  arrange(desc(LOG10P)) %>%
  select(CHROM, GENPOS, ID, A1FREQ, BETA, SE, CHISQ, LOG10P) %>%
  slice(1:20)
print(top_20)

# regional plot around the first hit
top_hit <- top_20[1, ]
region_range <- 250000
region <- gxe %>%
  filter(
    CHROM == top_hit$CHROM,
    GENPOS >= top_hit$GENPOS - region_range,
    GENPOS <= top_hit$GENPOS + region_range
  )

png(filename = snakemake@output[[2]], width = 800, height = 600)

plot(region$GENPOS, region$LOG10P,
     main = paste("Region around", top_hit$ID),
     xlab = "Position", ylab = "-log10(P)",
     pch = 20, col = "blue")
abline(v = top_hit$GENPOS, col = "red", lty = 2)



# QQ plot components + lambda GC
gxe_qq <- gxe_unique %>%
  arrange(P) %>%
  mutate(
    observed = LOG10P,
    expected = -log10(ppoints(n())),
    chisq = qchisq(pmax(1 - P, 1e-300), df = 1)
  )

lambda_gc <- median(gxe_qq$chisq, na.rm = TRUE) / qchisq(0.5, df = 1)
cat("Genomic control lambda:", round(lambda_gc, 3), "\n")

# QQ plot
ggplot(gxe_qq, aes(x = expected, y = observed)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = paste0("QQ plot - Regenie GxE (SNP × IPI), λGC = ", round(lambda_gc, 3)),
    x = "Expected -log10(P)", y = "Observed -log10(P)"
  ) +
  theme_bw()

ggsave(filename = snakemake@output[[3]], plot = last_plot(), width = 8, height = 6)

# significant SNPs
sig_5e6 <- gxe_unique %>% filter(LOG10P > -log10(5e-6))
sig_1e6 <- gxe_unique %>% filter(LOG10P > -log10(1e-6))

cat("SNPs with P < 5e-6:", nrow(sig_5e6), "\n")
cat("SNPs with P < 1e-6:", nrow(sig_1e6), "\n")

#getwd()
# write.table(sig_1e6, "top_gxe_ipi_p1e6.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(sig_5e6, "top_gxe_ipi_p5e6.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# write.table(top_20, "top20_gxe_ipi_regenie.txt", sep = "\t", quote = FALSE, row.names = FALSE)
# write.csv(top_20, file = "~/scratch/agnes/gene-IPI/summary/ipi_regenie_top_20_gxe_maf05.csv", row.names = FALSE)
write.csv(top_20, file = snakemake@output[[1]], row.names = FALSE)

summary(gxe_unique$BETA)
sd(gxe_unique$BETA)
summary(10^(-gxe_unique$LOG10P))
quantile(10^(-gxe_unique$LOG10P), probs = c(0.01, 0.05, 0.10))


table(gxe_unique$LOG10P > 30)

# lambda GC by MAF bins
gxe_unique %>%
  mutate(MAF_bin = cut(MAF, breaks = seq(0.05, 0.5, by = 0.05))) %>%
  group_by(MAF_bin) %>%
  summarise(lambda = median(qchisq(pmax(1 - 10^(-LOG10P), 1e-300), df = 1)) / qchisq(0.5, df = 1)) %>%
  ggplot(aes(x = MAF_bin, y = lambda)) +
  geom_bar(stat = "identity") +
  labs(title = "Lambda GC by MAF bin", y = expression(lambda), x = "MAF bin") +
  theme_bw()

#file.show("~/scratch/agnes/gene-IPI/top_20_gxe_maf05.csv")



# investigating the p-palue
ggplot(gxe_unique, aes(x = LOG10P)) +
  geom_histogram(bins = 60, fill = "lightblue", color = "black") +
  labs(title = "-log10(p) distribution (GxE SNP × IPI)",
       x = expression(-log[10](p)), y = "Count") +
  theme_bw()

ggplot(gxe_unique, aes(x = P)) +
  geom_histogram(bins = 60, fill = "lightblue", color = "black") +
  scale_x_log10() +
  labs(title = "Raw p-value distribution (GxE SNP × IPI)",
       x = "P-value (log10 scale)", y = "Count") +
  theme_bw()

ggplot(gxe_unique, aes(x = P)) +
  stat_ecdf(geom = "step", color = "lightblue") +
  scale_x_log10() +
  labs(title = "Cumulative distribution of p-values (GxE SNP × IPI)",
       x = "P-value (log10 scale)", y = "Cumulative fraction") +
  theme_bw()
