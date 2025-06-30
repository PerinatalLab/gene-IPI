# multiparous, IPI
library(data.table)
library(dplyr)
library(broom)
library(ggplot2)

df <- fread("results/singel-variants/final-SNP-pheno-multiparous.txt")
df$IPI_days <- df$IPI

snp_cols <- grep("^chr[0-9XYM]+_[0-9]+_[ACGT]+_[ACGT]+$", colnames(df), value = TRUE)

interaction_results <- data.frame()

for (snp in snp_cols) {
  formula_str <- paste0("SVLEN_UL_DG ~ IPI_days * as.factor(`", snp, "`)")
  model <- try(lm(as.formula(formula_str), data = df), silent = TRUE)
  
  if (!inherits(model, "try-error")) {
    res <- tidy(model) %>%
      filter(grepl("IPI_days:as.factor", term)) %>%
      mutate(SNP = snp) %>%
      select(SNP, term, estimate, std.error, p.value)
    
    interaction_results <- bind_rows(interaction_results, res)
  }
}

# p-values
interaction_results <- interaction_results %>%
  arrange(p.value)

# FDR correction Benjamini-Hochberg methods
interaction_results <- interaction_results %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  arrange(p.adj)

# Top SNPs by FDR < 0.05
top_snps_fdr <- interaction_results %>%
  filter(p.adj < 0.05) %>%
  pull(SNP) %>%
  unique()

print(interaction_results %>% filter(p.adj < 0.1))

# top SNPs (p < 0.1) checking
top_snps <- interaction_results %>%
  filter(p.value < 0.1) %>%
  pull(SNP) %>%
  unique()


for (snp in top_snps) {
  cat("  SNP:", snp, " \n")
  model <- lm(as.formula(paste0("SVLEN_UL_DG ~ IPI_days * as.factor(`", snp, "`)")), data = df)
  print(summary(model))  
  
  # plot
  p <- ggplot(df, aes(x = IPI_days, y = SVLEN_UL_DG, color = as.factor(df[[snp]]))) +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_minimal() +
    labs(
      title = "Gestational duration vs IPI_days by SNP dosage",
      subtitle = snp,
      x = "IPI (days)",
      y = "Gestational duration (days)",
      color = "Dosage"
    )
  print(p)
}

# ipi length
print(summary(df$IPI_days))
print(sd(df$IPI_days, na.rm = TRUE))
print(quantile(df$IPI_days, probs = c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), na.rm = TRUE))

# histogram IPI distribution
p_ipi <- ggplot(df, aes(x = IPI_days)) +
  geom_histogram(binwidth = 100, fill = "steelblue", color = "white") +
  labs(
    title = "IPI distribution by day)",
    x = "Interpregnancy interval (day)",
    y = "Frequency"
  ) +
  theme_minimal()
print(p_ipi)


# check..
snp_of_interest <- "chr3_82459425_C_T"

if (!(snp_of_interest %in% colnames(df))) stop("nn")

# IPI statistics by  SNP dosage
ipi_snp_summary <- df %>%
  group_by(Dosage = as.factor(!!sym(snp_of_interest))) %>%
  summarise(
    N = n(),
    Mean_IPI = mean(IPI_days, na.rm = TRUE),
    Median_IPI = median(IPI_days, na.rm = TRUE),
    SD_IPI = sd(IPI_days, na.rm = TRUE),
    Min_IPI = min(IPI_days, na.rm = TRUE),
    Max_IPI = max(IPI_days, na.rm = TRUE)
  )

print(ipi_snp_summary)

ggplot(df, aes(x = as.factor(df[[snp_of_interest]]), y = IPI_days)) +
  geom_boxplot(fill = "skyblue") +
  labs(
    title = paste("IPI lenhlt by dosage:", snp_of_interest),
    x = "SNP genotype (dosage)",
    y = "IPI (days)"
  ) +
  theme_minimal()


for (snp in top_snps) {
  if (!(snp %in% colnames(df))) next
  
  cat("\n==== SNP:", snp, "====\n")
  
  print(df %>%
          group_by(Dosage = as.factor(df[[snp]])) %>%
          summarise(
            N = n(),
            Mean_IPI = mean(IPI_days, na.rm = TRUE),
            Median_IPI = median(IPI_days, na.rm = TRUE),
            SD_IPI = sd(IPI_days, na.rm = TRUE)
          )
  )
  
  print(
    ggplot(df, aes(x = as.factor(df[[snp]]), y = IPI_days)) +
      geom_boxplot(fill = "lightgreen") +
      labs(
        title = paste("IPI distribution by SNP dosage:", snp),
        x = "Genotype",
        y = "IPI (days)"
      ) +
      theme_minimal()
  )
}



#### interaction effect by genotype
summary_table <- data.frame()

for (snp in top_snps) {
  model <- lm(as.formula(paste0("SVLEN_UL_DG ~ IPI_days * as.factor(`", snp, "`)")), data = df)
  coefs <- coef(summary(model))
  
  # IPI  (genotype 0)
  beta_base <- coefs["IPI_days", "Estimate"]
  
  # interaction by genotype
  for (geno_code in c(1, 2)) {
    interaction_rowname <- grep(paste0("IPI_days:as.factor\\(", snp, "\\)", geno_code), rownames(coefs), value = TRUE)
    
    if (length(interaction_rowname) == 1) {
      beta_int <- coefs[interaction_rowname, "Estimate"]
      p_val <- coefs[interaction_rowname, "Pr(>|t|)"]
      total_effect <- beta_base + beta_int
      effect_dir <- ifelse(total_effect > 0, "positive", ifelse(total_effect < 0, "negative", "zero"))
      
      summary_table <- rbind(summary_table, data.frame(
        SNP = snp,
        Genotype = geno_code,
        Base_IPI_effect = beta_base,
        Interaction = beta_int,
        Total_IPI_Effect = total_effect,
        Direction = effect_dir,
        p_value = signif(p_val, 3)
      ))
    }
  }
}


print(summary_table)


summary_table$Genotype <- as.factor(summary_table$Genotype)

ggplot(summary_table, aes(x = Genotype, y = Total_IPI_Effect, fill = Genotype)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  geom_text(aes(label = round(Total_IPI_Effect, 4)), vjust = -0.4, size = 3) +
  facet_wrap(~ SNP, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Genotype-specific IPI effect on gestational duration",
    subtitle = "Slope: change in gestational length per 1-day increase in IPI",
    x = "Genotype (dosage)",
    y = "Total IPI effect (β)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )




# IPI categories (ipi_cat)
ipi_cat <- df %>%
  mutate(
    IPI_Category = cut(
      IPI_days,
      breaks = c(180, 365, 730, 1460, Inf),
      labels = c("6-12 months", "1-2 years", "2-4 years", ">4 years"),
      include.lowest = TRUE
    )
  )

# IPI cat distribution
print(table(ipi_cat$IPI_Category, useNA = "ifany"))

# SNP × IPI cat, by top SNPs
for (snp in top_snps) {
  if (!(snp %in% colnames(ipi_cat))) next
  
  cat("\n==== SNP:", snp, "====\n")
  
  snp_tab <- ipi_cat %>%
    filter(!is.na(IPI_Category)) %>%
    group_by(Dosage = as.factor(ipi_cat[[snp]]), IPI_Category) %>%
    summarise(N = n(), .groups = "drop") %>%
    tidyr::complete(Dosage, IPI_Category, fill = list(N = 0)) %>%
    arrange(Dosage, IPI_Category)
  
  print(snp_tab)
  
  p_snp <- ggplot(snp_tab, aes(x = IPI_Category, y = N, fill = Dosage)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      title = paste("IPI categories by SNP genotype:", snp),
      x = "IPI category",
      y = "Count"
    ) +
    theme_minimal()
  
  print(p_snp)
}




