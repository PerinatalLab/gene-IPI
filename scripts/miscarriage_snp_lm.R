# miscarriage × SNP GxE interaction with plots -- by additive model
library(data.table)
library(dplyr)
library(ggplot2)
library(broom)

setwd("~/scratch/agnes/gene-IPI")

df <- fread("results/singel-variants/final-SNP-pheno-nulliparous.txt")
df$miscarriage_bin <- ifelse(df$miscarriage > 0, 1, 0)

snp_cols <- grep("^chr[0-9XYM]+_[0-9]+_[ACGT]+_[ACGT]+$", names(df), value = TRUE)

interaction_results <- data.frame()

# interactions (additive model)
for (snp in snp_cols) {
  df$SNP_dosage <- df[[snp]]  # additive coding: 0, 1, 2
  
  model <- try(lm(SVLEN_UL_DG ~ miscarriage_bin * SNP_dosage, data = df), silent = TRUE)
  
  if (!inherits(model, "try-error")) {
    res <- tidy(model) %>%
      filter(term == "miscarriage_bin:SNP_dosage") %>%
      mutate(SNP = snp) %>%
      select(SNP, term, estimate, std.error, p.value)
    
    interaction_results <- bind_rows(interaction_results, res)
  }
}

# FDR correction
interaction_results <- interaction_results %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  arrange(p.adj)

# top SNPs (FDR < 0.1)
top_snps <- interaction_results %>%
  filter(p.adj < 0.1) %>%
  pull(SNP) %>%
  unique()

cat("Top SNPs (FDR < 0.1):", length(top_snps), "\n")

# top SNPs – descriptive plots
for (snp in top_snps) {
  cat("\n\n=== SNP:", snp, "===\n")
  
  df$genotype_numeric <- df[[snp]]
  df$genotype_label <- factor(df$genotype_numeric, levels = c(0, 1, 2),
                              labels = c("Hom. Ref", "Het", "Hom. Alt"))
  
  df_plot <- df %>%
    filter(!is.na(SVLEN_UL_DG), !is.na(miscarriage_bin), !is.na(genotype_label))
  
  model <- lm(SVLEN_UL_DG ~ miscarriage_bin * genotype_label, data = df_plot)
  print(summary(model))
  
  # stat table
  stat_table <- df_plot %>%
    group_by(miscarriage_bin, genotype_label) %>%
    summarise(
      mean = mean(SVLEN_UL_DG),
      sd = sd(SVLEN_UL_DG),
      n = n(),
      .groups = "drop"
    )
  print(stat_table)
  
  # barplot (mean ± SD)
  p_bar <- ggplot(stat_table, aes(x = genotype_label, y = mean, fill = factor(miscarriage_bin))) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                  position = position_dodge(width = 0.8), width = 0.2) +
    labs(
      title = paste("Gestational duration by genotype and miscarriage:", snp),
      x = "Genotype",
      y = "Mean gestational duration (days)",
      fill = "Miscarriage"
    ) +
    scale_fill_manual(values = c("0" = "#a6cee3", "1" = "#fb9a99"),
                      labels = c("0" = "No", "1" = "Yes")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top")
  
  print(p_bar)
  
  # LM plot: regressions
  p_lm <- ggplot(df_plot, aes(x = miscarriage_bin, y = SVLEN_UL_DG, color = genotype_label)) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    geom_smooth(method = "lm", se = TRUE, aes(group = genotype_label), formula = y ~ x) +
    scale_x_continuous(breaks = c(0, 1), labels = c("No miscarriage", "Miscarriage")) +
    labs(
      title = paste("Interaction plot (lm):", snp),
      x = "Miscarriage",
      y = "Gestational duration (days)",
      color = "Genotype"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top")
  
  print(p_lm)
}

# summary table + barplot of dosage effects by miscarriage group
summary_table <- data.frame()

for (snp in top_snps) {
  df$SNP_dosage <- df[[snp]]
  model <- lm(SVLEN_UL_DG ~ miscarriage_bin * SNP_dosage, data = df)
  coefs <- coef(summary(model))
  
  # baseline: miscarriage_bin == 0
  beta_base <- coefs["SNP_dosage", "Estimate"]
  
  # interaction
  interaction_row <- "miscarriage_bin:SNP_dosage"
  if (interaction_row %in% rownames(coefs)) {
    beta_int <- coefs[interaction_row, "Estimate"]
    p_val <- coefs[interaction_row, "Pr(>|t|)"]
    total_effect <- beta_base + beta_int
    
    summary_table <- rbind(summary_table, data.frame(
      SNP = snp,
      Group = c("No miscarriage", "Miscarriage"),
      Effect = c(beta_base, total_effect),
      p_value = c(NA, signif(p_val, 3))
    ))
  }
}

summary_table$Group <- factor(summary_table$Group, levels = c("No miscarriage", "Miscarriage"))

# summary barplot (facet per SNP)
p_summary <- ggplot(summary_table, aes(x = Group, y = Effect, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  geom_text(aes(label = round(Effect, 3)), vjust = -0.5, size = 3) +
  facet_wrap(~ SNP, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Effect of SNP dosage on gestational duration",
    subtitle = "Stratified by miscarriage history (additive G×E model)",
    x = "Miscarriage group",
    y = "Effect estimate (days)"
  ) +
  scale_fill_manual(values = c("No miscarriage" = "#a6cee3", "Miscarriage" = "#fb9a99")) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  )

print(p_summary)
