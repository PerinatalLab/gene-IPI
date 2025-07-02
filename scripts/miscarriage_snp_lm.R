# miscarriage × SNP GxE interaction with plots
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

# top SNPs 
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

# summary barplot
plot_data <- interaction_results %>%
  filter(SNP %in% top_snps) %>%
  mutate(Dosage = case_when(
    grepl("as.factor\\(.*\\)1$", term) ~ "Het",
    grepl("as.factor\\(.*\\)2$", term) ~ "HomAlt",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(Dosage))

#
if (nrow(plot_data) > 0) {
  ggplot(plot_data, aes(x = SNP, y = estimate, fill = Dosage)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.9) +
    geom_errorbar(aes(
      ymin = estimate - std.error * 1.96,
      ymax = estimate + std.error * 1.96
    ), position = position_dodge(width = 0.7), width = 0.2) +
    labs(
      title = "Interaction effect: miscarriage × SNP genotype",
      x = NULL,
      y = "Effect estimate (days)",
      fill = "Genotype"
    ) +
    scale_fill_manual(
      values = c("Het" = "orange", "HomAlt" = "tomato"),
      labels = c("Het" = "Heterozygous (1)", "HomAlt" = "Homozygous Alt (2)")
    ) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
} else {
  message("none")
}





