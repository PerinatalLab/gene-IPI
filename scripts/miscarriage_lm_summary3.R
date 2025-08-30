# nulliparous summary lm, miscarriage plots 
library(data.table)
library(dplyr)
library(ggplot2)
library(broom)

setwd("~/scratch/agnes/gene-IPI")
df <- fread("results/singel-variants/final-SNP-pheno-nulliparous.txt")
df$miscarriage_bin <- ifelse(df$miscarriage > 0, 1, 0)

snp_cols <- grep("^chr[0-9XYM]+_[0-9]+_[ACGT]+_[ACGT]+$", names(df), value = TRUE)
interaction_results <- data.frame()

for (snp in snp_cols) {
  df$SNP_dosage <- df[[snp]]
  model <- try(lm(SVLEN_UL_DG ~ miscarriage_bin * SNP_dosage, data = df), silent = TRUE)
  if (!inherits(model, "try-error")) {
    res <- tidy(model) %>%
      filter(term == "miscarriage_bin:SNP_dosage") %>%
      mutate(SNP = snp) %>%
      select(SNP, term, estimate, std.error, p.value)
    interaction_results <- bind_rows(interaction_results, res)
  }
}

interaction_results <- interaction_results %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  arrange(p.adj)

top_snps <- interaction_results %>%
  filter(p.adj < 0.1) %>%
  pull(SNP) %>%
  unique()

cat("Top SNPs (FDR < 0.1):", length(top_snps), "\n")
if (length(top_snps)) print(interaction_results %>% filter(SNP %in% top_snps))

for (snp in top_snps) {
  cat("\n\n=== SNP:", snp, "===\n")
  df$genotype_numeric <- df[[snp]]
  df$genotype_label <- factor(df$genotype_numeric, levels = c(0, 1, 2),
                              labels = c("Hom. Ref", "Het", "Hom. Alt"))
  df_plot <- df %>%
    filter(!is.na(SVLEN_UL_DG), !is.na(miscarriage_bin), !is.na(genotype_label))
  model <- lm(SVLEN_UL_DG ~ miscarriage_bin * genotype_label, data = df_plot)
  print(summary(model))
  stat_table <- df_plot %>%
    group_by(miscarriage_bin, genotype_label) %>%
    summarise(mean = mean(SVLEN_UL_DG), sd = sd(SVLEN_UL_DG), n = n(), .groups = "drop")
  print(stat_table)
  p_bar <- ggplot(stat_table, aes(x = genotype_label, y = mean, fill = factor(miscarriage_bin))) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                  position = position_dodge(width = 0.8), width = 0.2) +
    labs(title = paste("Gestational duration by genotype and miscarriage:", snp),
         x = "Genotype", y = "Mean gestational duration (days)", fill = "Miscarriage") +
    scale_fill_manual(values = c("0" = "#a6cee3", "1" = "#fb9a99"),
                      labels = c("0" = "No", "1" = "Yes")) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top")
  print(p_bar)
  p_lm <- ggplot(df_plot, aes(x = miscarriage_bin, y = SVLEN_UL_DG, color = genotype_label)) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    geom_smooth(method = "lm", se = TRUE, aes(group = genotype_label), formula = y ~ x) +
    scale_x_continuous(breaks = c(0, 1), labels = c("No miscarriage", "Miscarriage")) +
    labs(title = paste("Interaction plot (lm):", snp),
         x = "Miscarriage", y = "Gestational duration (days)", color = "Genotype") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top")
  print(p_lm)
}

summary_table <- data.frame()
top_snps_chr <- top_snps[grepl("^chr(10|11)_", top_snps)]
cat("Summary plot SNPs (chr10/chr11):", ifelse(length(top_snps_chr)>0, paste(top_snps_chr, collapse=", "), "—"), "\n")

for (snp in top_snps_chr) {
  df$SNP_dosage <- df[[snp]]
  model <- lm(SVLEN_UL_DG ~ miscarriage_bin * SNP_dosage, data = df)
  coefs <- coef(summary(model))
  if (!("SNP_dosage" %in% rownames(coefs))) next
  beta_base <- coefs["SNP_dosage","Estimate"]
  if ("miscarriage_bin:SNP_dosage" %in% rownames(coefs)) {
    beta_int <- coefs["miscarriage_bin:SNP_dosage","Estimate"]
    p_val <- coefs["miscarriage_bin:SNP_dosage","Pr(>|t|)"]
    total_effect <- beta_base + beta_int
    summary_table <- rbind(summary_table,
                           data.frame(SNP = snp, Group = "No miscarriage", Effect = beta_base, p_value = NA),
                           data.frame(SNP = snp, Group = "Miscarriage", Effect = total_effect, p_value = signif(p_val, 3)))
  }
}

if (nrow(summary_table) > 0) {
  summary_table$Group <- factor(summary_table$Group, levels = c("No miscarriage","Miscarriage"))
  p_summary <- ggplot(summary_table, aes(x = Group, y = Effect, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.6) +
    geom_text(aes(label = round(Effect, 3)), vjust = -0.5, size = 3) +
    facet_wrap(~ SNP, scales = "free_y") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Effect of SNP dosage on gestational duration",
         subtitle = "Stratified by miscarriage history (chr10 & chr11 only)",
         x = "Miscarriage group", y = "Effect estimate (days)") +
    scale_fill_manual(values = c("No miscarriage" = "#a6cee3", "Miscarriage" = "#fb9a99")) +
    theme_minimal() +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))
  print(p_summary)


