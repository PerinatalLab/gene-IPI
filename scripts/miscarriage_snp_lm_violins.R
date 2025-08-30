# miscarriage × SNP GxE interaction with plots -- by additive model === plots --> 1 PDF file
library(data.table)
library(dplyr)
library(ggplot2)
library(broom)


collapse_mode <- "always" 
n_thresh <- 30             

pdf("results/summary/miscarriage_interaction_plots.pdf", width = 10, height = 6)
setwd("~/scratch/agnes/gene-IPI")

df <- fread("results/singel-variants/final-SNP-pheno-nulliparous.txt")
df$miscarriage_bin <- ifelse(df$miscarriage > 0, 1, 0)

snp_cols <- grep("^chr[0-9XYM]+_[0-9]+_[ACGT]+_[ACGT]+$", names(df), value = TRUE)
interaction_results <- data.frame()

# interactions (adittive, SNP-dosage 0/1/2)

for (snp in snp_cols) {
  df$SNP_dosage <- df[[snp]]  # 0/1/2
  model <- try(lm(SVLEN_UL_DG ~ miscarriage_bin * SNP_dosage, data = df), silent = TRUE)
  if (!inherits(model, "try-error")) {
    res <- tidy(model) %>%
      filter(term == "miscarriage_bin:SNP_dosage") %>%
      mutate(SNP = snp) %>%
      select(SNP, term, estimate, std.error, p.value)
    interaction_results <- bind_rows(interaction_results, res)
  }
}

# FDR
interaction_results <- interaction_results %>%
  mutate(p.adj = p.adjust(p.value, method = "fdr")) %>%
  arrange(p.adj)

# top SNPs (FDR < 0.1)
top_snps <- interaction_results %>%
  filter(p.adj < 0.1) %>% pull(SNP) %>% unique()
cat("Top SNPs (FDR < 0.1):", length(top_snps), "\n")


# descriptive and plots

if (length(top_snps) > 0) {
  for (snp in top_snps) {
    cat("\n\n=== SNP:", snp, "===\n")
    
    df$genotype_numeric <- df[[snp]]
    df$genotype_label   <- factor(df$genotype_numeric, levels = c(0,1,2),
                                  labels = c("Hom. Ref","Het","Hom. Alt"))
    df_plot <- df %>% filter(!is.na(SVLEN_UL_DG), !is.na(miscarriage_bin), !is.na(genotype_label))
    
    # categorical descriptive model
    model_cat <- lm(SVLEN_UL_DG ~ miscarriage_bin * genotype_label, data = df_plot)
    print(summary(model_cat))
    
    # summary stat
    stat_table <- df_plot %>%
      group_by(miscarriage_bin, genotype_label) %>%
      summarise(mean = mean(SVLEN_UL_DG), sd = sd(SVLEN_UL_DG), n = n(), .groups = "drop")
    print(stat_table)
    
  

    df_plot2 <- df_plot %>%
      mutate(
        miscarriage_label = factor(miscarriage_bin, levels = c(0,1),
                                   labels = c("No miscarriage","Miscarriage"))
      )
    
    if (collapse_mode == "auto") {
      cell_counts <- df_plot2 %>%
        group_by(genotype_label, miscarriage_label) %>%
        summarise(n = n(), .groups = "drop")
      need_collapse <- any(cell_counts$n < n_thresh)
    } else {
      need_collapse <- TRUE
    }
    
    if (need_collapse) {
      df_plot2 <- df_plot2 %>%
        mutate(geno_show = ifelse(genotype_numeric == 0, "Non-carrier", "Carrier (>=1 alt)"))
      df_plot2$geno_show <- factor(df_plot2$geno_show,
                                   levels = c("Non-carrier","Carrier (>=1 alt)"))
      x_lab <- "Genotype (collapsed)"
    } else {
      df_plot2 <- df_plot2 %>%
        mutate(geno_show = as.character(genotype_label))
      df_plot2$geno_show <- factor(df_plot2$geno_show,
                                   levels = c("Hom. Ref","Het","Hom. Alt"))
      x_lab <- "Genotype"
    }

    n_lab <- df_plot2 %>%
      group_by(miscarriage_label, geno_show) %>%
      summarise(n = n(), .groups = "drop")
    
    y_max <- max(df_plot2$SVLEN_UL_DG, na.rm = TRUE)
    y_min <- min(df_plot2$SVLEN_UL_DG, na.rm = TRUE)
    
 
    p_int <- tryCatch({
      anova(lm(SVLEN_UL_DG ~ miscarriage_label * geno_show, data = df_plot2))[
        "miscarriage_label:geno_show","Pr(>F)"]
    }, error = function(e) NA_real_)
    subtxt <- if (need_collapse)
      paste0("Collapsed: non-carrier (0) vs carrier (1/2)",
             if (!is.na(p_int)) paste0("   |   global interaction p=", signif(p_int, 3)) else "")
    else if (!is.na(p_int)) paste0("Global interaction p=", signif(p_int, 3)) else NULL
    
    p_violin <- ggplot(df_plot2, aes(x = geno_show, y = SVLEN_UL_DG, fill = geno_show)) +
      geom_violin(trim = FALSE, alpha = 0.7, width = 0.7) + 
      geom_boxplot(width = 0.12, outlier.shape = NA) +
      stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "white") +
      geom_text(data = n_lab,
                aes(x = geno_show, y = y_max * 1.008, label = paste0("n = ", n)),
                inherit.aes = FALSE, size = 3) +
      coord_cartesian(ylim = c(y_min, y_max * 1.04)) +
      scale_x_discrete(drop = FALSE) +
      facet_grid(. ~ miscarriage_label) +
      labs(
        title = paste("Gestational duration by miscarriage status and genotype:", snp),
        subtitle = subtxt,
        x = x_lab, y = "Gestational duration (days)",
        fill = "Genotype"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "top", strip.text = element_text(face = "bold"))
    
    print(p_violin)
    

  }
} else {
  message("nnincs top SNP (FDR < 0.1).")
}


# summary (dosage-effect by miscarriage-groups)

summary_table <- data.frame()
for (snp in top_snps) {
  df$SNP_dosage <- df[[snp]]
  model <- lm(SVLEN_UL_DG ~ miscarriage_bin * SNP_dosage, data = df)
  coefs <- coef(summary(model))
  beta_base <- coefs["SNP_dosage","Estimate"]
  interaction_row <- "miscarriage_bin:SNP_dosage"
  if (interaction_row %in% rownames(coefs)) {
    beta_int <- coefs[interaction_row,"Estimate"]
    p_val <- coefs[interaction_row,"Pr(>|t|)"]
    total_effect <- beta_base + beta_int
    summary_table <- rbind(summary_table, data.frame(
      SNP = snp,
      Group = c("No miscarriage","Miscarriage"),
      Effect = c(beta_base, total_effect),
      p_value = c(NA, signif(p_val, 3))
    ))
  }
}
summary_table$Group <- factor(summary_table$Group, levels = c("No miscarriage","Miscarriage"))

p_summary <- ggplot(summary_table, aes(x = Group, y = Effect, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  geom_text(aes(label = round(Effect, 3)), vjust = -0.5, size = 3) +
  facet_wrap(~ SNP, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Effect of SNP dosage on gestational duration",
    subtitle = "Stratified by miscarriage history",
    x = "Miscarriage group", y = "Effect estimate (days)"
  ) +
  theme_minimal() +
  theme(legend.position = "none", strip.text = element_text(face = "bold"))

print(p_summary)


fwrite(interaction_results, file = "results/summary/miscarriage_interaction_results.csv")
dev.off()
