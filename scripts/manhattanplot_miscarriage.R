library(data.table)
library(dplyr)
library(ggplot2)
library(grid)

pdf("results/summary/miscarriage_manhattan_plot.pdf", width = 8, height = 5)

setwd("~/scratch/agnes/gene-IPI/")
df <- fread("results/gwas/regenie/gxe_miscarriage_gd_gw_SVLEN_UL_DG.regenie")

df <- df %>%
  filter(TEST == "ADD-INT_SNPxmiscarriage") %>%
  mutate(
    MAF = ifelse(A1FREQ > 0.5, 1 - A1FREQ, A1FREQ)
  ) %>%
  filter(MAF >= 0.05) %>%
  rename(CHR = CHROM, POS = GENPOS) %>%
  filter(!is.na(LOG10P), !is.na(CHR), !is.na(POS))

df$CHR2 <- ifelse(as.numeric(as.character(df$CHR)) %% 2 == 0, "even", "odd")

df <- df %>%
  group_by(CHR) %>%
  summarise(chr_len = max(POS)) %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(CHR, tot) %>%
  left_join(df, by = "CHR") %>%
  arrange(CHR, POS) %>%
  mutate(BPcum = POS + tot)

axis_df <- df %>%
  group_by(CHR) %>%
  summarize(center = mean(BPcum))

gw_line <- -log10(5e-8)

label_df <- df %>%
  arrange(desc(LOG10P)) %>%
  slice_head(n = 3) %>%
  mutate(
    x_label = BPcum + c(-4e6, 0, 4e6),
    y_label = LOG10P + c(1.0, 1.4, 1.0)
  )

p <- ggplot(df, aes(x = BPcum, y = LOG10P, color = CHR2)) +
  geom_point(size = 0.3) +
  

  geom_text(
    data = label_df,
    aes(
      x = x_label,
      y = y_label,
      label = ID
    ),
    color = "black", size = 6
  ) +
  
  scale_color_manual(values = c("odd" = "#1f78b4", "even" = "#a6cee3")) +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR) +
  geom_hline(yintercept = gw_line, linetype = "dashed", color = "red") +
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(p))),
    title = "Manhattan plot – miscarriage × SNP G×E (MAF ≥ 0.05)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 12)
  )

print(p)

dev.off()
