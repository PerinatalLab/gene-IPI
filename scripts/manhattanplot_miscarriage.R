library(data.table)
library(dplyr)
library(ggplot2)

pdf("results/summary/miscarriage_manhattan_plot.pdf", width = 10, height = 5)

#setwd("~/scratch/agnes/gene-IPI/")
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

top_snps <- df %>% arrange(desc(LOG10P)) %>% slice_head(n = 3)

label_df <- top_snps %>%
  mutate(nudge = seq(0.3, 1.5, length.out = n()))

gw_line <- -log10(5e-8)

# Manhattan plot
label_df <- df %>%
  arrange(desc(LOG10P)) %>%
  slice_head(n = 5)

label_df$y_offset <- 0.3

duplicated_chrs <- label_df %>%
  count(CHR) %>%
  filter(n > 1) %>%
  pull(CHR)

for (chr in duplicated_chrs) {
  idx <- which(label_df$CHR == chr)
  
  if (length(idx) == 2 && all(idx %in% 1:2)) {
    label_df$y_offset[idx] <- c(0.5, 1.0)
  } else {
    label_df$y_offset[idx] <- seq(0.3, 0.8, length.out = length(idx))
  }
}

p <- ggplot(df, aes(x = BPcum, y = LOG10P, color = CHR2)) +
  geom_point(size = 0.3) +
  geom_text(data = label_df,
            aes(x = BPcum,
                y = LOG10P + y_offset,
                label = ID),
            color = "black", size = 2.5) +
  scale_color_manual(values = c("odd" = "gray40", "even" = "gray70")) +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR) +
  geom_hline(yintercept = gw_line, linetype = "dashed", color = "red") +
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(p))),
    title = "Manhattan plot – miscarriage × SNP G×E (MAF ≥ 0.05)"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

print(p)

# ggsave("results/summary/miscarriage_manhattan_plot_top5_maf05.png", plot = p, width = 10, height = 5, dpi = 300)
dev.off()