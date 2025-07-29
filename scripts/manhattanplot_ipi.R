library(data.table)
library(dplyr)
library(ggplot2)

pdf("results/summary/ipi_manhattan_plot.pdf", width = 8, height = 5)

#setwd('/mnt/scratch/agnes/gene-IPI/')
df <- fread("results/gwas/regenie/gxe_ipi_gd_gw_SVLEN_UL_DG.regenie")

df <- df %>%
  filter(TEST == "ADD-INT_SNPxIPI") %>%
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

top_snps <- df %>% arrange(desc(LOG10P)) %>% slice_head(n = 10)
df$label <- ifelse(df$ID %in% top_snps$ID, df$ID, NA)

gw_line <- -log10(5e-8)



# top 5 SNP
top_snps <- df %>% arrange(desc(LOG10P)) %>% slice_head(n = 3)

label_df <- top_snps %>%
  mutate(nudge = c(0.3, 0.6, 0.9)) 

# plot
p <- ggplot(df, aes(x = BPcum, y = LOG10P, color = CHR2)) +
  geom_point(size = 0.3) +
  geom_text(data = label_df,
            aes(label = ID, y = LOG10P + nudge),
            color = "black", size = 2.5) +
  scale_color_manual(values = c("odd" = "gray40", "even" = "gray70")) +
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR) +
  geom_hline(yintercept = gw_line, linetype = "dashed", color = "red") +
  labs(x = "Chromosome", y = expression(-log[10](italic(p))),
       title = "Manhattan plot – IPI × SNP G×E (MAF ≥ 0.05)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")
print(p)

 dev.off()