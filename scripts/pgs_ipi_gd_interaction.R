library(data.table)
library(ggplot2)
library(dplyr)

data <- data.table::fread("~/scratch/agnes/gene-IPI/results/final_phenotype/IPI_pgs_covariates.txt")
str(data)
summary(data)



# ipi groups
data$IPI_group <- cut(data$IPI, breaks = c(-Inf, 365, 1800, Inf), labels = c("Short", "Medium", "Long"))

data %>%
  group_by(IPI_group) %>%
  summarise(N = n())

summary(data)

ipi_counts <- data[, .N, by = IPI_group]
ipi_counts[, Percentage := round(100 * N / sum(N), 1)]
ipi_counts



# PGS standars.
data$PGS_std <- scale(data$PGS)


# ANOVA and non-parametric comparison
anova_result <- aov(PGS_std ~ IPI_group, data = data)
summary(anova_result)
kruskal.test(PGS_std ~ IPI_group, data = data)
TukeyHSD(anova_result)



# linear mod. with continuous IPI
model_ipi <- lm(PGS_std ~ IPI, data = data)
summary(model_ipi)

# linear models of PGS ~ IPI within each IPI_group
model_short  <- lm(PGS_std ~ IPI, data = filter(data, IPI_group == "Short"))
model_medium <- lm(PGS_std ~ IPI, data = filter(data, IPI_group == "Medium"))
model_long   <- lm(PGS_std ~ IPI, data = filter(data, IPI_group == "Long"))

summary(model_short)
summary(model_medium)
summary(model_long)

#plot grouped lm
ggplot(data, aes(x = IPI, y = PGS_std, color = IPI_group)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Effect of IPI on PGS within IPI groups",
    x = "IPI (days)",
    y = "PGS (standardised)",
    color = "IPI group"
  ) +
  theme_minimal()




# PGS distribution across IPI groups
ggplot(data, aes(x = IPI_group, y = PGS_std, fill = IPI_group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.2) +
  theme_minimal() +
  labs(title = "Distribution of PGS by IPI groups", x = "IPI group", y = "PGS (standardised)")

# PGS vs IPI (continuous)
ggplot(data, aes(x = IPI, y = PGS_std)) +
  geom_point(alpha = 0.3, color = "purple") +
  geom_smooth(method = "lm", se = TRUE, color = "black", size = 1.2) +
  theme_minimal() +
  labs(title = "Relationship between IPI and PGS", x = "IPI (days)", y = "PGS (standardised)")


ggsave("~/scratch/agnes/gene-IPI/results/plots/PGS_vs_IPI_continuous.png", width = 8, height = 6, dpi = 300)
ggsave("~/scratch/agnes/gene-IPI/results/plots/PGS_vs_IPI_boxplot.png", width = 8, height = 6, dpi = 300)

# 
theme_custom <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# interaction plot
ggplot(data, aes(x = PGS_std, y = SVLEN_UL_DG, color = IPI_group)) +
  geom_smooth(method = "lm", se = TRUE, size = 1.2) +
  scale_color_manual(values = c("Short" = "orange", "Medium" = "darkgreen", "Long" = "blue")) +
  labs(
    title = "Effect of PGS on gestational duration by IPI group",
    x = "Polygenic score (standardised)",
    y = "Gestational duration (days)",
    color = "IPI group"
  ) +
  theme_custom

ggsave("~/scratch/agnes/gene-IPI/results/plots/PGS_IPI_interaction_clean.png", width = 9, height = 6, dpi = 300)



# gestational duration x IPI + PGS
# Distribution of gestational duration
ggplot(data, aes(x = SVLEN_UL_DG)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of gestational Duration", x = "gestational Duration (days)", y = "Count")

# Distribution of IPI (days)
ggplot(data, aes(x = IPI)) +
  geom_histogram(bins = 50, fill = "darkgreen", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of IPI", x = "IPI (days)", y = "Count")

# Distribution of PGS
ggplot(data, aes(x = PGS_std)) +
  geom_histogram(bins = 50, fill = "purple", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of PGS (standardised)", x = "PGS", y = "Count")

# lm gestational duration x pgs x ipi
model_interaction <- lm(SVLEN_UL_DG ~ PGS_std * IPI_group, data = data)
summary(model_interaction)

# stratified
model_short <- lm(SVLEN_UL_DG ~ PGS_std, data = filter(data, IPI_group == "Short"))
model_medium <- lm(SVLEN_UL_DG ~ PGS_std, data = filter(data, IPI_group == "Medium"))
model_long <- lm(SVLEN_UL_DG ~ PGS_std, data = filter(data, IPI_group == "Long"))

summary(model_short)
summary(model_medium)
summary(model_long)

ggplot(data, aes(x = PGS_std, y = SVLEN_UL_DG, color = IPI_group)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(title = "Interaction: PGS × IPI group on gestational duration", x = "PGS (standardised)", y = "Gestational duration (days)")

####

