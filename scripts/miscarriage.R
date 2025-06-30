# 

library(data.table)
library(dplyr)
library(broom)
library(tidyr)
library(ggplot2)

setwd("/mnt/scratch/agnes/gene-IPI/")

nulli <- fread("results/final_phenotype/IPI_pgs_covariates_nulliparous.txt")
nulli= filter(nulli, FSTART== 1)

# 
nulli$miscar <- ifelse(nulli$SPABORT_12_5 > 0 | nulli$SPABORT_23_5 > 0 | nulli$DODFODTE_5 > 0 ,"previous miscarriage", "no prev misc")
table(nulli$miscar)

prop.table(table(nulli$miscar))



# PGS standars.
nulli$PGS_std <- scale(nulli$PGS)


# alap stat teszteket lefuttatni, ld jegyzet

# ANOVA and non-parametric comparison
anova_result <- aov(PGS_std ~ nulli$miscar, data = nulli)
summary(anova_result)



# linear mod. with continuous miscar
model_miscar <- lm(PGS_std ~ miscar, data = nulli)
summary(model_miscar)

model2 = lm(SVLEN_UL_DG ~ PGS_std + miscar, data = nulli)
summary(model2)

model3 = lm (SVLEN_UL_DG ~ PGS_std  * miscar, data = nulli)
summary(model3)

# interaction plot
ggplot(nulli, aes(x = PGS_std, y = SVLEN_UL_DG, color = miscar)) +
  geom_smooth(method = "lm", se = TRUE, size = 1.2) +
  scale_color_manual(values = c("blue", "green")) +
  labs(
    title = "Effect of PGS on gestational duration by miscarriage",
    x = "Polygenic score (standardised)",
    y = "Gestational duration (days)",
    color = "Miscarriage"
  ) +
  theme_minimal() 

ggsave("~/scratch/agnes/gene-IPI/results/plots/miscar_pgs.png", width = 9, height = 6, dpi = 300)

