# Run steps until ###### every session
#cleaned, QC, U-shape, review, ? validation ?
  
library(data.table)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(broom)

setwd('/mnt/scratch/agnes/gene-IPI/')

mfr <- fread("../PDB1724_MFR_541_v12.csv") %>%
  select(PREG_ID_1724, SVLEN_UL_DG, UL_TERMIN, PARITET_5, FAAR)

mapping <- fread("/mnt/scratch/agnes/parental_ID_to_PREG_ID.csv")

merged_data <- mfr %>%
  inner_join(mapping, by = "PREG_ID_1724") %>%
  mutate(
    UL_TERMIN = as.Date(UL_TERMIN, format = "%Y-%m-%d"),
    ConceptionDate = UL_TERMIN - days(280),
    DeliveryDate = ConceptionDate + days(SVLEN_UL_DG)
  )

# IPI - calculation
merged_data <- merged_data %>%
  arrange(M_ID_1724, DeliveryDate) %>%
  group_by(M_ID_1724) %>%
  mutate(
    IPI = as.numeric(difftime(ConceptionDate, lag(DeliveryDate), units = "days"))  #IPI in days
  ) %>%
  ungroup()


###QC: sanity checks
#invalid parity progression
invalid_parity <- merged_data %>%
  group_by(M_ID_1724) %>%
  filter(PARITET_5 < lag(PARITET_5, default = 0))

#skipped parities
skipped_parity <- merged_data %>%
  group_by(M_ID_1724) %>%
  filter(PARITET_5 - lag(PARITET_5, default = PARITET_5[1]) > 1)

#twin pregnancies (IPI < 60 days)
twin_pregnancies <- merged_data %>%
  filter(IPI < 60)

#excl invalid rows: neg IPI or implausible gestational duration
merged_data <- merged_data %>%
  filter(SVLEN_UL_DG >= 154 & SVLEN_UL_DG <= 310, IPI >= 0)

###proportion of preterm deliveries by IPI range
# IPI ranges and calculate preterm rates
u_shape_data <- merged_data %>%
  mutate(
    IPI_range = cut(
      IPI,
      breaks = c(0, 180, 365, 730, 1460, Inf),
      labels = c("<6 months", "6–12 months", "1–2 years", "2–4 years", ">4 years"),
      include.lowest = TRUE
    ),
    Preterm = ifelse(SVLEN_UL_DG < 259, "Preterm", "Term")  #preterm defined as <259 days
  ) %>%
  group_by(IPI_range, Preterm) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(IPI_range) %>%
  mutate(Proportion = Count / sum(Count))

write.csv(merged_data, "mfr_results_with_IPI.csv", row.names = FALSE)

# flagged cases for review
write.csv(invalid_parity, "invalid_parity_cases.csv", row.names = FALSE)
write.csv(skipped_parity, "skipped_parity_cases.csv", row.names = FALSE)
write.csv(twin_pregnancies, "twin_pregnancies.csv", row.names = FALSE)

cat("Processed data and flagged cases have been saved.\n")

### Visualisations
# U-shaped trend of preterm delivery rates
u_shape_data %>%
  filter(Preterm == "Preterm") %>%
  ggplot(aes(x = IPI_range, y = Proportion, group = 1)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(size = 3, color = "red") +
  labs(
    title = "U-Shaped trend of preterm delivery rates by IPI range",
    x = "IPI range",
    y = "Proportion of preterm deliveries"
  ) +
  theme_minimal()

#distribution of IPI
merged_data %>%
  ggplot(aes(x = IPI)) +
  geom_histogram(binwidth = 30, fill = "red", color = "black") +
  labs(
    title = "Distribution of interpregnancy interval (IPI)",
    x = "IPI (days)",
    y = "Frequency"
  ) +
  theme_minimal()


###review and validation
# flagged cases
View(invalid_parity)
View(skipped_parity)
View(twin_pregnancies)

#flagged cases
cat("Invalid Parity Cases:", nrow(invalid_parity), "\n")
cat("Skipped Parity Cases:", nrow(skipped_parity), "\n")
cat("Twin Pregnancy Cases (IPI < 60):", nrow(twin_pregnancies), "\n")

#distribution of twin pregnancies by parity
twin_pregnancies %>%
  count(PARITET_5) %>%
  arrange(desc(n))

#U-shape plot
ggsave("u_shape_preterm_rates.png", width = 8, height = 6)

#logistic regression
#preterm and IPI ranges
merged_data <- merged_data %>%
  mutate(
    Preterm = ifelse(SVLEN_UL_DG < 259, 1, 0),
    IPI_range = cut(
      IPI,
      breaks = c(0, 180, 365, 730, 1460, Inf),
      labels = c("<6 months", "6–12 months", "1–2 years", "2–4 years", ">4 years"),
      include.lowest = TRUE
    )
  )

logistic_model <- glm(
  Preterm ~ IPI_range,
  data = merged_data,
  family = binomial(link = "logit")
)

#summary
summary(logistic_model)

#visualisation
#extract odds ratios and confidence intervals
odds_ratios <- tidy(logistic_model, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = gsub("IPI_range", "", term))

# Forest plot
odds_ratios %>%
  ggplot(aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange(color = "blue", size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  #OR = 1 reference line
  coord_flip() +
  labs(
    title = "Odds ratios for preterm delivery by IPI range",
    x = "IPI range",
    y = "Odds ratio (log scale)"
  ) +
  scale_y_log10() +
  theme_minimal()

ggsave("forest_plot_odds_ratios.png", width = 8, height = 6)
