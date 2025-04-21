library(data.table)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(broom)

setwd('/mnt/scratch/agnes/gene-IPI/')

mfr <- fread("../PDB1724_MFR_541_v12.csv")
#%>%  select(PREG_ID_1724, SVLEN_UL_DG, UL_TERMIN, PARITET_5, FAAR)
mapping <- fread("/mnt/scratch/agnes/parental_ID_to_PREG_ID.csv")

merged_data <- mfr %>%
  inner_join(mapping, by = "PREG_ID_1724") %>%
  mutate(
    UL_TERMIN = as.Date(UL_TERMIN, format = "%Y-%m-%d"),
    ConceptionDate = UL_TERMIN - days(280),
    DeliveryDate = ConceptionDate + days(SVLEN_UL_DG)
  )

#merged data check
colnames(merged_data) # megnézni milyen változók vannak még ami fontos lehet, akár későbbi elemzésnél
dim(mfr)
dim(mapping)
dim(merged_data)

#checking for duplicates
merged_data %>%
  count(PREG_ID_1724) %>%
  filter(n > 1)

#keep one pregnancy from twin pregnancy
merged_data = group_by(merged_data, PREG_ID_1724, !is.na(FLERFODSEL)) %>% sample_n(1)

# Exclude pregnancies with parity 4 or more
merged_data <- filter(merged_data, PARITET_5 < 4)

#IPI calculation
merged_data <- merged_data %>%
  arrange(M_ID_1724, DeliveryDate) %>%
  group_by(M_ID_1724) %>%
  mutate(
    IPI = ifelse(PARITET_5 > 0, as.numeric(difftime(ConceptionDate, lag(DeliveryDate), units = "days")), NA)
  ) %>%
  ungroup

#QC
#date calculation check
head(merged_data %>% select(UL_TERMIN, ConceptionDate, DeliveryDate))

#looking for iff diff is non-zero for any rows
merged_data %>%
  mutate(diff = as.numeric(difftime(UL_TERMIN, DeliveryDate, units = "days"))) %>%
  filter(diff != 0)

#IPI // is there any neg value, nulliparous (PARITET_5 == 0) IPI = NA
summary(merged_data$IPI)
table(merged_data$IPI<0)

#xxchecking edge cases, if no rows = ok
merged_data %>%
  filter(IPI < 0 | PARITET_5 == 0 & !is.na(IPI)) %>% select(UL_TERMIN, ConceptionDate, DeliveryDate,
                                                            FAAR, M_ID_1724, PREG_ID_1724, IPI)

merged_data %>% filter(M_ID_1724 == "M052463") %>% select(UL_TERMIN, ConceptionDate, DeliveryDate,
                                                          FAAR, PREG_ID_1724, M_ID_1724, PARITET_5, IPI)
#exclude every neg. IPI 
merged_data <- filter(merged_data, IPI>0 | is.na(IPI))

#Invalid parity progression
invalid_parity <- merged_data %>%
  group_by(M_ID_1724) %>%
  filter(PARITET_5 < lag(PARITET_5, default = 0))

select(invalid_parity, UL_TERMIN, ConceptionDate, DeliveryDate,
       FAAR, PREG_ID_1724, M_ID_1724, PARITET_5, IPI)

merged_data %>% filter(M_ID_1724 == "M000554") %>% select(UL_TERMIN, ConceptionDate, DeliveryDate,
                                                          FAAR, PREG_ID_1724, M_ID_1724, PARITET_5, IPI)


invalid_parity %>%
  count(M_ID_1724)


#skipped parities
skipped_parity <- merged_data %>%
  group_by(M_ID_1724) %>%
  filter(PARITET_5 - lag(PARITET_5) != 1) %>% mutate(diff_parity= PARITET_5 - lag(PARITET_5))

# Make a list of ids where there is a jump in parity (backwards or front)
bad_ids= skipped_parity %>% 
  select(UL_TERMIN, ConceptionDate, DeliveryDate,
         FAAR, PREG_ID_1724, M_ID_1724, PARITET_5, IPI, diff_parity) %>% 
  filter(!is.na(diff_parity)) %>%
  pull(M_ID_1724)

# Exclude mothers with jumps in parity
merged_data <- filter(merged_data, !(M_ID_1724 %in% bad_ids))


skipped_parity %>%
  count(M_ID_1724)

# We remove mother with parity = 0
merged_data <- filter(merged_data, PARITET_5>0)


#remove all women who only have one pregnancy
merged_data <- filter(merged_data, !is.na (IPI))

#twin pregnancies
twin_pregnancies <- merged_data %>%
  filter(IPI < 60)

# ! = number of twin pregnancies
twin_pregnancies %>%
  count(PARITET_5)

#negative IPI and implausible gestational durations are excluded
merged_data <- merged_data %>%
  filter(SVLEN_UL_DG >= 154 & SVLEN_UL_DG <= 310, IPI >= 0)

#IPI ranges and preterm births ident
u_shape_data <- merged_data %>%
  mutate(
    IPI_range = cut(
      IPI,
      breaks = c(0, 180, 365, 730, 1460, Inf),
      labels = c("<6 months", "6–12 months", "1–2 years", "2–4 years", ">4 years"),
      include.lowest = TRUE
    ),
    Preterm = ifelse(SVLEN_UL_DG < 259, "Preterm", "Term")
  ) %>%
  group_by(IPI_range, Preterm) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(IPI_range) %>%
  mutate(Proportion = Count / sum(Count))

write.csv(merged_data, "mfr_results_with_IPI.csv", row.names = FALSE)
write.csv(invalid_parity, "invalid_parity_cases.csv", row.names = FALSE)
write.csv(skipped_parity, "skipped_parity_cases.csv", row.names = FALSE)
write.csv(twin_pregnancies, "twin_pregnancies.csv", row.names = FALSE)

cat("Invalid parity cases:", nrow(invalid_parity), "\n")
cat("Skipped parity cases:", nrow(skipped_parity), "\n")
cat("Twin pregnancy cases (IPI < 60):", nrow(twin_pregnancies), "\n")

#Visualisations
#U-shaped trend of preterm birth rates
u_shape_data %>%
  filter(Preterm == "Preterm") %>%
  ggplot(aes(x = IPI_range, y = Proportion, group = 1)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(size = 3, color = "red") +
  labs(
    title = "U-shaped trend of preterm delivery rates by IPI range",
    x = "IPI range",
    y = "Proportion of preterm deliveries"
  ) +
  theme_minimal()

ggsave("u_shape_preterm_rates.png", width = 8, height = 6)

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

ggsave("ipi_distribution.png", width = 8, height = 6)

#logistic regression
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

odds_ratios <- tidy(logistic_model, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = gsub("IPI_range", "", term))

#forest plot
odds_ratios %>%
  ggplot(aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange(color = "blue", size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    title = "Odds ratios for preterm delivery by IPI range",
    x = "IPI range",
    y = "Odds ratio (log scale)"
  ) +
  scale_y_log10() +
  theme_minimal()

ggsave("forest_plot_odds_ratios.png", width = 8, height = 6)

#principal component analysis (PCA)
#pca_data <- merged_data %>%
#select(SVLEN_UL_DG, IPI, PARITET_5) %>%
#drop_na()  #no missing values -- or

#pca_results <- prcomp(pca_data, scale. = TRUE)

#screeplot(pca_results, type = "lines", main = "Scree plot of PCA")

#biplot(pca_results, main = "PCA biplot")

#autoplot(pca_results, data = pca_data, colour = 'PARITET_5') +
# labs(title = "PCA of key variables")

#ggsave("pca_plot.png", width = 8, height = 6)

#PCA using ggplot2 only
pca_data <- merged_data %>%
  select(SVLEN_UL_DG, IPI, PARITET_5) %>%
  drop_na()  #missing values<

pca_results <- prcomp(pca_data, scale. = TRUE)

summary(pca_results)

#extr PCA scores / principal components
pca_scores <- as.data.frame(pca_results$x) %>%
  mutate(Individual = rownames(.))

#print PC1 and PC2
ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.7, color = "blue") +
  labs(
    title = "PCA of key variables (PC1 vs PC2)",
    x = "Principal component 1 (PC1)",
    y = "Principal component 2 (PC2)"
  ) +
  theme_minimal()
