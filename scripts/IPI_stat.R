# Quality control and plots: IPI vs foetus and maternal health, r-code for Snakemake 
library(data.table)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(tidyr)
library(broom)
library(usethis)


setwd('/mnt/scratch/agnes/gene-IPI/')

#mfr <- fread("../PDB1724_MFR_541_v12.csv")
mfr <- fread(snakemake@input[[1]])
#mapping <- fread("/mnt/scratch/agnes/parental_ID_to_PREG_ID.csv")
mapping <- fread(snakemake@input[[2]])
linkage <- fread(snakemake@input[[3]])
psam <- fread(snakemake@input[[4]])
psam= psam[ ,c("#FID", "IID" )]
#output for snakemake: filtered pregnancies
filtered_pregnancies <- (snakemake@output[[1]])
#dir.create("results/phenotype", recursive = TRUE, showWarnings = FALSE)
#write.csv(mfr, "results/phenotype/filtered_pregnancies.csv", row.names = FALSE)

# maternal records with mapping data
merged_data <- mfr %>%
  inner_join(mapping, by = "PREG_ID_1724") %>%
  mutate(
    UL_TERMIN = as.Date(UL_TERMIN, format = "%Y-%m-%d"),
    ConceptionDate = UL_TERMIN - days(280),
    DeliveryDate = ConceptionDate + days(SVLEN_UL_DG)
  )

#merge M ID to samples ID
merged_data <- merged_data %>%
	left_join(linkage, by = "M_ID_1724")

#merged FID and sample ID
merged_data <- merged_data %>%
        left_join(psam, by = c("SENTRIX_ID"="IID"))

# print duplicate pregnancies
duplicates <- merged_data %>%
  count(PREG_ID_1724) %>%
  filter(n > 1)
print(duplicates)

# twin pregnancies
twin_pregnancies <- merged_data %>%
  filter(FLERFODSEL == 1) 

# twin with the maximum gestational duration kept
twin_pregnancies <- twin_pregnancies %>%
  group_by(PREG_ID_1724) %>%
  filter(max(SVLEN_UL_DG) == SVLEN_UL_DG) %>%
  sample_n(1)

# removing twins from the original data
merged_data <- merged_data %>%
  filter(is.na(FLERFODSEL))

# selected twin pregnancies back
merged_data <- bind_rows(merged_data, twin_pregnancies)

# filtering out pregnancies with parity ≥ 4
merged_data <- merged_data %>%
  filter(PARITET_5 < 4)

# removing rows with missing gestational age
cleaned_data <- merged_data %>%
  filter(!is.na(SVLEN_UL_DG))
print(cleaned_data)

# IPI calculation 
cleaned_data <- cleaned_data %>%
  arrange(M_ID_1724, DeliveryDate) %>%
  group_by(M_ID_1724) %>%
  mutate(
    IPI = as.numeric(difftime(DeliveryDate, lag(DeliveryDate), units = "days"))
  ) %>%
  ungroup()

# removing negative IPI values
cleaned_data <- cleaned_data %>%
  mutate(IPI = ifelse(IPI < 0, NA, IPI))

# filtering out unrealistic IPI values (< 6 months)
cleaned_data <- cleaned_data %>%
  filter(!is.na(IPI), IPI >= 180)

summary(cleaned_data$IPI)

# IPI_category for plots
cleaned_data <- cleaned_data %>%
  mutate(
    Preterm = ifelse(SVLEN_UL_DG < 259, 1, 0),
    IPI_Category = cut(
      IPI,
      breaks = c(180, 365, 730, 1460, Inf),
      labels = c("6-12 months", "1-2 years", "2-4 years", ">4 years"),
      include.lowest = TRUE
    )
  )

# Delivery_type category
cleaned_data <- cleaned_data %>%
  mutate(
    Delivery_Type = case_when(
      FSTART == 1 ~ "Spontaneous",
      FSTART == 2 ~ "Induced",
      FSTART == 3 ~ "Cesarean",
      TRUE ~ NA_character_
    )
  )

# IPI visualisation
# histogram of IPI distribution
p1 = ggplot(cleaned_data, aes(x = IPI)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  labs(title = "Distribution of IPI", x = "IPI (days)", y = "Count") +
  theme_minimal()

#ggsave("ipi_distribution_corrected.png", width = 8, height = 6)
ggsave(snakemake@output[[2]], p1, width = 8, height = 6)

# Density plot of IPI
p2 = ggplot(cleaned_data, aes(x = IPI)) +
  geom_density(fill = "lightgreen", alpha = 0.5) +
  labs(title = "Density plot of IPI", x = "IPI (days)", y = "Density") +
  theme_minimal()

#ggsave("density_ipi.png", width = 8, height = 6) 
ggsave(snakemake@output[[3]], p2, width = 8, height = 6)

# U-shape trend
ipi_vs_preterm <- cleaned_data %>%
  group_by(IPI_Category) %>%
  summarise(Proportion_Preterm = mean(Preterm, na.rm = TRUE))

p3 = ggplot(ipi_vs_preterm, aes(x = IPI_Category, y = Proportion_Preterm, group = 1)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(size = 3, color = "red") +
  labs(
    title = "U-Shaped trend of preterm birth by IPI",
    x = "IPI Category",
    y = "Proportion of preterm births"
  ) +
  theme_minimal()

#ggsave("u_shape_preterm_rates.png", width = 8, height = 6)
ggsave(snakemake@output[[4]], p3, width = 8, height = 6)

# delivery type vs IPI
ipi_vs_delivery <- cleaned_data %>%
  group_by(IPI_Category, Delivery_Type) %>%
  summarise(Count = n(), .groups = "drop")

p4 = ggplot(ipi_vs_delivery, aes(x = IPI_Category, y = Count, fill = Delivery_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Spontaneous" = "cyan4", "Induced" = "orange", "Cesarean" = "purple")) +
  labs(title = "Delivery type by IPI", x = "IPI category", y = "Count") +
  theme_minimal()

#ggsave("delivery_type_by_IPI.png", width = 8, height = 6)
ggsave(snakemake@output[[5]], p4, width = 8, height = 6)

# stillbirth vs live birth by IPI
# birth outcomes categories
cleaned_data <- cleaned_data %>%
  mutate(
    Birth_Outcome = case_when(
      !is.na(DODFODTE_5) & DODFODTE_5 >= 1 ~ "Stillbirth",
      !is.na(LEVENDEFODTE_5) & LEVENDEFODTE_5 >= 1 ~ "Live birth",
      TRUE ~ NA_character_
    )
  )

# removing missing birth outcomes for plotting
cleaned_data2 <- cleaned_data %>%
  filter(!is.na(Birth_Outcome))

p5 = ggplot(cleaned_data2, aes(x = IPI_Category, fill = Birth_Outcome)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("Stillbirth" = "red", "Live birth" = "blue")) +
  labs(title = "Stillbirth vs. live birth by IPI", x = "IPI category", y = "Proportion") +
  theme_minimal()

#ggsave("stillbirth_vs_livebirth_by_IPI.png", width = 8, height = 6)
ggsave(snakemake@output[[6]], p5, width = 8, height = 6)

write.csv(cleaned_data,snakemake@output[[1]], row.names = FALSE)

cleaned_data= cleaned_data[, c("#FID", "SENTRIX_ID")]
names(cleaned_data) = c("#FID", "IID")
fwrite(cleaned_data, snakemake@output[[7]], row.names = FALSE, sep= '\t')
