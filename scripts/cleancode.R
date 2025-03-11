#  IPI vs maternal health, IPI vs foetus health

library(data.table)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(tidyr)
library(broom)

library(usethis)

#setwd('/mnt/scratch/agnes/gene-IPI/')

#mfr <- fread("../PDB1724_MFR_541_v12.csv")
#mapping <- fread("/mnt/scratch/agnes/parental_ID_to_PREG_ID.csv")
mfr = fread(snakemake@input[[1]])
mapping = fread(snakemake@input[[2]])

merged_data <- mfr %>%
  inner_join(mapping, by = "PREG_ID_1724") %>%
  mutate(
    UL_TERMIN = as.Date(UL_TERMIN, format = "%Y-%m-%d"),
    ConceptionDate = UL_TERMIN - days(280),
    DeliveryDate = ConceptionDate + days(SVLEN_UL_DG)
  )

# duplicate pregnancies
duplicates <- merged_data %>%
  count(PREG_ID_1724) %>%
  filter(n > 1)

print(duplicates)

# twin pregnancies
twin_pregnancies <- merged_data %>%
  filter(FLERFODSEL == 1) 

merged_data <- merged_data  %>%
  filter(is.na(FLERFODSEL
  ))

# don't remove, keep max gestational duration
# keep one of the two randomly
twin_pregnancies = group_by(twin_pregnancies, PREG_ID_1724) %>% 
  filter(max(SVLEN_UL_DG) == SVLEN_UL_DG) %>% 
  sample_n(1)

merged_data <- bind_rows(merged_data, twin_pregnancies)

# to filter out pregnancies with parity ≥ 4
merged_data <- merged_data %>%
  filter(PARITET_5 < 4)


# to check for missing gestational age
# missing_gestational_age <- twin_pregnancies %>%
# filter(is.na(SVLEN_UL_DG))

# print(missing_gestational_age)

cleaned_data <- merged_data %>%
  filter(!is.na(SVLEN_UL_DG))
print(cleaned_data)



# CALCULATION of IPI
cleaned_data <- cleaned_data %>%
  arrange(M_ID_1724, DeliveryDate) %>%
  group_by(M_ID_1724) %>%
  mutate(
    IPI = as.numeric(difftime(DeliveryDate, lag(DeliveryDate), units = "days"))
  ) %>%
  ungroup()


# plot distribution of IPI
ggplot(cleaned_data, aes(x = IPI)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  labs(title = "Distribution of IPI (corrected)", x = "IPI (days)", y = "Count") +
  theme_minimal()

cleaned_data= mutate(cleaned_data, IPI = ifelse(IPI < 0, NA, IPI))  # to remove negative IPI values
summary(cleaned_data$IPI)

# filter out unrealistic IPI values (< 6 months)
cleaned_data <- cleaned_data %>%
  filter(!is.na(IPI), !is.na(UL_TERMIN), IPI >= 180) %>%
  filter(!is.na(PARITET_5))

#write.csv(filtered_pregnancies, "filtered_pregnancies.csv", row.names = FALSE)
#write.csv(filtered_pregnancies, snakmake@output[[1]])

# colnames(merged_data)

#print(filtered_pregnancies)

# select(IPI, UL_TERMIN, PARITET_5) %>% print()


# to check summary statistics
summary(cleaned_data$IPI)
table(is.na(cleaned_data$IPI))  # to ensure no missing IPI values
table(is.na(cleaned_data$PARITET_5))  # to ensure no missing parities



# VISUALISATION -- IPI

# IPI distribution
p1 = ggplot(cleaned_data, aes(x = IPI)) +
  geom_histogram(bins = 30, fill = "lightblue", color = "black") +
  labs(title = "Distribution of IPI (corrected)", x = "IPI (days)", y = "Count") +
  theme_minimal()


ggsave(snakemake@output[[2]], p1, width = 8, height = 6
)
mean_IPI <- mean(cleaned_data$IPI, na.rm = TRUE)
median_IPI <- median(cleaned_data$IPI, na.rm = TRUE)
cat("Mean IPI:", mean_IPI, "\n")
cat("Median IPI:", median_IPI, "\n")
# the mean IPI is greater than the median, the distribution is right (positively) skewed

# Density plot of IPI
p2 = ggplot(cleaned_data, aes(x = IPI)) +
  geom_density(fill = "lightgreen", alpha = 0.5) +
  labs(title = "Density plot of IPI", x = "IPI (days)", y = "Density") +
  theme_minimal()

# Boxplot of IPI by parity
p3 = ggplot(cleaned_data, aes(x = factor(PARITET_5), y = IPI)) +
  geom_boxplot(fill = "coral", alpha = 0.7) +
  labs(title = "IPI by parity", x = "Parity", y = "IPI (days)")

ggsave(snakemake@output[[3]], p2, width = 8, height = 6)
ggsave(snakemake@output[[4]], p3, width =8, height = 6)
table(cleaned_data$PARITET_5)
# summary(cleaned_data$IPI) # min> smallest IPI value, Q1> 25%, median Q2, 3rdQ 75, max> largest ipi value, excluding outliers

# U-shaped trend of preterm birth rates by IPI
# categorising IPI and proportion of preterm deliveries
cleaned_data <- cleaned_data %>%
  mutate(
    Preterm = ifelse(SVLEN_UL_DG < 259, 1, 0),
    IPI_range = cut(
      IPI,
      breaks = c(0, 180, 365, 730, 1460, Inf),
      labels = c("<6 months", "6–12 months", "1–2 years", "2–4 years", ">4 years"),
      include.lowest = TRUE
    )
  ) %>%
  group_by(IPI_range) %>%
  summarise(Proportion_Preterm = mean(Preterm, na.rm = TRUE))

# plot> U-shaped trend, relationship between the IPI range and proportion of preterm birth
# DUBUG> to exclude NAs from the plot
p4=ggplot(cleaned_data %>% filter(!is.na(IPI_range) & !is.na(Proportion_Preterm)), aes(x = IPI_range, y = Proportion_Preterm, group = 1)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(size = 3, color = "red") +
  labs(
    title = "U-shaped trend of preterm delivery rates by IPI range",
    x = "IPI Range",
    y = "Proportion of preterm deliveries"
  ) +
  theme_minimal()

#ggsave("u_shape_preterm_rates.png", width = 8, height = 6)
ggsave(snakemake@output[[5]],p4, width = 8, height = 6)



# GESTATIONAL DURATION, filter spontanious or not, preterm/term, plot by IPI


ipi_vs_preterm <- cleaned_data %>%
  group_by(IPI_Category, Birth_Category) %>%
  summarise(Count = n()) %>%
  ungroup()

P5=ggplot(ipi_vs_preterm, aes(x = IPI_Category, y = Count, fill = Birth_Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Preterm" = "red", "Term" = "lightblue")) + 
  labs(title = "Gestational duration by IPI", 
       x = "IPI category", 
       y = "Count", 
       fill = "Birth category") +
  theme_minimal()

#ggsave("gestational_duration_by_IPI_colored.png", width = 8, height = 6)
ggsave(snakemake@output[[6]], P5, width = 8, height = 6)

# delivery type by IPI 
ipi_vs_delivery_type <- cleaned_data %>%
  group_by(IPI_Category, Delivery_Type) %>%
  summarise(Count = n()) %>%
  ungroup()

P6=ggplot(ipi_vs_delivery_type, aes(x = IPI_Category, y = Count, fill = Delivery_Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Spontaneous" = "cyan4", "Induced" = "orange", "Cesarean" = "purple")) +
  labs(title = "Delivery type by IPI", 
       x = "IPI category", 
       y = "Count", 
       fill = "Delivery Type") +
  theme_minimal()

#ggsave("delivery_type_by_IPI.png", width = 8, height = 6)
ggsave(snakemake@output[[7]], P6, width = 8, height =6)

#Proportion of delivery type by IPI"
P7=ggplot(cleaned_data, aes(x = IPI_Category, fill = Delivery_Type)) +
  geom_bar(position = "fill") + 
  scale_fill_manual(values = c("Spontaneous" = "cyan4", "Induced" = "orange", "Cesarean" = "purple")) +
  labs(title = "Proportion of delivery type by IPI", 
       x = "IPI category", 
       y = "Proportion", 
       fill = "Delivery type") +
  theme_minimal()

#ggsave("induced_vs_spontaneous_by_IPI.png", width = 8, height = 6)
ggsave(snakemake@output[[8]], P7, width = 8, height = 6)

# Proportion of preterm vs. term births by IPI
P8=ggplot(cleaned_data, aes(x = IPI_Category, fill = Birth_Category)) +
  geom_bar(position = "fill") +  # Proportion-based stacked bar plot
  scale_fill_manual(values = c("Preterm" = "red", "Term" = "lightblue")) +
  labs(title = "Proportion of preterm vs. term births by IPI", 
       x = "IPI category", 
       y = "Proportion", 
       fill = "Birth category") +
  theme_minimal()

#ggsave("preterm_vs_term_by_IPI.png", width = 8, height = 6)
ggsave(sankemake@output[[9]], P8, width = 8, height = 6)

print(ipi_vs_delivery_type)

# Delivery type by IPI category (separated by preterm vs. term)
P9=ggplot(cleaned_data, aes(x = IPI_Category, fill = Delivery_Type)) +
  geom_bar(position = "fill") +
  facet_wrap(~ Birth_Category) +  
  scale_fill_manual(values = c("Spontaneous" = "cyan4", "Induced" = "orange", "Cesarean" = "purple")) +
  labs(title = "Delivery type by IPI category (separated by preterm vs. term)",
       x = "IPI category", 
       y = "Proportion", 
       fill = "Delivery type") +
  theme_minimal()

#ggsave("delivery_type_by_IPI_and_preterm.png", width = 10, height = 6)
ggsave(snakemake@output[[10]], P9, width = 8, height = 6)

# Distribution of IPI by the mode of delivery
P10=ggplot(cleaned_data, aes(x = Delivery_Type, y = IPI, fill = Delivery_Type)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Spontaneous" = "cyan4", "Induced" = "orange", "Cesarean" = "purple")) +
  labs(title = "Distribution of IPI by the mode of delivery",
       x = "Delivery type", 
       y = "IPI (days)") +
  theme_minimal()

#ggsave("IPI_by_delivery_type.png", width = 8, height = 6)
ggsave(snakemake@output[[11]], P10, width = 8, height = 6)



cleaned_data <- cleaned_data %>%
  mutate(
    IPI_Category = cut(
      IPI,
      breaks = c(180, 365, 730, 1460, Inf),
     labels = c("6-12 months", "1-2 years", "2-4 years", ">4 years"),
      include.lowest = TRUE
    )
  )

# IPI vs c-section
cleaned_data <- cleaned_data %>%
  mutate(
    Delivery_Type = case_when(
      FSTART == 1 ~ "Spontaneous",
      FSTART == 2 ~ "Induced",
      FSTART == 3 ~ "Cesarean",
      TRUE ~ NA_character_
    )
  )

cleaned_data <- cleaned_data %>%
  mutate(
    IPI_Category = cut(
      IPI,
      breaks = c(180, 365, 730, 1460, Inf),  
      labels = c("6-12 months", "1-2 years", "2-4 years", ">4 years"),
      include.lowest = TRUE
    )
  )

ipi_vs_csection <- cleaned_data %>%
  filter(!is.na(Delivery_Type)) %>%  # Ensure no missing delivery types
  mutate(C_Section = ifelse(Delivery_Type == "Cesarean", 1, 0)) %>%
  group_by(IPI_Category) %>%
  summarise(Proportion_C_Section = mean(C_Section, na.rm = TRUE))

P11=ggplot(ipi_vs_csection, aes(x = IPI_Category, y = Proportion_C_Section, fill = IPI_Category)) +
  geom_col() +
  labs(title = "C-section rate by IPI",
       x = "IPI category",
       y = "Proportion of C-section deliveries") +
  scale_fill_manual(values = c("6-12 months" = "coral1",
                               "1-2 years" = "coral1",
                               "2-4 years" = "coral1",
                               ">4 years" = "coral1")) +
  theme_minimal()

#ggsave("c_section_by_IPI_category.png", width = 8, height = 6)
ggsave(snakemake@output[[12]], P11, width = 8, height =6)

print(ipi_vs_csection)


# statistical analysis

# Chi-square test for IPI vs. preterm birth
cleaned_data <- cleaned_data %>%
  mutate(
    Birth_Category = case_when(
      SVLEN_UL_DG < 259 ~ "Preterm",
      SVLEN_UL_DG >= 259 ~ "Term",
      TRUE ~ NA_character_
    )
  )

ipi_vs_preterm_table <- table(cleaned_data$IPI_Category, cleaned_data$Birth_Category)
chisq_test_preterm <- chisq.test(ipi_vs_preterm_table)
print(chisq_test_preterm)

prop.table(ipi_vs_preterm_table, margin = 1)  # proportion of preterm births in each IPI category

kruskal_test_preterm <- kruskal.test(IPI ~ Birth_Category, data = cleaned_data)
print(kruskal_test_preterm)

# Chi-square test for IPI vs. delivery type
ipi_vs_delivery_table <- table(cleaned_data$IPI_Category, cleaned_data$Delivery_Type)
chisq_test_delivery <- chisq.test(ipi_vs_delivery_table)
print(chisq_test_delivery)


prop.table(ipi_vs_delivery_table, margin = 1)




# analysis of stillbirth vs live births by IPI

colnames(merged_data)[grepl("DODFODTE|LEVENDEFODTE", colnames(merged_data))]



# binary coding for birth outcome (stillbirth vs. live birth)
cleaned_data <- cleaned_data %>%
  mutate(
    Birth_Outcome = case_when(
      DODFODTE_5 >= 1 ~ "Stillbirth",
      LEVENDEFODTE_5 >= 1 ~ "Live birth",
      TRUE ~ NA_character_
    )
  )

# frequency of stillbirth vs. live birth
table(cleaned_data$Birth_Outcome, useNA = "always")

# contingency table for stillbirth vs. live birth by IPI
ipi_vs_birth_table <- table(cleaned_data$IPI_Category, cleaned_data$Birth_Outcome)
print(ipi_vs_birth_table)

# proportions
prop_table_birth <- prop.table(ipi_vs_birth_table, margin = 1)
print(prop_table_birth)

# Chi-square test for association between IPI and birth outcome
chisq_test_birth <- chisq.test(ipi_vs_birth_table)
print(chisq_test_birth)

# Cramér's V (Effect Size)
chisq_stat <- chisq_test_birth$statistic  # chi-square value
n_total <- sum(ipi_vs_birth_table)  # sample size
r <- nrow(ipi_vs_birth_table)  # no. of IPI categories
c <- ncol(ipi_vs_birth_table)  # no. of birth outcomes

# Cramér’s V calculation
cramers_v <- sqrt(chisq_stat / (n_total * min(r - 1, c - 1)))
cat("Cramér's V:", cramers_v, "\n")


# Visualisation: proportion of stillbirths vs. live births by IPI
P12=ggplot(cleaned_data, aes(x = IPI_Category, fill = Birth_Outcome)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("Stillbirth" = "red", "Live birth" = "blue")) +
  labs(title = "Proportion of stillbirths vs. live births by IPI",
       x = "IPI category",
       y = "Proportion",
       fill = "Birth outcome") +
  theme_minimal()

#ggsave("stillbirth_vs_livebirth_by_IPI.png", width = 8, height = 6)
ggsave(snakemake@output[[13]], P12, width = 8, height = 6)
write.csv(cleaned_data, snakemake@output[[1]])
