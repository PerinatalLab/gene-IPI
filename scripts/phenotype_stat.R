# Quality control: iatrogen deliveries excluded from multiparous group, IPI definition, r-code for Snakemake 
library(data.table)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(tidyr)
library(broom)
library(usethis)


print(paste0('Parity = ', snakemake@wildcards[['parity']]))

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

# remove outlier gestational duration

mfr <- filter(mfr, SVLEN_UL_DG >= 154, SVLEN_UL_DG< 308)

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

if (snakemake@wildcards[['parity']] == 'multiparous') {
  
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
  
  cleaned_data <- cleaned_data %>% filter(FSTART == 1)

  
} else if (snakemake@wildcards[['parity']] == 'nulliparous') {
  
  merged_data <- merged_data %>%
    filter(PARITET_5 < 1)
  
  cleaned_data <- merged_data %>%
    filter(!is.na(SVLEN_UL_DG))
  
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
    filter(!is.na(DODFODTE_5), !is.na(SPABORT_12_5), !is.na(SPABORT_23_5))
  
}

cleaned_data <- cleaned_data[!is.na(cleaned_data$'#FID'), ]

write.csv(cleaned_data,snakemake@output[[1]], row.names = FALSE)

cleaned_data= cleaned_data[, c("#FID", "SENTRIX_ID")]
names(cleaned_data) = c("#FID", "IID")
fwrite(cleaned_data, snakemake@output[[2]], row.names = FALSE, sep= '\t')
