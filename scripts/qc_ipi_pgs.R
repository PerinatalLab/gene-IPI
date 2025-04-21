# merging pgs and phenotype dataset, one ID/mother

library(data.table)
library(ggplot2)
library(dplyr)

setwd('/mnt/scratch/agnes/gene-IPI/')

#pheno <- fread("results/phenotype/filtered_pregnancies.csv")
pheno <- fread(snakemake@input[[1]])
#pgs <- fread("results/pgs/ipi_pgs.sscore")
pgs <- fread(snakemake@input[[2]])
names(pgs)=c("IID", "PGS")

# ceu ancestor kept only
#ceu <- fread("/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/pca/moba_genotypes_2024.12.03_ceu_core_ids", header = FALSE)
ceu <- fread(snakemake@input[[3]])
names(ceu)= c("FID", "IID")

pheno_pgs=inner_join(pheno, pgs, by= c("SENTRIX_ID"="IID"))
pheno_pgs=filter(pheno_pgs, SENTRIX_ID %in% ceu$IID)

pheno_pgs_unique <- pheno_pgs %>%
  group_by(SENTRIX_ID) %>%
  slice(1) %>%
  ungroup()


#fwrite(pheno_pgs_unique, "results/phenotype/pheno_pgs_unique.csv")
fwrite(pheno_pgs_unique, snakemake@output[[1]])


#check 
#head(pheno_pgs_unique)

# check duplicates
#any(duplicated(pheno_pgs_unique$SENTRIX_ID))

# check a random mother_ID
#mother_id <- pheno_pgs_unique %>%
 # filter(PREG_ID_1724 == 105579) %>%
  #pull(M_ID_1724)

#pheno_pgs %>%
 # filter(M_ID_1724 == mother_id) %>%
  #select(PREG_ID_1724, M_ID_1724, SENTRIX_ID)

#table(table(pheno_pgs_unique$SENTRIX_ID))




       
