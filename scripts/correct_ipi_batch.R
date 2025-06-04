library(data.table)

input_file <- "/home/a24agnju/scratch/agnes/gene-IPI/results/final_phenotype/IPI_pgs_covariates_multiparous.txt"
output_file <- "/home/a24agnju/scratch/agnes/gene-IPI/results/final_phenotype/IPI_pgs_covariates_multiparous_corrected.txt"
pheno <- fread(input_file)

# IPI batch correction
model <- lm(IPI ~ batch, data = pheno)
pheno$IPI_corrected <- resid(model)

fwrite(pheno, output_file, sep = "\t")
