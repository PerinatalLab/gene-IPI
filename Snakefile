# libraries for data handling
import pandas as pd
import numpy as np
import json

# defining the two parity groups: nulliparous (first pregnancy) and multiparous
parity_names = ['nulliparous', 'multiparous']
rule all:
    input:
        # Phenotype files according to parity
        expand("results/phenotype/filtered_pregnancies_{parity}.csv", parity=parity_names),
        # Genotype files according to parity (st)
        expand("results/gwas/moba_common_qc_ipi_{parity}.pgen", parity=parity_names),
        expand("results/gwas/moba_common_qc_ipi_{parity}.pvar", parity=parity_names),
        expand("results/gwas/moba_common_qc_ipi_{parity}.psam", parity=parity_names),
        # Genomewide files for multiparous
        "results/gwas/moba_common_qc_ipi_multiparous_genomewide.pgen",
        "results/gwas/moba_common_qc_ipi_multiparous_genomewide.pvar",
        "results/gwas/moba_common_qc_ipi_multiparous_genomewide.psam",
        # Polygenic scores
        expand("results/pgs/ipi_pgs_{parity}.sscore", parity=parity_names),
        # CEU-qc phenotype + PGS 
        expand("results/phenotype/pheno_pgs_unique_{parity}.csv", parity=parity_names),
        # final phenotype + covariates        
#        expand("results/final_phenotype/IPI_pgs_covariates_{parity}.txt", parity=parity_names),
        # batch corrected phenotype and covariates
        "results/final_phenotype/IPI_pgs_covariates_multiparous_corrected.txt",
        # GxE genomewide output for multiparous only
        "results/gwas/gxe_ipi_gd_gw.SVLEN_UL_DG.glm.linear",
        "results/gwas/regenie/gxe_ipi_gd_gw_SVLEN_UL_DG.regenie"


# rule to clean phenotype data and create ID lists for genotype filtering
rule cleaned_data:
    input:
        "/mnt/scratch/agnes/PDB1724_MFR_541_v12.csv",
        "/mnt/scratch/agnes/parental_ID_to_PREG_ID.csv",
        "/mnt/work/p1724/v12/linkage_Mother_PDB1724.csv",
        "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common.psam"
    output:
        "results/phenotype/filtered_pregnancies_{parity}.csv",
        "results/phenotype/IDs_extract_{parity}.txt"
    script:
        "scripts/phenotype_stat.R"

# rule to perform SNP-level and sample-level quality control on genotype data
rule snp_qc:
    input:
        psam = "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common.psam",
        pvar = "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common.pvar",
        pgen = "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common.pgen",
        keep = "results/phenotype/IDs_extract_{parity}.txt",
        snplist = "SNP_to_extract.txt"
    output:
        pgen = "results/gwas/moba_common_qc_ipi_multiparous_genomewide_{parity}.pgen",
        pvar = "results/gwas/moba_common_qc_ipi_multiparous_genomewide_{parity}.pvar",
        psam = "results/gwas/moba_common_qc_ipi_multiparous_genomewide_{parity}.psam"
    params:
        pfile = "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common",
        results = "results/gwas/moba_common_qc_ipi_{parity}"
    log:
        "results/gwas/moba_common_qc_ipi_{parity}.log"
    shell:
        """
        plink2 \
          --pfile {params.pfile} \
          --keep {input.keep} \
          --extract {input.snplist} \
          --maf 0.01 \
          --geno 0.01 \
          --mind 0.01 \
          --hwe 1e-6 midp \
          --make-pgen \
          --out {params.results} > {log} 2>&1
        """

# rule to calculate polygenic scores for each sample
rule calculate_pgs:
    input:
        weights = "variant_weights.txt",
        pgen = "results/gwas/moba_common_qc_ipi_{parity}.pgen",
        pvar = "results/gwas/moba_common_qc_ipi_{parity}.pvar",
        psam = "results/gwas/moba_common_qc_ipi_{parity}.psam"
    output:
        sscore = "results/pgs/ipi_pgs_{parity}.sscore",
        log = "results/pgs/ipi_pgs_{parity}.log"
    params:
        pfile = "results/gwas/moba_common_qc_ipi_{parity}",
        out = "results/pgs/ipi_pgs_{parity}"
    shell:
        """
        mkdir -p results/pgs
        plink2 \
          --pfile {params.pfile} \
          --score {input.weights} 1 2 3 cols=scoresums \
          --nonfounders \
          --out {params.out} > {output.log} 2>&1
        """
# function to identify and remove related individuals based on kinship coefficients
def selectUnrelated(input_kin, df, x):
        kin= pd.read_csv(input_kin, header= 0, sep= '\t')
        kin= kin.loc[kin.Kinship > 0.125, :] # filtered out related pairs (kinship > ~1st cousins)
        kin= kin.loc[kin.ID1.isin(x.values), :]
        kin= kin.loc[kin.ID2.isin(x.values), :]
        kin= kin.loc[:, ['ID1','ID2','Kinship']]
        kin_temp= kin.copy()
        kin_temp.columns= ['ID2', 'ID1', 'Kinship']
        kin_temp= pd.concat([kin_temp, kin], ignore_index=True)
        kin_temp['n']= kin_temp.groupby('ID1')['ID1'].transform('count')
        kin_temp['nn']= kin_temp.groupby('ID2')['ID2'].transform('count')
        kin_temp.sort_values(by=['n', 'nn'], inplace= True)
        to_keep= list()
        for i in range(0, len(kin_temp.index)):
                if kin_temp.iloc[i, 0] in kin_temp.iloc[0:i, 1].values:
                        kin_temp.iloc[i, 1]= "X"
                else:
                        to_keep.append(kin_temp.iloc[i, 0])
        to_remove= [i for i in kin_temp.ID1 if i not in to_keep]
        to_remove= list(set(to_remove))
        remove= pd.DataFrame({'FID': to_remove})
        remove['IID']= remove.FID
        return remove

# rule to restrict analysis to CEU ancestry and merge with PGS
rule keep_CEU:
    input:
        pheno = "results/phenotype/filtered_pregnancies_{parity}.csv",
        pgs = "results/pgs/ipi_pgs_{parity}.sscore",
        ceu_ids = "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/pca/moba_genotypes_2024.12.03_ceu_core_ids"
    output:
        "results/phenotype/pheno_pgs_unique_{parity}.csv"
    script:
        "scripts/qc_ipi_pgs.R"

# rule to remove related individuals, merge phenotypes with PCA and batch info
rule remove_related:
    input:
        "results/phenotype/pheno_pgs_unique_{parity}.csv",
        "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/kinship/moba_genotypes_2024.12.03.kin0",
        "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/pca/moba_genotypes_2024.12.03.pcs",
        "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/batch/moba_genotypes_2024.12.03_batches",
        "resources/batches.json"
    output:
        "results/final_phenotype/IPI_pgs_covariates_{parity}.txt",
	"results/final_phenotype/IPI_GWAS_gxe_{parity}.txt"
    run:
        d = pd.read_csv(input[0])
        remove = selectUnrelated(input[1], d, d.SENTRIX_ID)
        d = d[~d.SENTRIX_ID.isin(remove.IID.values)]
        pcs = pd.read_csv(input[2], sep='\t')[['IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6']]
        d = pd.merge(d, pcs, right_on='IID', left_on='SENTRIX_ID')
        batches = pd.read_csv(input[3], sep='\t', header=0)
        batches.columns = ['IID', 'batch']
        d = pd.merge(d, batches, left_on='SENTRIX_ID', right_on='IID')
        with open(input[4], 'r') as fp:
            batches_dict = json.load(fp)
        d['batch'].replace(batches_dict, inplace=True)
        d['batch'] = d['batch'].str.replace(' ', '_')
        d.to_csv(output[0], sep='\t', index=False)
	d= d[['#FID', 'IID_x', 'SVLEN_UL_DG', 'IPI', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'batch', 'PARITET_5']]
	d.columns= ['FID', 'IID', 'SVLEN_UL_DG', 'IPI', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'batch', 'PARITET_5']
	d.to_csv(output[1], sep= '\t', header= True, index= False)

# rule to remove batch effect from IPI and create IPI_corrected file
rule correct_ipi_batch:
    input:
        "results/final_phenotype/IPI_pgs_covariates_multiparous.txt"
    output:
        "results/final_phenotype/IPI_pgs_covariates_multiparous_corrected.txt"
    script:
        "scripts/correct_ipi_batch.R"

# rule to generate figures summarising IPI and birth outcomes
rule figures_IPI:
    input:
        "/mnt/scratch/agnes/PDB1724_MFR_541_v12.csv",
        "/mnt/scratch/agnes/parental_ID_to_PREG_ID.csv",
        "/mnt/work/p1724/v12/linkage_Mother_PDB1724.csv",
        "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common.psam"
    output:
        "results/plots/phenotype/ipi_distribution_corrected.png",
        "results/plots/phenotype/ipi_density_IPI.png",
        "results/plots/phenotype/ipi_ushape.png",
        "results/plots/phenotype/ipi_delivery_type.png",
        "results/plots/phenotype/ipi_stillbirth_vs_livebirth.png"
    script:
        "scripts/figures_IPI_stat.R"

# rule to prepare genome-wide GxE dataset for multiparous pregnancies
rule prepare_gxe_genomewide_data:
    input:
        pfile = "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common.pgen",
        keep = "results/phenotype/IDs_extract_multiparous.txt"
    output:
        pgen = "results/gwas/moba_common_qc_ipi_multiparous_genomewide.pgen",
        pvar = "results/gwas/moba_common_qc_ipi_multiparous_genomewide.pvar",
        psam = "results/gwas/moba_common_qc_ipi_multiparous_genomewide.psam"
    log:
        "results/gwas/moba_common_qc_ipi_multiparous_genomewide.log"
    params:
        params_pfile = "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common"
    shell:
        """
        plink2 \
          --pfile {params.params_pfile} \
          --keep {input.keep} \
          --maf 0.01 \
          --geno 0.01 \
          --mind 0.01 \
          --hwe 1e-6 midp \
          --make-pgen \
          --out results/gwas/moba_common_qc_ipi_multiparous_genomewide > {log} 2>&1
        """

# rule to run genome-wide GxE interaction analysis using PLINK2
rule gxe_interaction_ipi_parameters12_gw:
    input:
        pgen = "results/gwas/moba_common_qc_ipi_multiparous_genomewide.pgen",
        pvar = "results/gwas/moba_common_qc_ipi_multiparous_genomewide.pvar",
        psam = "results/gwas/moba_common_qc_ipi_multiparous_genomewide.psam",
        keep = "results/phenotype/IDs_extract_multiparous.txt",
        pheno = "results/final_phenotype/IPI_GWAS_gxe_multiparous.txt",
        high_qual = "high_qual_snps.txt",
    output:
        "results/gwas/gxe_ipi_gd_gw.SVLEN_UL_DG.glm.linear",
        "results/gwas/gxe_ipi_gd_gw.log"
    params:
        pfile = "results/gwas/moba_common_qc_ipi_multiparous_genomewide",
        out = "results/gwas/gxe_ipi_gd_gw"
    threads: 10
    shell:
        """
          plink2 \
          --pfile {params.pfile} \
          --keep {input.keep} \
          --pheno {input.pheno} \
          --pheno-name SVLEN_UL_DG \
          --covar {input.pheno} \
          --covar-name IPI,PC1,PC2,PC3,PC4,PARITET_5 \
          --covar-variance-standardize \
          --glm interaction cols=+a1freq\
          --parameters 1-8 \
          --extract {input.high_qual} \
          --threads {threads} \
          --max-alleles 2 \
          --out {params.out}
        """

# rule to run genome-wide GxE interaction analysis using Regenie
rule gxe_interaction_ipi_parameters12_gw_regenie:
    input:
        pgen = "/mnt/scratch/moba/HDGB-MoBaGenetics/2024.12.03/geno/moba_genotypes_2024.12.03_common_no_multiallelic_joined.pgen",
        pvar = "/mnt/scratch/moba/HDGB-MoBaGenetics/2024.12.03/geno/moba_genotypes_2024.12.03_common_no_multiallelic_joined.pvar",
        psam = "/mnt/scratch/moba/HDGB-MoBaGenetics/2024.12.03/geno/moba_genotypes_2024.12.03_common_no_multiallelic_joined.psam",
        keep = "results/phenotype/IDs_extract_multiparous.txt",
        pheno = "results/final_phenotype/IPI_GWAS_gxe_multiparous.txt",
        high_qual = "high_qual_snps.txt",
    output:
        "results/gwas/regenie/gxe_ipi_gd_gw_SVLEN_UL_DG.regenie",
        "results/gwas/regenie/gxe_ipi_gd_gw.log"
    params:
        pfile = "/mnt/scratch/moba/HDGB-MoBaGenetics/2024.12.03/geno/moba_genotypes_2024.12.03_common_no_multiallelic_joined",
        out = "results/gwas/regenie/gxe_ipi_gd_gw"
    shell:
        """
        ./regenie_v4.1.gz_x86_64_Linux \
        --step 2 \
        --pgen {params.pfile} \
        --covarFile {input.pheno} \
        --phenoFile {input.pheno} \
        --keep {input.keep} \
        --phenoCol SVLEN_UL_DG \
        --interaction IPI \
        --covarColList IPI,PC1,PC2,PC3,PC4,PARITET_5,batch \
        --catCovarList batch \
        --extract {input.high_qual} \
        --bsize 1000 \
        --threads 20 \
        --ignore-pred \
        --out {params.out} \
        --verbose \
        --no-condtl
        """

# this sends a message to Agnes:: did she save the world or not?
onsuccess:
    shell("""
        curl -X POST -H 'Content-type: application/json' \
        --data '{{"text":"Hurray! Snakemake pipeline completed successfully!"}}' \
        https://hooks.slack.com/services/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
    """)

onerror:
    shell("""
        curl -X POST -H 'Content-type: application/json' \
        --data '{{"text":"Snakemake pipeline failed."}}' \
        https://hooks.slack.com/services/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
    """)
