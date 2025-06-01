# libraries for data handling
import pandas as pd
import numpy as np
import json

# defining the two parity groups: nulliparous (first pregnancy) and multiparous
parity_names= ['nulliparous', 'multiparous']

# the final target files that the pipeline should produce
rule all:
    input:
        # final phenotype files per parity
        expand("results/phenotype/filtered_pregnancies_{parity}.csv", parity= parity_names),
        # GWAS genotype data (after QC)
        "results/gwas/moba_common_qc_ipi.pgen",
        "results/gwas/moba_common_qc_ipi.pvar",
        "results/gwas/moba_common_qc_ipi.psam",
        # polygenic scores
        expand("results/pgs/ipi_pgs_{parity}.sscore", parity= parity_names),
        # phenotype merged with PGS per parity
        expand("results/phenotype/pheno_pgs_unique_{parity}.csv", parity= parity_names),
        # final analysis-ready phenotype with covariates per parity
        expand("results/final_phenotype/IPI_pgs_covariates_{parity}.txt", parity= parity_names)

# rule to clean phenotype data and create ID lists for genotype filtering
rule cleaned_data:
    input:
        # raw pregnancy data, parental linkage data, mother-child linkage file, and genotype sample info
        "/mnt/scratch/agnes/PDB1724_MFR_541_v12.csv",
        "/mnt/scratch/agnes/parental_ID_to_PREG_ID.csv",
        "/mnt/work/p1724/v12/linkage_Mother_PDB1724.csv",
        "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common.psam" 
    output:
        # filtered phenotype per parity and sample ID files
        "results/phenotype/filtered_pregnancies_{parity}.csv",
        "results/phenotype/IDs_extract_{parity}.txt"
    script:
         # r script to clean and filter phenotypic data
        "scripts/phenotype_stat.R"

# rule to perform SNP-level and sample-level quality control on genotype data
rule snp_qc:
    input:
        # original genotype files
        psam = "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common.psam",
        pvar = "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common.pvar",
        pgen = "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common.pgen",
        # IDs to keep and list of SNPs to extract
        keep = "results/phenotype/IDs_extract_{parity}.txt",
        snplist = "SNP_to_extract.txt"
    output:
        # QCed genotype files
        pgen = "results/gwas/moba_common_qc_ipi_{parity}.pgen",
        pvar = "results/gwas/moba_common_qc_ipi_{parity}.pvar",
        psam = "results/gwas/moba_common_qc_ipi_{parity}.psam"
    params:
        # base genotype file location and output path
        pfile= "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common",
        results= "results/gwas/moba_common_qc_ipi_{parity}"
    log:
        # log file capturing plink2 output and errors
        "results/gwas/moba_common_qc_ipi_{parity}.log"
    shell:
        # running plink2 with genotype and SNP QC thresholds
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
        weights = "variant_weights.txt", # SNP weights for PGS calculation
        pfile = expand("results/gwas/moba_common_qc_ipi_{{parity}}.{ext}", ext= ['pgen', 'pvar', 'psam'])
    output:
        sscore = "results/pgs/ipi_pgs_{parity}.sscore", # output score file
        log = "results/pgs/ipi_pgs_{parity}.log"
    params:
        pfile = "results/gwas/moba_common_qc_ipi_{parity}",
        out = "results/pgs/ipi_pgs_{parity}"
    shell:
        """
        mkdir -p results/pgs # ensuring output folder exists
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

# rule to restrict analysis to European ancestry (CEU) individuals and merge phenotypes with PGS
rule keep_CEU:
    input:
        pheno = "results/phenotype/filtered_pregnancies_{parity}.csv",
        pgs = "results/pgs/ipi_pgs_{parity}.sscore",
        ceu_ids = "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/pca/moba_genotypes_2024.12.03_ceu_core_ids"
    output:
        "results/phenotype/pheno_pgs_unique_{parity}.csv"
    script:
         # r script that filters for CEU samples and merges phenotypes with PGS scores
        "scripts/qc_ipi_pgs.R"

# rule to remove related individuals, merge phenotypes with PCA and batch info
rule remove_related:
        'Concat pheno files, and add PCA.'
        input:
                'results/phenotype/pheno_pgs_unique_{parity}.csv', # CEU-filtered phenotype+PGS
                '/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/kinship/moba_genotypes_2024.12.03.kin0', # Kinship matrix
                '/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/pca/moba_genotypes_2024.12.03.pcs', # principal components
		'/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/batch/moba_genotypes_2024.12.03_batches', # batch info
                'resources/batches.json' # batch harmonisation dictionary
        output:
                'results/final_phenotype/IPI_pgs_covariates_{parity}.txt'
        run:
                # loading phenotype+PGS data
                d= pd.read_csv(input[0], header= 0, sep= ',')
                # removing related individuals
                remove= selectUnrelated(input[1], d, d.SENTRIX_ID)
                d= d.loc[~d.SENTRIX_ID.isin(remove.IID.values), : ]
                # merging principal components
                pcs= pd.read_csv(input[2], header= 0, sep= '\t')
		pcs= pcs[['IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6']]
                d= pd.merge(d, pcs, right_on= 'IID', left_on= 'SENTRIX_ID')
		# merging batch information
                batches= pd.read_csv(input[3], sep= '\t', header= 0)
		batches.columns= ['IID', 'batch']
		d= pd.merge(d, batches, left_on= 'SENTRIX_ID', right_on= 'IID')
                # harmonising batch names
                with open(input[4], 'r') as fp:
                        batches_dict= json.load(fp)
                d['batch'].replace(batches_dict, inplace= True)
                d['batch']= d['batch'].str.replace(' ', '_')
                # final analysis-ready phenotype
                d.to_csv(output[0], sep= '\t', header= True, index= False)

# rule to generate figures summarising IPI and birth outcomes
rule figures_IPI:
    input:
        "/mnt/scratch/agnes/PDB1724_MFR_541_v12.csv",
        "/mnt/scratch/agnes/parental_ID_to_PREG_ID.csv",
        "/mnt/work/p1724/v12/linkage_Mother_PDB1724.csv",
        "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common.psam"
    output:
    # output various descriptive plots about IPI and outcomes
        "results/plots/phenotype/ipi_distribution_corrected.png",
        "results/plots/phenotype/ipi_density_IPI.png",
        "results/plots/phenotype/ipi_ushape.png",
        "results/plots/phenotype/ipi_delivery_type.png",
        "results/plots/phenotype/ipi_stillbirth_vs_livebirth.png",
    script:
    # r script that generates descriptive plots
        "scripts/figures_IPI_stat.R"


# innen csekkolni
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
	    params_pfile= "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common"
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

rule gxe_interaction_ipi_parameters12_gw:
    input:
        pgen = "results/gwas/moba_common_qc_ipi_multiparous_genomewide.pgen",
        pvar = "results/gwas/moba_common_qc_ipi_multiparous_genomewide.pvar",
        psam = "results/gwas/moba_common_qc_ipi_multiparous_genomewide.psam",
        keep = "results/phenotype/IDs_extract_multiparous.txt",
        pheno = "results/final_phenotype/plink_ready_IPI_pgs_covariates_multiparous.txt"
    output:
        "results/gwas/gxe_ipi_gd_gw.SVLEN_UL_DG.glm.linear",
        "results/gwas/gxe_ipi_gd_gw.log"
    params:
	pfile = "results/gwas/moba_common_qc_ipi_multiparous_genomewide"
    threads: 4
    shell:
        """
        plink2 \
          --pfile {params.pfile} \
          --keep {input.keep} \
          --pheno {input.pheno} \
          --pheno-name SVLEN_UL_DG \
          --covar {input.pheno} \
          --covar-name IPI,PGS,PC1,PC2,PC3,PC4,PC5,PC6 \
          --covar-variance-standardize \
          --glm interaction \
          --parameters 1-9 \
          --threads {threads} \
          --out results/gwas/gxe_ipi_gd_gw
        """
