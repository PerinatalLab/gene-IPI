import pandas as pd
import numpy as np
import json

rule all:
    input:
        "results/phenotype/filtered_pregnancies.csv",
        "results/gwas/moba_common_qc_ipi.pgen",
        "results/gwas/moba_common_qc_ipi.pvar",
        "results/gwas/moba_common_qc_ipi.psam",
        "results/pgs/ipi_pgs.sscore",
        "results/phenotype/pheno_pgs_unique.csv",
        "results/final_phenotype/IPI_pgs_covariates.txt"

rule cleaned_data:
    input:
        "/mnt/scratch/agnes/PDB1724_MFR_541_v12.csv",
        "/mnt/scratch/agnes/parental_ID_to_PREG_ID.csv",
        "/mnt/work/p1724/v12/linkage_Mother_PDB1724.csv",
        "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common.psam" 
    output:
        "results/phenotype/filtered_pregnancies.csv",
        "results/plots/phenotype/ipi_distribution_corrected.png",
        "results/plots/phenotype/ipi_density_IPI.png",
        "results/plots/phenotype/ipi_ushape.png",
        "results/plots/phenotype/ipi_delivery_type.png",
        "results/plots/phenotype/ipi_stillbirth_vs_livebirth.png",
        "results/phenotype/IDs_extract.txt"
    script:
        "scripts/IPI_stat.R"

rule snp_qc:
    input:
        psam = "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common.psam",
        pvar = "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common.pvar",
        pgen = "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common.pgen",
        keep = "results/gwas/keep_samples.txt",
        snplist = "SNP_to_extract.txt"
    output:
        pgen = "results/gwas/moba_common_qc_ipi.pgen",
        pvar = "results/gwas/moba_common_qc_ipi.pvar",
        psam = "results/gwas/moba_common_qc_ipi.psam"
    params:
        pfile= "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common",
        results= "results/gwas/moba_common_qc_ipi"
    log:
        "results/gwas/moba_common_qc_ipi.log"
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

rule calculate_pgs:
    input:
        weights = "variant_weights.txt",
        pfile = expand("results/gwas/moba_common_qc_ipi.{ext}", ext= ['pgen', 'pvar', 'psam'])
    output:
        sscore = "results/pgs/ipi_pgs.sscore",
        log = "results/pgs/ipi_pgs.log"
    params:
        pfile = "results/gwas/moba_common_qc_ipi",
        out = "results/pgs/ipi_pgs"
    shell:
        """
        mkdir -p results/pgs
        plink2 \
          --pfile {params.pfile} \
          --score {input.weights} 1 2 3 cols=scoresums \
          --nonfounders \
          --out {params.out} > {output.log} 2>&1
        """

def selectUnrelated(input_kin, df, x):
        kin= pd.read_csv(input_kin, header= 0, sep= '\t')
        kin= kin.loc[kin.Kinship > 0.125, :]
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

rule keep_CEU:
    input:
        pheno = "results/phenotype/filtered_pregnancies.csv",
        pgs = "results/pgs/ipi_pgs.sscore",
        ceu_ids = "/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/pca/moba_genotypes_2024.12.03_ceu_core_ids"
    output:
        "results/phenotype/pheno_pgs_unique.csv"
    script:
        "scripts/qc_ipi_pgs.R"

rule remove_related:
        'Concat pheno files, and add PCA.'
        input:
                'results/phenotype/pheno_pgs_unique.csv',
                '/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/kinship/moba_genotypes_2024.12.03.kin0',
                '/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/pca/moba_genotypes_2024.12.03.pcs',
		'/mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/batch/moba_genotypes_2024.12.03_batches',
                'resources/batches.json'
        output:
                'results/final_phenotype/IPI_pgs_covariates.txt'
        run:
                d= pd.read_csv(input[0], header= 0, sep= '\t')
                remove= selectUnrelated(input[1], d, d.IID)
                d= d.loc[~d.IID.isin(remove.IID.values), : ]
                pcs= pd.read_csv(input[2], header= 0, sep= '\t')
		pcs= pcs[['IID', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6']]
                d= pd.merge(d, pcs, right_on= 'IID', left_on= 'SENTRIX_ID')
		batches= pd.read_csv(input[3], sep= '\t', header= 0)
		batches.columns= ['IID', 'batch']
		d= pd.merge(d, batches, left_on= 'SENTRIX_ID', right_on= 'IID')
                with open(input[4], 'r') as fp:
                        batches_dict= json.load(fp)
                d['batch'].replace(batches_dict, inplace= True)
                d['batch']= d['batch'].str.replace(' ', '_')
                d.to_csv(output[0], sep= '\t', header= True, index= False)

