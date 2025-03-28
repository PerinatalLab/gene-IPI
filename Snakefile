rule all:
    input:
        "results/phenotype/filtered_pregnancies.csv",
        "results/gwas/moba_common_qc_ipi.pgen",
        "results/gwas/moba_common_qc_ipi.pvar",
        "results/gwas/moba_common_qc_ipi.psam"
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
    log:
        "results/gwas/moba_common_qc_ipi.log"
    shell:
        """
        plink2 \
          --pfile /mnt/archive/moba/geno/HDGB-MoBaGenetics/2024.12.03/moba_genotypes_2024.12.03_common \
          --keep {input.keep} \
          --extract {input.snplist} \
          --maf 0.01 \
          --geno 0.01 \
          --mind 0.01 \
          --hwe 1e-6 midp \
          --make-pgen \
          --out results/gwas/moba_common_qc_ipi > {log} 2>&1
        """
