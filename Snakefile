rule all: 
    "collect the results"
    input: 
        "results/phenotype/filtered_pregnancies.csv"

rule cleaned_data: 
    "clean IPI and parity data"
    input:
        "/mnt/scratch/agnes/PDB1724_MFR_541_v12.csv", 
        "/mnt/scratch/agnes/parental_ID_to_PREG_ID.csv"
    output:
        "results/phenotype/filtered_pregnancies.csv",
        "results/plots/phenotype/ipi_distribution_corrected.png",
        "results/plots/phenotype/ipi_density_IPI.png",
        "results/plots/phenotype/ipi_parity_boxplot.png",
        "results/plots/phenotype/ipi_ushape.png",
        "results/plots/phenotype/ipi_gestational_duration.png",
        "results/plots/phenotype/ipi_delivery_type.png",
        "results/plots/phenotype/ipi_induced_vs_spont.png",
        "results/plots/phenotype/ipi_preterm_vs_term.png",
        "results/plots/phenotype/ipi_delivery_type_by_IPI_and_preterm.png",
        "results/plots/phenotype/ipi_distribution_mode_of_delivery.png",
        "results/plots/phenotype/ipi_c_section.png",        
        "results/plots/phenotype/ipi_stillbirth_vs_livebirth.png"
    script:
        "scripts/cleancode.R"
