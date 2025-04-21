Gene-IPI Project:

This repository contains the Snakefile workflow and R scripts developed for the master's thesis project, titled: CAN LONG PREGNANCY INTERVALS RESET THE EFFECT OF MATERNAL GENES?

Project summary: the project explores whether long interpregnancy intervals (IPI) can reset maternal genetic effects on gestational duration and birth outcomes. Using genotype and phenotype data, this analysis integrates SNP-based GWAS, polygenic scores (PGS), and phenotype filtering to test the resetting hypothesis.

Repository structure:

- `Snakefile`: workflow to automate genotype QC, PGS calculation, and phenotype preparation.
- `scripts/`: R scripts used for phenotype cleaning, SNP QC, and visualisation.
- `.gitignore`: specifies files and folders excluded from Git tracking (e.g., large data outputs).
- `.ipynb_checkpoints/`: auto-generated notebook checkpoints (ignored).

Workflow overview:
1. Filtering pregnancy phenotype datasets based on parity groups (nulliparous vs multiparous).
2. Performing SNP QC and extract selected variants for analysis.
3. Calculating polygenic scores (PGS) stratified by parity.
4. Merging phenotype and PGS data, filter for unrelated individuals, and prepare final analysis datasets.
5. Generating summary figures for IPI distributions and related outcomes.

Data files (raw phenotypes, genotypes, results) are not included in this repository.


Ágnes Judit Juhász
Master's Thesis Project — Molecular Biotechnology / Major in Biosciences  
University of Skövde & Sahlgrenska Academy, University of Gothenburg

Supervisors: PhD Pol Solé-Navais and Karin Ytterberg
