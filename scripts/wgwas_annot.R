
# GWAS annotation - biomaRt + clusterProfiler
# because of I could not use biomaRT and clusterProfiler on the server, I worked on my desktop

library(biomaRt)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

# RSIDs (top 20 SNP)
file_path <- file.choose()  # the file name's on the server: top20_rsids.txt 
rsids <- readLines(file_path)

# SNP annotation - Ensembl variation // asia mirror is faster
ensembl_snp <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp", mirror = "asia")

annot_snp <- getBM(
  attributes = c("refsnp_id", "chr_name", "chrom_start", "allele", 
                 "consequence_type_tv", "ensembl_gene_stable_id"),
  filters = "snp_filter",
  values = rsids,
  mart = ensembl_snp
)

# Ensembl genes
ensembl_gene <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")

gene_names <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = unique(annot_snp$ensembl_gene_stable_id),
  mart = ensembl_gene
)

# annotation// merging
annot_final <- merge(
  annot_snp,
  gene_names,
  by.x = "ensembl_gene_stable_id",
  by.y = "ensembl_gene_id",
  all.x = TRUE
)


print(head(annot_final, 20))

write.table(
  annot_final,
  file = "~/Desktop/top20_rsids_annotation.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# annotated SNPs (with gene names only)
annotated_snps <- annot_final %>%
  filter(!is.na(external_gene_name))

print(annotated_snps)

write.table(
  annotated_snps,
  file = "~/Desktop/top20_snps_with_genes.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# unique genes
gene_list <- unique(na.omit(annotated_snps$external_gene_name))

writeLines(gene_list, con = "~/Desktop/top_genes_for_enrichment.txt")

write.csv(
  data.frame(Gene = gene_list),
  file = "~/Desktop/top_genes_for_enrichment.csv",
  row.names = FALSE,
  quote = FALSE
)

# GO enrichment analysis
# 
genes <- read.table("~/Desktop/top_genes_for_enrichment.txt", stringsAsFactors = FALSE)[,1]

# --> SYMBOL → ENTREZ ID conversion
entrez_ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# --> GO (biological process)
ego <- enrichGO(
  gene = entrez_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

head(ego)
dotplot(ego, showCategory = 10)