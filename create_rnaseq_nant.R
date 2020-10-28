# Run after create_data_clinical_nant.R because it uses to filter samples
library(dplyr)
library(tidyr)
library(stringr)
library(biomaRt)

# Set WD to this file's dir to find appropriate
#Sys.setenv(R_CONFIG_ACTIVE = "nantomics")

# config <- config::get()

dataset <- "c:/Users/abhmalat/OneDrive - Indiana University/RI cBioPortal/PEDS_brrhelm"
clinical <- "c:/Users/abhmalat/OneDrive - Indiana University/RI cBioPortal/PEDS_brrhelm/matching_results_2020_1005a"
study <- "c:/Users/abhmalat/OneDrive - Indiana University/cbio_PEDSNantomics/study"

setwd(dataset)

raw_genes <- read.csv(file="nant_expr_tpmlog2n.csv", sep=",", header = TRUE, na = "NA")
colnames(raw_genes)[1] <- "Hugo_Symbol"

ensembl <- useEnsembl(host = 'grch37.ensembl.org',
                   biomart = 'ENSEMBL_MART_ENSEMBL',
                   dataset = 'hsapiens_gene_ensembl')

hugo <- raw_genes$Hugo_Symbol
entrez <- getBM(c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values = hugo, ensembl)

entrez[, "hgnc_symbol"] <- as.factor(entrez[, "hgnc_symbol"])

foo <- unique(arrange(entrez, hgnc_symbol, entrezgene_id))

names(foo) <- c("Hugo_Symbol", "Entrez_Gene_Id")
genes_with_entrez <- left_join(raw_genes, foo, by = "Hugo_Symbol")

# Since we're not using headers, update row 1 with names for Hugo and Entrez columns
#genes_with_entrez[1, "Hugo_Symbol"] <- "Hugo_Symbol"
#genes_with_entrez[1, "Entrez_Gene_Id"] <- "Entrez_Gene_Id"

genes_final <- genes_with_entrez %>% dplyr::select(Hugo_Symbol, Entrez_Gene_Id, everything())

genes_header <- as.data.table(names(genes_final))
genes_header$V1 <- ifelse(startsWith(genes_header$V1, "X"), substr(genes_header$V1, 2, nchar(genes_header$V1)), genes_header$V1)

sample_names <- data.table(colnames(genes_final))
sample_names$id <- ifelse(startsWith(sample_names$V1, "X"), substr(sample_names$V1, 2, nchar(sample_names$V1)), sample_names$V1)
# Find samples present in sample_final
sample_names$valid <- ifelse(sample_names$id %in% sample_final$sample_id | sample_names$id %in% c("Hugo_Symbol", "Entrez_Gene_Id"), TRUE, FALSE)

# Filter out sample columns that are not present in the data_patient_samples
genes <- rbind(transpose(sample_names[, c(3,2)]), genes_final, use.names = FALSE)
genes <- subset(genes, select = which(genes[1, ] == TRUE), use.names = FALSE)
genes <- genes[2:nrow(genes), ]

setwd(study)
rnaDataFile <- "data_expression_rnaseq.txt"

if(file.exists(rnaDataFile))
  file.remove(rnaDataFile)
file.create(rnaDataFile)

# Write to file
write.table(as.data.frame(genes), rnaDataFile, sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE, na = "NA")

print("RNA Expression data file completed")

rnaMetaFile <- "meta_expression_rnaseq.txt"

f <- file(rnaMetaFile)
writeLines(c(
  "cancer_study_identifier: peds_nantomics",
  "genetic_alteration_type: MRNA_EXPRESSION",
  "datatype: CONTINUOUS",
  "stable_id: rna_seq_mrna",
  "show_profile_in_analysis_tab: false",
  "profile_name: mRNA expression",
  "profile_description: RNA Gene Expression Log2TPM [Continuous]",
  paste("data_filename: ", rnaDataFile)
), f
)
close(f)
print("RNA Expression metafile completed")

setwd('case_lists')

sampleIds <- genes[1, 3:ncol(genes)]
rnaCaseListFile <- "cases_rna_seq_mrna.txt"

if (file.exists(rnaCaseListFile)) {
  file.remove(rnaCaseListFile)
}
file.create(rnaCaseListFile)

f <- file(rnaCaseListFile)
writeLines(c(
  "cancer_study_identifier: peds_nantomics",
  "stable_id: peds_nantomics_rna_seq_mrna",
  "case_list_name: RNA Seq Log2TPM",
  "case_list_description: RNA Seq Log2TPM [Continuous]",
  paste("case_list_ids: ", paste(sampleIds, collapse = '\t'))
), f
)
close(f)
print("Case lists completed")

#setwd(home)