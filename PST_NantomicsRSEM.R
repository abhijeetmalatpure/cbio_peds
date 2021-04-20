# Title     : Run tximport to get gene level RNA Seq data, merge with previous rnaseq file
# Objective : TODO
# Created by: abhmalat
# Created on: 4/5/2021

library(dplyr)
library(stringr)
library(biomaRt)
library(data.table)
library(tximport)
library(stringr)


# rsem <- read.csv("c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS/data/nantomics/rsem/rsem_combined.tsv", header=TRUE, sep="\t")
#
# rsem_file <- c("c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS/data/nantomics/rsem/rsem_combined.tsv")

read_tsv <- function(file) {
  return(read.csv(file, sep="\t", header = TRUE))
}

# names(rsem)
# txi.rsem <- tximport(rsem_file, 'rsem', FALSE, FALSE,
#                      geneIdCol = "gene", txIdCol = "transcript", abundanceCol = "TPM",
#                      countsCol = "expected_count", lengthCol = "effective_length", importer = read_tsv)

rsem_files <- file.path(list.files("c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS/data/nantomics/rsem", pattern = "*rsem.txt"))

get_sample <- function (x) {
  rsem_file <- unlist(strsplit(unlist(strsplit(x, "_"))[6], "\\."))[1]
  return(metadata[metadata$somatic_tumor_rna == rsem_file, "ContrastUUID"])
}

metadata <- read.csv("data/nantomics/peds_vcf_metadata.csv", sep=',', header=TRUE, na.strings = "") %>%
  dplyr::select(c("ContrastUUID", "somatic_tumor_rna"))  %>%
  dplyr::filter(!is.na(somatic_tumor_rna))

names(rsem_files)<- unlist(lapply(rsem_files, get_sample))

setwd("c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS/data/nantomics/rsem")
txi.rsem <- tximport(rsem_files, 'rsem', FALSE, FALSE,
                     geneIdCol = "gene_id", txIdCol = "transcript", abundanceCol = "TPM",
                     countsCol = "expected_count", lengthCol = "effective_length", importer = read_tsv)

# Write to file. Set col.names to NA to capture row names in it's own column
write.table(as.data.frame(txi.rsem$counts), "rsem_counts.txt", sep="\t",
             col.names = NA, quote = FALSE, append = FALSE, na = "")

setwd("c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS")
nant_rna <- read.csv("data/nantomics/rsem/rsem_counts.txt", sep="\t", header = FALSE)
nant_rna$V1[1] <- "Hugo_Symbol"

# This is the second phase - combine old rnaseq data with the nant rsem rna
ensembl <- useMart(host = 'grch37.ensembl.org',
                    biomart = 'ENSEMBL_MART_ENSEMBL',
                    dataset = 'hsapiens_gene_ensembl')

Entrez_Gene_Ids <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = 'hgnc_symbol',
                          values=nant_rna$V1, ensembl)

colnames(Entrez_Gene_Ids) <- c("V1", "Entrez_Gene_Id")
nant_rna_entrez <- left_join(nant_rna, Entrez_Gene_Ids, by="V1")
nant_rna_entrez$Entrez_Gene_Id[1] <- "Entrez_Gene_Id"

# Remove NA entrez_id rows, remove Hugo symbol column, move Entrez_id column to first position
nant_rna_entrez <- nant_rna_entrez %>%
  dplyr::filter(!is.na(Entrez_Gene_Id)) %>%
  dplyr::select(-c("V1")) %>%
  dplyr::select(Entrez_Gene_Id, everything())

write.table(as.data.frame(nant_rna_entrez), "data/nantomics/rsem/rsem_counts_entrez.txt", sep="\t",
             col.names = FALSE, row.names = FALSE, quote = FALSE, append = FALSE, na = "")

# Run  python git/datahub-study-curation-tools/generate_Zscores/NormalizeExpressionLevels_allsampleref.py
#       -i ~/rsem_counts_entrez.txt -o ~/rsem_counts_entrez_Zscore.txt -l
# before proceeding below

old_rna <- read.csv("data_rna_expression.txt", sep="\t", header = FALSE)

names(old_rna)[1] <- "Entrez_Gene_Id"

new_rna <- full_join(old_rna, nant_rna_entrez, by = "Entrez_Gene_Id")

write.table(as.data.frame(new_rna), "data_rna_expression_2.txt", sep="\t",
             col.names = FALSE, row.names = FALSE, quote = FALSE, append = FALSE, na = "")

print("RSEM incorporated into old rnaseq_expression file")