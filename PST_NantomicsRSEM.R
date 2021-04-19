# Title     : Run tximport to get gene level RNA Seq data
# Objective : TODO
# Created by: abhmalat
# Created on: 4/5/2021

library(dplyr)
library(stringr)
library(biomaRt)
library(data.table)
library(tximport)
library(stringr)



rsem <- read.csv("c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS/data/nantomics/rsem/rsem_combined.tsv", header=TRUE, sep="\t")

rsem_file <- c("c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS/data/nantomics/rsem/rsem_combined.tsv")

read_tsv <- function(file) {
  return(read.csv(file, sep="\t", header = TRUE))
}

names(rsem)
txi.rsem <- tximport(rsem_file, 'rsem', FALSE, FALSE,
                     geneIdCol = "gene", txIdCol = "transcript", abundanceCol = "TPM",
                     countsCol = "expected_count", lengthCol = "effective_length", importer = read_tsv)


rsem_files <- file.path(list.files("c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS/data/nantomics/rsem", pattern = "*rsem.txt"))

get_sample <- function (x) {
   return(unlist(strsplit(unlist(strsplit(x, "_"))[6], "\\."))[1])
}

names(rsem_files)<- unlist(lapply(rsem_files, get_sample))

txi.rsem <- tximport(rsem_files, 'rsem', FALSE, FALSE,
                     geneIdCol = "gene_id", txIdCol = "transcript", abundanceCol = "TPM",
                     countsCol = "expected_count", lengthCol = "effective_length", importer = read_tsv)

write.table(as.data.frame(txi.rsem$counts), "rsem_counts.txt", sep="\t",
             col.names = TRUE, row.names = TRUE,
             quote = FALSE, append = FALSE, na = "")

setwd("c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS")
nant_rna <- read.csv("data/nantomics/rsem/rsem_counts.txt", sep="\t", header = FALSE)

ensembl <- useMart(host = 'grch37.ensembl.org',
                    biomart = 'ENSEMBL_MART_ENSEMBL',
                    dataset = 'hsapiens_gene_ensembl')

Entrez_Gene_Ids <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = 'hgnc_symbol',
                          values=nant_rna$V1, ensembl)

colnames(Entrez_Gene_Ids) <- c("V1", "Entrez_Gene_Id")
nant_rna_entrez <- left_join(nant_rna, Entrez_Gene_Ids, by="V1")
nant_rna_entrez$Entrez_Gene_Id[1] <- "Entrez_Gene_Id"

old_rna <- read.csv("data_rna_expression.txt", sep="\t", header = FALSE)
nant_rna_entrez <- nant_rna_entrez[!is.na(nant_rna_entrez$Entrez_Gene_Id),]
nant_rna_entrez <- nant_rna_entrez %>% dplyr::select(-c("V1")) %>% dplyr::select(Entrez_Gene_Id, everything())

names(old_rna)[1] <- "Entrez_Gene_Id"

new_rna <- rbindlist(list(old_rna, nant_rna_entrez), fill = TRUE)

new_rna <- full_join(old_rna, nant_rna_entrez, by = "Entrez_Gene_Id")