library(dplyr)
library(tidyr)
library(stringr)
library(biomaRt)

dataset <- "c:/Users/abhmalat/OneDrive - Indiana University/RI cBioPortal/PEDS_brrhelm/TREEHOUSE"
cbio_study <- "c:/Users/abhmalat/OneDrive - Indiana University/cbio_PEDS/study"

setwd(dataset)

raw_genes <- read.csv(file="TumorCompendium_v11_PolyA_hugo_log2tpm_58581genes_2020-04-09.tsv",
                     sep = "\t", header = FALSE, na = "NA")

colnames(raw_genes)[1] <- "Hugo_Symbol"

ensembl <- useMart(host = 'grch37.ensembl.org',
                   biomart = 'ENSEMBL_MART_ENSEMBL',
                   dataset = 'hsapiens_gene_ensembl')

entrez <- getBM(c("hgnc_symbol", "entrezgene_id"), filters = "hgnc_symbol", values = raw_genes$Hugo_Symbol[2:nrow(raw_genes)], ensembl)

entrez[, "hgnc_symbol"] <- as.factor(entrez[, "hgnc_symbol"])

foo <- unique(arrange(entrez, hgnc_symbol, entrezgene_id))

names(foo) <- c("Hugo_Symbol", "Entrez_Gene_Id")
genes_with_entrez <- left_join(raw_genes, foo, by = "Hugo_Symbol")

# Since we're not using headers, update row 1 with names for Hugo and Entrez columns
genes_with_entrez[1, "Hugo_Symbol"] <- "Hugo_Symbol"
genes_with_entrez[1, "Entrez_Gene_Id"] <- "Entrez_Gene_Id"

genes_final <- genes_with_entrez %>% dplyr::select(Hugo_Symbol, Entrez_Gene_Id, everything())

setwd(cbio_study)
rnaDataFile <- "data_expression_rnaseq.txt"

if(file.exists(rnaDataFile))
  file.remove(rnaDataFile)
file.create(rnaDataFile)

# Write to file
write.table(genes_final, rnaDataFile, sep="\t", col.names = FALSE, row.names = FALSE,
            quote = FALSE, append = TRUE, na = "NA")

print("RNA Expression data file completed")

rnaMetaFile <- "meta_expression_rnaseq.txt"

f <- file(rnaMetaFile)
writeLines(c(
  "cancer_study_identifier: peds_treehouse",
  "genetic_alteration_type: MRNA_EXPRESSION",
  "datatype: CONTINUOUS",
  "stable_id: rna_seq_mrna",
  "show_profile_in_analysis_tab: false",
  "profile_name: mRNA expression",
  "profile_description: RNA Gene Expression Log2TPM",
  paste("data_filename: ", rnaDataFile)
), f
)
close(f)
print("RNA Expression metafile completed")


#raw_clinical <- read.csv(file = "clinical_TumorCompendium_v11_PolyA_2020-04-09.tsv",
#                         sep = "\t", header = TRUE, na.strings = c("N/A", "", "unavailable"))
#
#clinical <- raw_clinical
## colnames(clinical) <- toupper(colnames(clinical))
#
## Set site_donor_id as substring of th_sampleid before last _
#clinical$site_donor_id <- ifelse(is.na(clinical$site_donor_id) & startsWith(clinical$site_id, "TH"),
#                                 substr(clinical$th_sampleid, 0, regexpr("_[^_]*$", clinical$th_sampleid) - 1),
#                                 clinical$site_donor_id)
#
## find location of third underscore for each sample_id
## Turs out its 17 for 'Target' site_id and 13 for TCGA uniformly
##clinical$dashes <- str_locate_all(clinical$th_sampleid, "-")
#
#clinical$site_donor_id <- ifelse(is.na(clinical$site_donor_id) & tolower(clinical$site_id) == "target",
#                                 substr(clinical$th_sampleid, 0, 17 - 1),
#                                 clinical$site_donor_id)
#
#clinical$site_donor_id <- ifelse(is.na(clinical$site_donor_id) & tolower(clinical$site_id) == "tcga",
#                                 substr(clinical$th_sampleid, 0, 13 - 1),
#                                 clinical$site_donor_id)
#
## Update verbose values of yes to Yes
#clinical$pedaya <- ifelse(clinical$pedaya == "Yes, age < 30 years", "Yes", clinical$pedaya)
#
## clinical <- unique(clinical[, !(names(clinical) %in% "th_sampleid")])
## Final clinical data DF
#header <- data.frame(site_donor_id = c("#Patient ID", "#Patient ID", "#STRING", "#1", "PATIENT_ID"),
#                     age_at_dx = c("Age at Diagnosis", "Age at Diagnosis", "NUMBER", "1", "AGE"),
#                     gender = c("Gender", "Gender", "STRING", "1", "GENDER"),
#                     site_id = c("Site ID", "Site ID", "STRING", "1", "SITE_ID"),
#                     pedaya = c("Pediatric Adolescent and Young Adult", "Pediatric Adolescent and Young Adult", "STRING", "1", "PEDAYA"))
#
#clinical_final <- rbind(header, unique(clinical %>% dplyr::select(c("site_donor_id", "age_at_dx",
#                                                                    "pedaya", "gender", "site_id"))))
#
#setwd(cbio_code)
#clinicalFile <- "data_clinical_patient.txt"
#
#if(file.exists(clinicalFile))
#  file.remove(clinicalFile)
#file.create(clinicalFile)
#
## Write to file
#write.table(clinical_final, clinicalFile, sep="\t", col.names = FALSE, row.names = FALSE,
#            quote = FALSE, append = TRUE, na = "NA")
#
#clinicalMetaFile <- "meta_clinical_patient.txt"
#
#f <- file(clinicalMetaFile)
#writeLines(c(
#  "cancer_study_identifier: peds_treehouse",
#  "genetic_alteration_type: CLINICAL",
#  "datatype: PATIENT_ATTRIBUTES",
#  paste("data_filename: ", clinicalFile)
#), f
#)
#close(f)
#print("Clinical completed")
# Write to file cancer types and subtypes
#write.table(unique(samples[, c("disease", "CANCER_TYPE")]), "cancer_types.tsv", sep="\t", col.names = TRUE, row.names = FALSE,
#            quote = FALSE, append = TRUE, na = "")

#
#
#samples_genes <- as.data.table(colnames(raw_genes[2:ncol(raw_genes)]))
#
## Get samples with blank site_donor_id
#no_site_id <- raw_clinical %>% dplyr::filter(is.na(raw_clinical$site_donor_id))
#
## Set site_donor_id as substring of th_sampleid before last _
#no_site_id$site_donor_id <- substr(no_site_id$th_sampleid, 0, regexpr("_[^_]*$", no_site_id$th_sampleid) - 1)
#
#df_clinical <- rbind(raw_clinical %>% dplyr::filter(is.na(raw_clinical$site_donor_id)))
