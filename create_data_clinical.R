# Something is up with Patient SJINF068_1:
# Patient SJINF068_1 being added for the first time. Apparently this patient was not in the samples file, or the samples file is not yet loaded (should be loaded before this one); 1x

library(dplyr)
library(tidyr)
library(stringr)
library(biomaRt)
library(data.table)

dataset <- "c:/Users/abhmalat/OneDrive - Indiana University/RI cBioPortal/PEDS_brrhelm/TREEHOUSE"
cbio_code <- "c:/Users/abhmalat/OneDrive - Indiana University/cbio_PEDSTreehouse/study"

setwd(dataset)

raw_clinical <- read.csv(file = "clinical_TumorCompendium_v11_PolyA_2020-04-09.tsv",
                         sep = "\t", header = TRUE, na.strings = c("N/A", "", "unavailable"))

clinical <- raw_clinical
# colnames(clinical) <- toupper(colnames(clinical))

# Transform patient/sample names - replace offending chars and create mapping
clinical$sample_id <- clinical$th_sampleid
clinical$patient_id <- clinical$site_donor_id

# Set site_donor_id as substring of th_sampleid before last _
clinical$patient_id <- ifelse(is.na(clinical$patient_id) & startsWith(clinical$site_id, "TH"),
                                 substr(clinical$th_sampleid, 0, regexpr("_[^_]*$", clinical$th_sampleid) - 1),
                                 clinical$patient_id)

# find location of third underscore for each sample_id
# Turs out its 17 for 'Target' site_id and 13 for TCGA uniformly
#clinical$dashes <- str_locate_all(clinical$th_sampleid, "-")

clinical$patient_id <- ifelse(is.na(clinical$patient_id) & tolower(clinical$site_id) == "target",
                                 substr(clinical$th_sampleid, 0, 17 - 1),
                                 clinical$patient_id)

clinical$patient_id <- ifelse(is.na(clinical$site_donor_id) & tolower(clinical$site_id) == "tcga",
                                 substr(clinical$th_sampleid, 0, 13 - 1),
                                 clinical$patient_id)

clinical$patient_id <- gsub(" |/", "_", clinical$patient_id)
clinical$patient_id <- gsub("^icgc__", "", clinical$patient_id)
clinical$patient_id <- gsub("[()]", "", clinical$patient_id)
clinical$patient_id <- gsub("_withJunctionsOnGenome_dupsFlagged", "", clinical$patient_id)
clinical$patient_id <- gsub("_RNA-seq_PAIRED_ICGC", "", clinical$patient_id)
clinical$patient_id <- gsub("_RNA-Seq_paired-end_ICGC", "", clinical$patient_id)
clinical$patient_id <- gsub("RNA-Seq_paired-end_MB-", "", clinical$patient_id)
clinical$patient_id <- gsub("_RNA_tumor_ICGC", "", clinical$patient_id)

# Update verbose values of yes to Yes
clinical$pedaya <- ifelse(clinical$pedaya == "Yes, age < 30 years", "Yes", clinical$pedaya)

# clinical <- unique(clinical[, !(names(clinical) %in% "th_sampleid")])
# Final clinical data DF
header <- data.frame(patient_id = c("#Patient ID", "#Patient ID", "#STRING", "#1", "PATIENT_ID"),
                     site_donor_id = c("Site Donor ID", "Site Donor ID", "STRING", "1", "SITE_DONOR_ID"),
                     age_at_dx = c("Age at Diagnosis", "Age at Diagnosis", "NUMBER", "1", "AGE"),
                     gender = c("Gender", "Gender", "STRING", "1", "GENDER"),
                     site_id = c("Site ID", "Site ID", "STRING", "1", "SITE_ID"),
                     pedaya = c("Pediatric Adolescent and Young Adult", "Pediatric Adolescent and Young Adult", "STRING", "1", "PEDAYA"))


clinical_final <- rbind(header, unique(clinical %>% dplyr::select(c("patient_id", "site_donor_id", "age_at_dx",
                                                                    "gender", "site_id", "pedaya"))))

setwd(cbio_code)
clinicalFile <- "data_clinical_patient.txt"

if(file.exists(clinicalFile))
  file.remove(clinicalFile)
file.create(clinicalFile)

# Write to file
write.table(clinical_final, clinicalFile, sep="\t", col.names = FALSE, row.names = FALSE,
            quote = FALSE, append = TRUE, na = "NA")

clinicalMetaFile <- "meta_clinical_patient.txt"

f <- file(clinicalMetaFile)
writeLines(c(
  "cancer_study_identifier: peds_treehouse",
  "genetic_alteration_type: CLINICAL",
  "datatype: PATIENT_ATTRIBUTES",
  paste("data_filename: ", clinicalFile)
), f
)
close(f)
print("Clinical completed")


samples <- clinical[, c("patient_id", "sample_id", "th_sampleid", "site_donor_id", "disease", "dataset_accession", "site_sampleid")]

samples$dataset_accession <- ifelse(samples$dataset_accession == 'unavailable', NA, samples$dataset_accession)

# Get Waqas' mapping doc
setwd(dataset)
cancer_types <- read.csv('cancer_type.csv', sep=',', header = TRUE, na.strings = "")

cancer_types <- cancer_types[,c("CANCER_SUBTYPES", "CANCER_TYPE", "CANCER_TYPE_Tumor.Site")]
names(cancer_types) <- gsub("CANCER_TYPE_Tumor.Site","METASTATIC_SITE", names(cancer_types))

samples <- left_join(samples, cancer_types, by= c("disease" = "CANCER_SUBTYPES"))

# Capitalize first letter
samples[,c("disease", "CANCER_TYPE", "METASTATIC_SITE")] <- sapply(samples[,c("disease", "CANCER_TYPE", "METASTATIC_SITE")], stringr::str_to_title)

header <- data.frame(patient_id = c("#Patient ID", "#Patient ID", "#STRING", "#1", "PATIENT_ID"),
                     sample_id = c("Sample ID", "Sample ID", "STRING", "1", "SAMPLE_ID"),
                     th_sampleid = c("Treehouse Sample ID", "Treehouse Sample ID", "STRING", "1", "TH_SAMPLE_ID"),
                     CANCER_TYPE = c("Cancer Type", "Cancer Type", "STRING", "1", "CANCER_TYPE"),
                     disease = c("Cancer Sub-type", "Cancer Sub-type", "STRING", "1", "CANCER_TYPE_DETAILED"),
                     METASTATIC_SITE = c("Site of Tumor", "Site of Tumor", "STRING", "1", "METASTATIC_SITE"),
                     site_sampleid = c("Site Sample ID", "Site Sample ID", "STRING", "1", "SITE_SAMPLE_ID"),
                     dataset_accession = c("Dataset Accession", "Dataset Accession", "STRING", "1", "DATASET_ACCESSION"))

sample_final <- rbind(header, samples %>% dplyr::select(c("patient_id", "sample_id", "th_sampleid", "CANCER_TYPE", "disease", "METASTATIC_SITE", "site_sampleid", "dataset_accession")))

setwd(cbio_code)
samplesFile <- "data_clinical_sample.txt"

if(file.exists(samplesFile))
  file.remove(samplesFile)
file.create(samplesFile)

# Write to file
write.table(sample_final, samplesFile, sep="\t", col.names = FALSE, row.names = FALSE,
            quote = FALSE, append = TRUE, na = "NA")

samplesMetafile <- "meta_clinical_sample.txt"

f <- file(samplesMetafile)
writeLines(c(
  "cancer_study_identifier: peds_treehouse",
  "genetic_alteration_type: CLINICAL",
  "datatype: SAMPLE_ATTRIBUTES",
  paste("data_filename: ", samplesFile)
), f
)
close(f)

print("Samples data completed")

# Write to file cancer types and subtypes
#write.table(unique(samples[, c("disease", "CANCER_TYPE")]), "cancer_types.tsv", sep="\t", col.names = TRUE, row.names = FALSE,
#            quote = FALSE, append = TRUE, na = "")

#
#raw_genes <- read.csv(file="TumorCompendium_v11_PolyA_hugo_log2tpm_58581genes_2020-04-09.tsv",
#                     sep = "\t", header = TRUE, na = "NA", nrows = 1000)
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

setwd(dataset)
cancer_types <- read.csv('cancer_type.csv', sep=',', header = TRUE, na.strings = "")
cancer_types <- subset(cancer_types, select = c("CANCER_SUBTYPES", "CANCER_TYPE"))