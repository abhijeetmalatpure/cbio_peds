library(dplyr)
library(tidyr)
library(stringr)
library(biomaRt)

dataset <- "TREEHOUSE"
cbio_code <- "cbio_PEDS"

setwd(dataset)

raw_clinical <- read.csv(file="clinical_TumorCompendium_v11_PolyA_2020-04-09.tsv", 
                     sep = "\t", header = TRUE, na.strings = c("N/A", ""))

clinical <- raw_clinical
# colnames(clinical) <- toupper(colnames(clinical))

# Set site_donor_id as substring of th_sampleid before last _
clinical$site_donor_id <- ifelse(is.na(clinical$site_donor_id) & startsWith(clinical$site_id, "TH"), 
                                 substr(clinical$th_sampleid, 0, regexpr("_[^_]*$", clinical$th_sampleid) - 1), 
                                 clinical$site_donor_id)

# find location of third underscore for each sample_id
# Turs out its 17 for 'Target' site_id and 13 for TCGA uniformly
#clinical$dashes <- str_locate_all(clinical$th_sampleid, "-")

clinical$site_donor_id <- ifelse(is.na(clinical$site_donor_id) & tolower(clinical$site_id) == "target", 
                                 substr(clinical$th_sampleid, 0, 17 - 1), 
                                 clinical$site_donor_id)

clinical$site_donor_id <- ifelse(is.na(clinical$site_donor_id) & tolower(clinical$site_id) == "tcga", 
                                 substr(clinical$th_sampleid, 0, 13 - 1), 
                                 clinical$site_donor_id)

# Update verbose values of yes to Yes
clinical$pedaya <- ifelse(clinical$pedaya == "Yes, age < 30 years", "Yes", clinical$pedaya)

clinical <- unique(clinical[, !(names(clinical) %in% "th_sampleid")])
# Final clinical data DF
header <- data.frame(site_donor_id = c("#Patient ID", "#Patient ID", "#STRING", "#1", "PATIENT_ID"),
                          age_at_dx = c("Age at Diagnosis", "Age at Diagnosis", "NUMBER", "1", "AGE"),
                          gender = c("Gender", "Gender", "STRING", "1", "GENDER"),
                          site_id = c("Site ID", "Site ID", "STRING", "1", "AGE"),
                          pedaya = c("Pediatric Adolescent and Young Adult", "Pediatric Adolescent and Young Adult", "STRING", "1", "PEDAYA"))

clinical_final <- rbind(header, clinical %>% dplyr::select(c("site_donor_id", "age_at_dx", 
                                                         "pedaya", "gender", "site_id")))

setwd(cbio_code)
clinicalFile <- "data_clinical_patient.txt"

if(file.exists(clinicalFile))
  file.remove(clinicalFile)
file.create(clinicalFile)

# Write to file
write.table(clinical_final, clinicalFile, sep="\t", col.names = FALSE, row.names = FALSE, 
            quote = FALSE, append = TRUE, na = "NA")


raw_genes <- read.csv(file="TumorCompendium_v11_PolyA_hugo_log2tpm_58581genes_2020-04-09.tsv", 
                     sep = "\t", header = TRUE, na = "NA", nrows = 1000)

samples_genes <- as.data.table(colnames(raw_genes[2:ncol(raw_genes)]))

# Get samples with blank site_donor_id
no_site_id <- raw_clinical %>% dplyr::filter(is.na(raw_clinical$site_donor_id))

# Set site_donor_id as substring of th_sampleid before last _
no_site_id$site_donor_id <- substr(no_site_id$th_sampleid, 0, regexpr("_[^_]*$", no_site_id$th_sampleid) - 1)

df_clinical <- rbind(raw_clinical %>% dplyr::filter(is.na(raw_clinical$site_donor_id)))
