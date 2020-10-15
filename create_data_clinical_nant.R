library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(config)

# Set WD to this file's dir to find appropriate
#Sys.setenv(R_CONFIG_ACTIVE = "nantomics")
#
#config <- config::get()

dataset <- "c:/Users/abhmalat/OneDrive - Indiana University/RI cBioPortal/PEDS_brrhelm"
clinical <- "c:/Users/abhmalat/OneDrive - Indiana University/RI cBioPortal/PEDS_brrhelm/mapping"
study <- "c:/Users/abhmalat/OneDrive - Indiana University/cbio_PEDSNantomics/study"

setwd(clinical)

# Header file was separate
clinical_header <- read.csv(file = "dt_matched_patient_info_pst_nantomics_2020_1014a_pst_files_per_patient_listing_field_names.tsv", sep=",", header=FALSE)

raw_clinical <- read.csv(file = "dt_matched_patient_info_pst_nantomics_2020_1014a_pst_files_per_patient_listing.tsv",
                         sep = "|", header = FALSE, na.strings = c("N/A", ""))
names(raw_clinical) <- clinical_header$V1

# Final clinical data DF
header <- data.frame(nantomics_patient_uuid = c("#Nantomics Patient ID", "#Nantomics Patient ID", "#STRING", "#1", "PATIENT_ID"),
                     phc_patient_uuid = c("PHC Patient ID", "PHC Patient ID", "STRING", "1", "PHC_PATIENT_ID"),
                     phc_subject_id = c("Patient Subject ID", "Patient Subject ID in REDCap", "STRING", "1", "PHC_SUBJECT_ID"))

clinical_final <- rbind(header, unique(raw_clinical %>% dplyr::select(c("nantomics_patient_uuid", "phc_patient_uuid", "phc_subject_id"))))

setwd(study)
clinicalFile <- "data_clinical_patient_new.txt"

if(file.exists(clinicalFile))
  file.remove(clinicalFile)
file.create(clinicalFile)

# Write to file
write.table(clinical_final, clinicalFile, sep="\t", col.names = FALSE, row.names = FALSE,
            quote = FALSE, append = TRUE, na = "NA")

clinicalMetaFile <- "meta_clinical_patient.txt"

f <- file(clinicalMetaFile)
writeLines(c(
  "cancer_study_identifier: peds_nantomics",
  "genetic_alteration_type: CLINICAL",
  "datatype: PATIENT_ATTRIBUTES",
  paste("data_filename: ", clinicalFile)
), f
)
close(f)
print("Clinical completed")

# Get the first 8 chars from contrast id. These are sample IDs
raw_clinical$sample_id <- substr(raw_clinical$nantomics_contrast_id, 1, 8)

setwd(dataset)

sample_names <- data.table(names(read.csv(file="nant_expr_tpmlog2n.csv", sep=",", header = TRUE, stringsAsFactors=FALSE, nrows = 0)))

names(sample_names) <- "sample_id"
# Remove the first X
sample_names <- sample_names[2:nrow(sample_names)]

# Remove the X prefix that R added to columns starting with a numeral
sample_names$sample_id <- ifelse(startsWith(sample_names$sample_id, "X"), substr(sample_names$sample_id, 2, nchar(sample_names$sample_id)), sample_names$sample_id)

sample_names <- left_join(sample_names,
                          unique(raw_clinical[, c("sample_id", "nantomics_patient_uuid", "nantomics_barcode", "nantomics_contrast_id")]),
                          by = "sample_id")

header <- data.frame(nantomics_patient_uuid = c("#Nantomics Patient ID", "#Nantomics Patient ID", "#STRING", "#1", "PATIENT_ID"),
                     sample_id = c("Sample ID", "Sample ID", "STRING", "1", "SAMPLE_ID"),
                     nantomics_barcode = c("Nantomics Barcode", "Nantomics Barcode", "STRING", "1", "BARCODE"),
                     nantomics_contrast_id = c("Nantomics Contrast ID", "Nantomics Contrast ID", "STRING", "1", "CONTRAST_ID"))

sample_final <- rbind(header,sample_names[complete.cases(sample_names), ])

setwd(study)
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
  "cancer_study_identifier: peds_nantomics",
  "genetic_alteration_type: CLINICAL",
  "datatype: SAMPLE_ATTRIBUTES",
  paste("data_filename: ", samplesFile)
), f
)
close(f)

print("Samples data completed")
