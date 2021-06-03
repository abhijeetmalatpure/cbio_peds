library(dplyr)
library(tidyr)
library(stringr)
#library(biomaRt)
library(data.table)


# Set dir to cbio_PEDS
dataset <- "c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS"
setwd(dataset)

# Read patients into memory
patient <- read.csv("data_clinical_patient_ALL_PST.txt", sep = "\t", header = FALSE, na.strings = c("N/A", "", "unavailable"))

old_patient <- read.csv("data_clinical_patient_formatted.txt", sep = "\t", header = FALSE, na.strings = c("N/A", "", "unavailable"))

# Read samples into memory
samples <- read.csv("data_clinical_sample_ALL_PST.txt", sep = "\t", header = FALSE, na.strings = c("N/A", "", "unavailable"))
samples_unique <- unique(samples)

old_samples <- read.csv("data_clinical_sample_formatted.txt", sep = "\t", header = FALSE, na.strings = c("N/A", "", "unavailable"))

# Remove samples which have names that are too long
samples_final <- samples_unique[str_length(samples_unique$V2) < 50, ]

samples_final <- rbindlist(list(samples_final, old_samples %>% dplyr::filter(!(V2 %in% samples_final$V2))), fill = TRUE)
#
write.table(as.data.frame(samples_final), "data_clinical_sample_formatted.txt", sep="\t", col.names = FALSE, row.names = FALSE,
            quote = FALSE, append = FALSE, na = "NA")

#write.table(samples_final, "data_clinical_sample_formatted.txt", sep="\t", col.names = FALSE, row.names = FALSE,
#            quote = FALSE, append = FALSE, na = "NA")

# Read Labtest into memory, get additional headers, combine both, remove duplicate
labtest_raw <- read.csv("data_timeline_lab_test_ALL_PST.txt", sep = "\t", header = FALSE, na.strings = c("N/A", "", "unavailable"))
lt_header <- read.csv("data_extra_header_timeline_lab_test_ALL_PST.txt", sep = "\t", header = FALSE, na.strings = c("N/A", "", "unavailable"))
lt_colnames <- cbind(labtest_raw[1,1:7], lt_header[1,])
labtest <- read.csv("data_timeline_lab_test_ALL_PST.txt", sep = "\t", header = FALSE,
                    col.names = toupper(lt_colnames), na.strings = c("N/A", "", "unavailable"))
labtest <- unique(labtest)

# Data engineering

# Tobacco smoking
cat_never_smoked <- c("Never (less than 100 in lifetime)", "Never Smoker", "Never Smoked")
patient[1:5, "V6"] <- c('Tobacco Smoking', 'Tobacco Smoking', 'STRING', 1, 'TOBACCO_SMOKING')
smoking <- labtest %>% dplyr::filter(TEST == 'Tobacco: Smoking History') %>% dplyr::group_by(PATIENT_ID, TEST)
# If value is any of the never categories,
smoking$boolNever <- ifelse(smoking$RESULT %in% cat_never_smoked, 0, 1)
neverSmoked <- smoking[, c("PATIENT_ID", "boolNever")] %>% summarise(never = mean(boolNever))

patient$V6[patient$V1 %in% neverSmoked$PATIENT_ID] <- 'Never Smoked'

# Remove all Tobacco smoking tests associated with neverSmoked patients
labtest <- labtest %>% filter(!(PATIENT_ID %in% neverSmoked$PATIENT_ID & TEST == 'Tobacco: Smoking History'))

# Alcohol use
alcohol <- labtest %>% dplyr::filter(TEST == 'Alcohol Use') %>% dplyr::group_by(PATIENT_ID, TEST)



# Cause of death
death_cause <- labtest[labtest$TEST == 'Cause of Death', c("PATIENT_ID", "RESULT")]
colnames(death_cause)[2] <- "V7"
patient <- left_join(patient, death_cause, by = c('V1' = 'PATIENT_ID'))
patient[1:5, "V7"] <- c('Cause of Death', 'Cause of Death', 'STRING', 1, 'CAUSE_OF_DEATH')
labtest <- labtest %>% filter(!(PATIENT_ID %in% death_cause$PATIENT_ID & TEST == 'Cause of Death'))


# unique(labtest[labtest$TEST == 'Deceased', 'RESULT']) is Yes,
# so just make sure patient status are Deceased for those patients and then remove labtests
death <- labtest[labtest$TEST == 'Deceased', c("PATIENT_ID", "RESULT")]
death <- left_join(death, patient[6:nrow(patient), c("V1", "V2")], by = c('PATIENT_ID' = 'V1'))

if(unique(death$V2) == 'DECEASED') {
  labtest <- labtest %>% filter(!(PATIENT_ID %in% death$PATIENT_ID & TEST == 'Deceased'))
}


# Write patients table to file
write.table(patient, "data_clinical_patient_formatted.txt", sep="\t", col.names = FALSE, row.names = FALSE, quote = FALSE, append = FALSE, na = "NA")


# # This did not work.
# # Set sample id as P<PatientID>_D<StartDate> for labtest entries because no sample_ids are known
# labtest$SPECIMEN_REFERENCE_NUMBER <- paste0("P", labtest$PATIENT_ID, "_D", labtest$START_DATE)
# samples_lt <- unique(labtest[2:nrow(labtest), c('PATIENT_ID', 'SPECIMEN_REFERENCE_NUMBER')])
#
# lt_transpose <- labtest[2:nrow(labtest), c('TEST', 'RESULT', 'PATIENT_ID', 'SPECIMEN_REFERENCE_NUMBER')]
# lt_transpose$TESTCODE <- toupper(gsub(' ', '_', gsub( '[^[:alnum:][:space:]]','', lt_transpose$TEST)))
#
# # Transpose matrix
# tests <- data.frame(matrix(ncol= length(unique(lt_transpose$TESTCODE)), nrow = length(unique(lt_transpose$SPECIMEN_REFERENCE_NUMBER))+5))
# tests$SPECIMEN_REFERENCE_NUMBER <- c(samples[1:5, 'V2'], unique(lt_transpose$SPECIMEN_REFERENCE_NUMBER))
# colnames(tests)[1:ncol(tests)-1] <- unique(lt_transpose$TESTCODE)
#
# # For each row in lt_transpose, move values into tests dataframe
# for (row in 1:nrow(lt_transpose)) {
#   specimen <- lt_transpose[row, "SPECIMEN_REFERENCE_NUMBER"]
#   tst <- lt_transpose[row, "TESTCODE"]
#   result <- lt_transpose[row, "RESULT"]
#   tests[tests$SPECIMEN_REFERENCE_NUMBER == specimen, tst] <- result
# }
#
# # Fix first five rows
# tests[1:2, ] <- colnames(tests)
# tests[3, ] <- 'STRING'
# tests[4, ] <- 1
# tests[5, ] <- colnames(tests)
# tests[1:5, 'SPECIMEN_REFERENCE_NUMBER'] <- samples[1:5, 'V2']

# Append samples from labtest to samples data table
# samples_final <- full_join(samples_, unique(labtest[2:nrow(labtest), c('SPECIMEN_REFERENCE_NUMBER', 'PATIENT_ID')]), by = c('V2' = 'SPECIMEN_REFERENCE_NUMBER'))
# samples_final$V1 <- ifelse(is.na(samples_final$V1), samples_final$PATIENT_ID, samples_final$V1)

# samples_final <- left_join(samples_, tests, by = c('V2' = 'SPECIMEN_REFERENCE_NUMBER'))

# samples_final$V1 <- ifelse(is.na(samples_final$V1),
#                            substr(samples_final$V2, 2, str_locate(samples_final$V2, '_') -1),
#                            samples_final$V1)

# Write samples, labtest to file
write.table(labtest[2:nrow(labtest),],"data_timeline_lab_test_formatted.txt", sep="\t", col.names = TRUE, row.names = FALSE,
            quote = FALSE, append = FALSE, na = "")

# Read treatment, replace non-numeric start dates with blanks, remove rows with NA patient ids, drop duplicates
treatment <- read.csv("data_timeline_treatment_ALL_PST.txt", sep = "\t", header = TRUE, na.strings = c("N/A", "", "unavailable"))
#treatment$X_DATE <- ifelse(is.integer(treatment$START_DATE), treatment$START_DATE, "")
#treatment <- treatment[!is.na(treatment$PATIENT_ID), ]
treatment <- unique(treatment)
#treatment$SPECIMEN_REFERENCE_NUMBER <- paste0("P", treatment$PATIENT_ID, "_D", treatment$START_DATE)
treatment$DOSE <- substr(treatment$DOSAGE, 0, str_locate(treatment$DOSAGE, ', ')-1)

# Write to file
write.table((treatment), "data_timeline_treatment_formatted.txt", sep="\t", col.names = TRUE, row.names = FALSE,
            quote = FALSE, append = FALSE, na = "")

