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

write.table(patient, "data_clinical_patient_formatted.txt", sep="\t", col.names = FALSE, row.names = FALSE,
            quote = FALSE, append = FALSE, na = "NA")
# Read samples into memory
samples <- read.csv("data_clinical_sample_ALL_PST.txt", sep = "\t", header = FALSE, na.strings = c("N/A", "", "unavailable"))
samples_unique <- unique(samples)

# Remove samples which have names that are too long
samples_final <- samples_unique[str_length(samples_unique$V2) < 50, ]

#
write.table(samples_final, "data_clinical_sample_formatted.txt", sep="\t", col.names = FALSE, row.names = FALSE,
            quote = FALSE, append = FALSE, na = "NA")

#write.table(samples_final, "data_clinical_sample_formatted.txt", sep="\t", col.names = FALSE, row.names = FALSE,
#            quote = FALSE, append = FALSE, na = "NA")

# Read Labtest into memory, get additional headers, combine both, remove duplicate
labtest_raw <- read.csv("data_timeline_lab_test_ALL_PST.txt", sep = "\t", header = FALSE, na.strings = c("N/A", "", "unavailable"))
lt_header <- read.csv("data_extra_header_timeline_lab_test_ALL_PST.txt", sep = "\t", header = FALSE, na.strings = c("N/A", "", "unavailable"))
lt_colnames <- cbind(labtest_raw[1,1:6], lt_header[1,])
labtest <- read.csv("data_timeline_lab_test_ALL_PST.txt", sep = "\t", header = FALSE,
                    col.names = toupper(lt_colnames), na.strings = c("N/A", "", "unavailable"))
labtest <- unique(labtest)

# Set sample id as P<PatientID>_D<StartDate> for labtest entries because no sample_ids are known
labtest$SPECIMEN_REFERENCE_NUMBER <- paste0("P", labtest$PATIENT_ID, "_D", labtest$START_DATE)
samples_lt <- unique(labtest[2:nrow(labtest), c('PATIENT_ID', 'SPECIMEN_REFERENCE_NUMBER')])

lt_transpose <- labtest[2:nrow(labtest), c('TEST', 'RESULT', 'PATIENT_ID', 'SPECIMEN_REFERENCE_NUMBER')]
lt_transpose$TESTCODE <- toupper(gsub(' ', '_', gsub( '[^[:alnum:][:space:]]','', lt_transpose$TEST)))

# Transpose matrix
tests <- data.frame(matrix(ncol= length(unique(lt_transpose$TESTCODE)), nrow = length(unique(lt_transpose$SPECIMEN_REFERENCE_NUMBER))+5))
tests$SPECIMEN_REFERENCE_NUMBER <- c(samples[1:5, 'V2'], unique(lt_transpose$SPECIMEN_REFERENCE_NUMBER))
colnames(tests)[1:ncol(tests)-1] <- unique(lt_transpose$TESTCODE)

# For each row in lt_transpose, move values into tests dataframe
for (row in 1:nrow(lt_transpose)) {
  specimen <- lt_transpose[row, "SPECIMEN_REFERENCE_NUMBER"]
  tst <- lt_transpose[row, "TESTCODE"]
  result <- lt_transpose[row, "RESULT"]
  tests[tests$SPECIMEN_REFERENCE_NUMBER == specimen, tst] <- result
}

# Fix first five rows
tests[1:2, ] <- colnames(tests)
tests[3, ] <- 'STRING'
tests[4, ] <- 1
tests[5, ] <- colnames(tests)
tests[1:5, 'SPECIMEN_REFERENCE_NUMBER'] <- samples[1:5, 'V2']

# Append samples from labtest to samples data table
# samples_final <- full_join(samples_, unique(labtest[2:nrow(labtest), c('SPECIMEN_REFERENCE_NUMBER', 'PATIENT_ID')]), by = c('V2' = 'SPECIMEN_REFERENCE_NUMBER'))
# samples_final$V1 <- ifelse(is.na(samples_final$V1), samples_final$PATIENT_ID, samples_final$V1)

samples_final <- left_join(samples_, tests, by = c('V2' = 'SPECIMEN_REFERENCE_NUMBER'))

samples_final$V1 <- ifelse(is.na(samples_final$V1),
                           substr(samples_final$V2, 2, str_locate(samples_final$V2, '_') -1),
                           samples_final$V1)

# Write samples, labtest to file
write.table(labtest[2:nrow(labtest),],"data_timeline_lab_test_formatted.txt", sep="\t", col.names = TRUE, row.names = FALSE,
            quote = FALSE, append = FALSE, na = "")

# Read treatment, replace non-numeric start dates with blanks, remove rows with NA patient ids, drop duplicates
treatment <- read.csv("data_timeline_treatment_ALL_PST.txt", sep = "\t", header = TRUE, na.strings = c("N/A", "", "unavailable"))
treatment$START_DATE <- ifelse(is.integer(as.integer(treatment$START_DATE)), treatment$START_DATE, "")
treatment <- treatment[!is.na(treatment$PATIENT_ID), ]
treatment <- unique(treatment)
#treatment$SPECIMEN_REFERENCE_NUMBER <- paste0("P", treatment$PATIENT_ID, "_D", treatment$START_DATE)
treatment$DOSE <- substr(treatment$DOSAGE, 0, str_locate(treatment$DOSAGE, ', ')-1)

# Write to file
write.table((treatment), "data_timeline_treatment_formatted.txt", sep="\t", col.names = TRUE, row.names = FALSE,
            quote = FALSE, append = FALSE, na = "")

