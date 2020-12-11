library(dplyr)
library(tidyr)
library(stringr)
library(biomaRt)
library(data.table)

# Send file names to combine together. They need to have identical colnames to be combined
combine_data <- function (filenames) {
  merge_df <- data.table()

  for (i in filenames) {
    # Set df_name to the vendor name
    df <- read.csv(i, sep = "\t", header = TRUE, na.strings = c("N/A", "", "unavailable"))

    if (nrow(merge_df) == 0) {
      merge_df <- df
    }
    else {
      # Check for headers match
      if(identical(colnames(merge_df), colnames(df))) {
        merge_df <- rbind(merge_df, df[5:nrow(df), ])
      }
    }
  }
  return(merge_df)
}

# Write to a file and create metaFile
write_to_file <- function(data, dataFile, metaFile, datatype) {
  write.table(data, dataFile, sep="\t", col.names = FALSE, row.names = FALSE,
            quote = FALSE, append = FALSE, na = "NA")

  if(file.exists(metaFile))
    file.remove(metaFile)
  file.create(metaFile)

  f <- file(metaFile)
  writeLines(c(
    "cancer_study_identifier: pst_peds",
    "genetic_alteration_type: CLINICAL",
    paste("datatype: ", datatype),
    paste("data_filename: ", dataFile)
  ), f
  )
  close(f)
}

# Set dir to cbio_PEDS
dataset <- "c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS"
setwd(dataset)

# Get list of file names for clinical, sample and timeline_lab_test
clinical_filenames <- list.files(pattern = "data_clinical_patient_.*")
sample_filenames <- list.files(pattern = "data_clinical_sample_.*")
labtest_filenames <- list.files(pattern = "data_timeline_lab_test_.*")

#Combine clinical data from all vendors
clinical <- combine_data(clinical_filenames)

write_to_file(clinical,"data_clinical_patient.txt",
              "meta_clinical_patient.txt", "CLINICAL_ATTRIBUTES")
print("Clinical completed")


# Remove samples which dont have patient_id or sample_id or both
sample <- combine_data(sample_filenames)
sample <- sample[complete.cases(sample),]

dataFile <- "data_clinical_sample.txt"
metaFile <- "meta_clinical_sample.txt"


write_to_file(sample,"data_clinical_sample.txt",
              "meta_clinical_sample.txt", "SAMPLE_ATTRIBUTES")
print("Samples completed")

#
##Combine timeline labtest data from all vendors
#labtest <- combine_data(labtest_filenames)
#
#write_to_file(labtest,"data_timeline_lab_test.txt",
#              "meta_timeline_lab_test.txt", "TIMELINE")
#print("Timeline Data completed")

