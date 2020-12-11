library(dplyr)
library(tidyr)
library(stringr)
library(biomaRt)
library(data.table)
source('pst_merge_clinical_files.R')

# Send file names to combine together. They need to have identical colnames to be combined
combine_timeline_data <- function (filenames) {
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
        merge_df <- rbind(merge_df, df)
      }
    }
  }
  return(merge_df)
}

# Set dir to cbio_PEDS
dataset <- "c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS"
setwd(dataset)

# Get list of file names for timeline_lab_test
labtest_filenames <- list.files(pattern = "data_timeline_lab_test_.*")

#Combine timeline labtest data from all vendors
labtest <- combine_timeline_data(labtest_filenames)

write_to_file(labtest,"data_timeline_lab_test.txt",
              "meta_timeline_lab_test.txt", "TIMELINE")
print("Timeline Data completed")

