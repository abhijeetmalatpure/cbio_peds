library(dplyr)
library(stringr)
library(data.table)
library(maftools)

# New version that sequentially reads and merges as it proceeds
combine_maf_files <- function (filenames) {
  for (i in seq_along(filenames)) {
    print(paste(i,"-", "Reading", filenames[i]))
    maf <- read_maf_file(filenames[i])#, vc_nonSyn = vcNames)
    if (!is.na(maf)) {
      if (i == 1) {
        final_maf <- maf
        print("Reading next file to merge into first..")
      }
      else {
        print(paste(filenames[i], "read correctly. Merging into the previous: "))
        final_maf <- merge_maf_files(final_maf, maf)

        # Write results to file after every 10 MAF file merges
        if(i%%10 == 0) {
          print("Writing intermediate results to file")
          write.table(as.data.frame(final_maf@data), "pst_nantomics_germline.maf", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE, append = FALSE, na = "")
        }
      }
    }
    else {
      failed_maf <- append(failed_maf, filenames[i])
    }
  }
  return(list("merged"= final_maf, "failed" = failed_maf))
}


# Read MAF file using read.maf. Return NaN if it errors out
read_maf_file <- function(x) {
  maf <- data.table::fread(x)
  vcNames <- unique(maf$Variant_Classification)
  tryCatch(read.maf(x, vc_nonSyn = vcNames),
           error = function(x) {print(paste(x, ": Error reading file. Skipping")); NaN})
}


# Merge two MAF files into one. Return first if merging error occurs
merge_maf_files <- function(maf_merged, maf_b) {
  tryCatch(merge_mafs(c(maf_merged, maf_b)),
           error = function(x) {print(paste(": Error merging this file:", x, "Skipping")); maf_merged})
}

mafdir <- "/N/project/phi_asha_archive/peds_pst/nantomics/vcf2mafConversion/maf/germline"
setwd(mafdir)

germlinemafs <- list.files(mafdir, pattern = "*\\.germline.maf$", recursive = FALSE)
germlinemafs
# Call combine_maf_files
processed <- combine_maf_files(germlinemafs)

# Save merged MAF file and failed list of MAFs
merged <- processed$merged
failed <- processed$failed

merged@data$Mutation_Status <- as.character(merged@data$Mutation_Status)
merged@data$Mutation_Status <- "Germline" # Or "Germline" for germline

merged@data$Tumor_Sample_Barcode <- merged@data$Matched_Norm_Sample_Barcode

# Write merged MAF file
write.table(as.data.frame(merged@data), "pst_nantomics_germline.maf", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE, append = FALSE, na = "")

# Save failed maf file list
write.table(as.data.frame(failed), "pst_nantomics_germline_failed.txt", col.names = FALSE, row.names = FALSE, sep = "\n", quote = FALSE, append = FALSE, na = "")

mutationCL <- ("cases_sequenced.txt")

f <- file(mutationCL)
writeLines(c(
  "cancer_study_identifier: PST_PEDS_2020",
  "stable_id: PST_PEDS_2020_sequenced",
  "case_list_name: Sequenced Tumors",
  "case_list_description: All Sequenced Tumors",
  paste("case_list_ids: ", paste(unique(merged@data$Tumor_Sample_Barcode), collapse = '\t'))
  ), f
)
print('PEDS Nantomics Germline Mutation case list file completed.')
close(f)
