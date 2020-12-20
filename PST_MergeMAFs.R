library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(maftools)


combine_maf_files <- function (filenames) {
  maf_files <- list()
  for (i in seq_along(filenames)) {
    # Set maf_files[i] to the content of filenames[i]
    maf_files[i] <- read.maf(filenames[i])
  }
  return(merge_mafs(maf_files))
}

mafdir <- "c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS/data/ashion/enhanced"
setwd(mafdir)

mafs <- list.files(mafdir, pattern = "*.maf", recursive = FALSE)

combined_mafs <- combine_maf_files(mafs)

dataset <- "c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS"
setwd(dataset)

# write.mafSummary(somatic_combined, "data_mutation_somatic_enhanced_combined.maf")
mutationFile <- "data_mutations_combined_enhanced.maf"

write.table(as.data.frame(combined_mafs@data), mutationFile,
            col.names = TRUE, row.names = FALSE, sep = "\t",
            quote = FALSE, append = FALSE, na = "NA")

mutationCL <- ("case_lists/cases_sequenced.txt")

f <- file(mutationCL)
writeLines(c(
  "cancer_study_identifier: PST_PEDS_2020",
  "stable_id: PST_PEDS_2020_sequenced",
  "case_list_name: Sequenced Tumors",
  "case_list_description: All Sequenced Tumors",
  paste("case_list_ids: ", paste(unique(combined_mafs@data$Tumor_Sample_Barcode), collapse = '\t'))
  ), f
)
print('Mutation case list file completed.')
close(f)

mutationMetafile <- "meta_mutations_combined.txt"

f <- file(mutationMetafile)
writeLines(c(
 "cancer_study_identifier: PST_PEDS_2020",
 "genetic_alteration_type: MUTATION_EXTENDED",
 "datatype: MAF",
 "stable_id: mutations",
 "show_profile_in_analysis_tab: true",
 "profile_name: Germline and somatic mutation data",
 "profile_description: Mutations",
  paste("data_filename: ", mutationFile)
), f
)
close(f)
print("Mutations metafile completed")
