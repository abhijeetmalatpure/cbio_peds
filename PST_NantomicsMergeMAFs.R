library(dplyr)
#library(tidyr)
library(stringr)
library(data.table)
library(maftools)
# library(doParallel)
#
# cores <- detectCores()
# cl <- makeCluster(cores[1]/2)
# registerDoParallel(cl)

combine_maf_files <- function (filenames) {
  maf_files <- list()
  x <- 1
  for (i in seq_along(filenames)) {
    # Set maf_files[i] to the content of filenames[i]
    #print(paste(i, filenames[i]))
    print(paste(i, ": ", filenames[i]))
    maf_files[x] <- read.maf(filenames[i])#, vc_nonSyn = vcNames)
    x <- x + 1
  }
  return(merge_mafs(maf_files))
}

mafdir <- "/N/project/phi_asha_archive/peds_pst/nantomics/vcf2mafConversion/maf"
mafdir <- "c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS/data/nantomics/maf"
setwd(mafdir)

# Add samples from germline to samples file. Match patient with the same as it's somatic sample
# vcf_metadata <- read.csv('vcf_metadata.csv', sep = ',', header = TRUE)
#mafs <- list.files(mafdir, pattern = "*.maf$", recursive = FALSE)
#somatic_mafs <- list.files(mafdir, pattern = "*somatic.maf$", recursive = FALSE)
germline_mafs <- list.files(mafdir, pattern = "data_mutations_g.*", recursive = FALSE)

# maf_files <- c("data_mutations_somatic_combined.maf",
#                "data_mutations_germline_combined.maf",
#                "data_mutations_combined_enhanced.maf")

#combined_mafs <- combine_maf_files(mafs)
#somatic_combined_mafs <- combine_maf_files(somatic_mafs)
germline_combined_mafs <- combine_maf_files(germline_mafs)
#combined_mafs <- combine_maf_files(maf_files)

# write.mafSummary(somatic_combined, "data_mutation_somatic_enhanced_combined.maf")
mutationFile <- "data_mutations_g6.maf"
print(paste("Writing to a final maf file:", mutationFile))


write.table(as.data.frame(germline_combined_mafs@data), mutationFile,
            col.names = TRUE, row.names = FALSE, sep = "\t",
            quote = FALSE, append = FALSE, na = "NA")

mutationCL <- ("cases_sequenced_g6.txt")

f <- file(mutationCL)
writeLines(c(
  "cancer_study_identifier: PST_PEDS_2020",
  "stable_id: PST_PEDS_2020_sequenced",
  "case_list_name: Sequenced Tumors",
  "case_list_description: All Sequenced Tumors",
  paste("case_list_ids: ", paste(unique(germline_combined_mafs@data$Tumor_Sample_Barcode), collapse = '\t')),
  paste("germline_case_list_ids: ", paste(unique(germline_combined_mafs@data$Matched_Norm_Sample_Barcode), collapse = '\t'))
  ), f
)
print('Mutation case list g5 file completed.')
close(f)

# mutationMetafile <- "meta_mutations_combined.txt"
#
# f <- file(mutationMetafile)
# writeLines(c(
#  "cancer_study_identifier: PST_PEDS_2020",
#  "genetic_alteration_type: MUTATION_EXTENDED",
#  "datatype: MAF",
#  "stable_id: mutations",
#  "show_profile_in_analysis_tab: true",
#  "profile_name: Germline and somatic mutation data",
#  "profile_description: Mutations",
#   paste("data_filename: ", mutationFile)
# ), f
# )
# close(f)
# print("Mutations metafile completed")


#
# # Add samples from germline to samples file. Match patient with the same as it's somatic sample
# setwd(mafdir)
# vcf_metadata <- read.csv('vcf_metadata.csv', sep = ',', header = TRUE)
#
# dataset <- "c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS"
# setwd(dataset)
#
# samples <- read.csv('data_clinical_sample_formatted.txt', sep = '\t', header = FALSE)
#
# vcf_metadata_merged <- left_join(vcf_metadata, samples %>% dplyr::select(c(V1,V2, V6)), by = c("OrderBarcode"  = "V2"))
# vcf_metadata_merged$V3 <- NA
# vcf_metadata_merged$V4 <- NA
# vcf_metadata_merged$V5 <- NA
#
# normal_samples <- vcf_metadata_merged %>% dplyr::select(c(V1, normal_id, V3, V4, V5, V6))
# colnames(normal_samples)[2] <- 'V2'
#
# tumor_samples <- vcf_metadata_merged %>% dplyr::select(c(V1, tumor_id, V3, V4, V5, V6))
# colnames(tumor_samples)[2] <- 'V2'
#
# samples_final <- do.call("rbind", list(samples, normal_samples %>% dplyr::filter(!is.na(V1)), tumor_samples %>% dplyr::filter(!is.na(V1))))
# samples_final <- unique(samples_final)
# write.table(samples_final, "data_clinical_sample_formatted.txt", col.names = FALSE,
#             row.names = FALSE, sep = "\t", quote = FALSE, append = FALSE, na = "NA")
