library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(maftools)


combine_maf_files <- function (filenames) {
  maf_files <- list()
  failed_maf <- list()
  x <- 1
  for (i in seq_along(filenames)) {
      maf <- read_maf_file(filenames[i])#, vc_nonSyn = vcNames)
      if (!is.na(maf)) {
        maf_files[x] <- maf
        x <- x + 1 }
      else {
        failed_maf <- append(failed_maf, filenames[i])
      }
  }
  print("=======================================================================")
  return(list("merged"= merge_mafs(maf_files), "failed" = failed_maf))
}


read_maf_file <- function(x) {
  maf <- data.table::fread(x)
  vcNames <- unique(maf$Variant_Classification)
  tryCatch(read.maf(x, vc_nonSyn = vcNames),
           error = function(x) {print(paste(x, ": Error reading file. Skipping")); NaN})
}

#mafdir <- "c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS/data/ashion/enhanced"
#mafdir <- "c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS/data/nantomics/maf"
mafdir <- "/N/project/phi_asha_archive/peds_pst/ashion/maf/somatic"
setwd(mafdir)

mafs <- list.files(mafdir, pattern = "*.maf$", recursive = FALSE)
somaticmafs <- list.files(mafdir, pattern = "*\\.somatic.maf$", recursive = FALSE)
germlinemafs <- list.files(mafdir, pattern = "*\\.germline.maf$", recursive = FALSE)

combined_mafs <- combine_maf_files(mafs)
somatic_combined <- combine_maf_files(somaticmafs)
germline_combined <- combine_maf_files(germlinemafs)

# Save merged MAF file and failed list of MAFs
peds_somatic_merged <- somatic_combined$merged
peds_somatic_failed <- somatic_combined$failed

peds_somatic_merged@data$Mutation_Status <- as.character(peds_somatic_merged@data$Mutation_Status)
peds_somatic_merged@data$Mutation_Status <- "Somatic" # Or "Germline" for germline
# Write merged MAF file
write.table(as.data.frame(peds_somatic_merged@data), "peds_ashion_somatic.maf", col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE, append = FALSE, na = "")

# Save failed maf file list
write.table(as.data.frame(peds_somatic_failed), "peds_foundation_somatic_failed.txt", col.names = FALSE, row.names = FALSE, sep = "\n",quote = FALSE, append = FALSE, na = "NA")

dataset <- "c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS"
setwd(dataset)

# Add samples from germline to samples file. Match patient with the same as it's somatic sample
vcf_metadata <- read.csv('vcf_metadata.csv', sep = '\t', header = TRUE) %>% filter(vcftype == 'germline') %>% dplyr::select(c(sample_id, normal))

samples <- read.csv('data_clinical_sample_formatted.txt', sep = '\t', header = FALSE)

vcf_metadata <- left_join(vcf_metadata, (samples %>% dplyr::select(c(V1, V2))), by = c("sample_id"="V2"))
vcf_metadata$V6 <- 'ASHION'
colnames(vcf_metadata)[1] <- "V2"
vcf_metadata[, c('V3', 'V4', 'V5')] <- NA

samples_final <- rbind(samples_final, vcf_metadata %>% dplyr::select(c('V1', 'V2', 'V3', 'V4', 'V5', 'V6')))
samples_final <- unique(samples_final)
write.table(samples_final, "data_clinical_sample_formatted.csv", col.names = FALSE,
            row.names = FALSE, sep = "\t", quote = FALSE, append = FALSE, na = "NA")


# write.mafSummary(somatic_combined, "data_mutation_somatic_enhanced_combined.maf")
mutationFile <- "data_mutations_combined_enhanced.maf"

write.table(as.data.frame(som_data), "data_somatic_mutations_combined_enhanced.maf", col.names = TRUE,
            row.names = FALSE, sep = "\t", quote = FALSE, append = FALSE, na = "NA")

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
