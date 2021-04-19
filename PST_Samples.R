library(dplyr)
library(tidyr)
library(stringr)
library(data.table)


setwd("c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS")
smpl <- read.csv("data_clinical_sample_formatted.txt", sep = "\t", header = FALSE)
smpl$V6 <- str_trim(smpl$V6)
smpl[smpl$V6 == 'ASHION', ]
a_smpl <- smpl[smpl$V6 == 'ASHION', ]
n_smpl <- smpl[smpl$V6 == 'NANTOMICS', ]
f_smpl <- smpl[smpl$V6 == 'FOUNDATION', ]

a_map <- read.csv("data/ashion/ashion_vcf_metadata.csv", header=TRUE, sep = "\t")
final <- smpl

a_patientid <- left_join(unique(a_map[a_map$vcftype=="somatic", c("tumor", "normal", "vcftype")]), a_smpl[,c("V1","V2")], by = c("tumor"="V2"))
colnames(a_patientid) <- c("tumor", "V2", "V6", "V1")
a_patientid <- a_patientid[, c("V1", "V2", "V6")]
a_patientid$V3 <- ""
a_patientid$V4 <- ""
a_patientid$V5 <- ""
rbind(final, a_patientid)
final <- rbind(final, a_patientid)
final[final$V6=="somatic", "V6"] <- "ASHION"
unique(final$V6)

final[393:432, c("V3", "V4", "V5")] <- NA
write.table(final, "data_clinical_patient_formatted.txt", sep="\t", col.names = FALSE, row.names = FALSE, append=FALSE, quote = FALSE, na = "")
final <- rbind(final[1,], unique(final[2:nrow(final)]))
final <- rbind(final[1,], unique(final[2:nrow(final),]))
f_map <- read.csv("data/foundation/vcf_metadata_foundation_somatic.tsv", sep="\t", header = TRUE)
setdiff(f_map$sample_id, final$V2)
n_map <- read.csv("data/nantomics/vcf_metadata.csv", sep=",", header = TRUE)
setdiff(n_map$tumor_id, final$V2)
setdiff(n_map$normal_id, final$V2)
write.table(final, "data_clinical_sample_formatted.txt", sep="\t", col.names = FALSE, row.names = FALSE, append=FALSE, quote = FALSE, na = "")
germline <- c(unique(a_map$normal), unique(n_map$normal_id))
somatic <- c(unique(a_map$tumor), unique(n_map$tumor_id), unique(f_map$sample_id))


mutationCL <- ("case_lists/cases_germline.txt")

f <- file(mutationCL)
writeLines(c(
  "cancer_study_identifier: PST_PEDS_2020",
  "stable_id: PST_PEDS_2020_germline",
  "case_list_name: Germline Samples",
  "case_list_description: Germline samples [ASHION, NANTOMICS]",
  paste("case_list_ids: ", paste(unique(germ), collapse = '\t'))
  ), f
)
print('PEDS Germline samples case list file completed.')
close(f)

mutationCL <- ("case_lists/cases_somatic.txt")

f <- file(mutationCL)
writeLines(c(
  "cancer_study_identifier: PST_PEDS_2020",
  "stable_id: PST_PEDS_2020_somatic",
  "case_list_name: Somatic Samples",
  "case_list_description: Somatic samples [ASHION, FOUNDATION, NANTOMICS]",
  paste("case_list_ids: ", paste(unique(somatic), collapse = '\t'))
  ), f
)
print('PEDS Somatic Mutation case list file completed.')
close(f)