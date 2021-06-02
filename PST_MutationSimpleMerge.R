# Title     : Merge MAF files into one using bindlist
# Objective : TODO
# Created by: abhmalat
# Created on: 4/22/2021

library(dplyr)
library(stringr)
library(data.table)
library(maftools)

combine_maf_files <- function(filenames) {
  final_maf <- c()
  for(f in filenames) {
      print(paste("Reading", f))
      maf <- read.csv(f, sep="\t", header = TRUE, comment.char = "#")
    if(length(final_maf)==0) {
      final_maf <- maf
    }
    else {
      final_maf <- rbindlist(list(final_maf, maf), fill = TRUE)
      remove(maf)
    }
    gc()
  }
  return(final_maf)
}

mafdir <- "/N/project/phi_asha_archive/peds_pst/foundation/maf/somatic"
setwd(mafdir)

#germlinemafs <- list.files(mafdir, pattern = "*\\.germline.maf$", recursive = FALSE)
pattrn <- "*\\.somatic.maf"#"indels"#"snvs"# "somaticSV"
mafs <- file.path(list.files(mafdir, pattern = pattrn, recursive = FALSE))
mafs
# Call combine_maf_files
processed <- combine_maf_files(mafs)

gc()

processed$Mutation_Status <- as.character(processed$Mutation_Status)
processed$Mutation_Status <- "Somatic" # Or "Germline" for germline
processed$Center <- as.character(processed$Center)
processed$Center <- "FOUNDATION"

#processed$NCBI_Build <- "GRCh37"

write.table(as.data.frame(processed),
            paste0("peds_foundation_somatic_simple.maf"), col.names = TRUE, row.names = FALSE, sep = "\t",quote = FALSE, append = FALSE, na = "")
