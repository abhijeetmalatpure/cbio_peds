# Title     : Merge MAF files into one using bindlist. Separate script for Nantomics large files
# Objective : TODO
# Created by: abhmalat
# Created on: 4/22/2021

library(dplyr)
library(stringr)
library(data.table)
library(maftools)

# Get the header row for the maf file
get_header <- function(file) {
  print(paste("Reading header for ", file))
  maf <- read.csv(file, sep="\t", header = TRUE, comment.char = "#", nrows = 2) %>%
    dplyr::select(-c("all_effects"))
  #print(colnames(maf))
  return(colnames(maf))
}


combine_maf_files <- function(filenames, header) {
  # Empty maf file with header = final_header.
  empty_maf <- transpose(data.table(header))
  colnames(empty_maf) <- header
  empty_maf <- empty_maf[-1,]

  for(f in seq_along(filenames)) {
    print(paste("Reading ", filenames[f]))
    maf <- read.csv(filenames[f], sep="\t", header = TRUE, comment.char = "#") %>%
    dplyr::filter(Hugo_Symbol != "Unknown")
    # maf$Mutation_Status <- as.character(maf$Mutation_Status)
    # maf$Mutation_Status <- "Germline" # Or "Germline" for germline
    # maf$Center <- as.character(maf$Center)
    # maf$Center <- "NANTOMICS"

    #  rbindlist each maf file with this empty maf
    combined_maf <- rbindlist(list(empty_maf, maf), fill = TRUE)

    if(f==1) {
      # If first file, header=True, append = False
      write.table(as.data.frame(combined_maf), paste0("peds_nant_germline_simple2.maf"),
                  col.names = TRUE,  row.names = FALSE, sep = "\t", quote = FALSE, append = FALSE, na = "")
    }
    else {
      # For others, append to existing final file with header=False, append=True
      write.table(as.data.frame(combined_maf), paste0("peds_nant_germline_simple2.maf"),
                  col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE, append = TRUE,  na = "")
    }
    remove(maf)
    remove(combined_maf)
    gc()
  }
  print("File created.")
}

#mafdir <- "/N/project/phi_asha_archive/peds_pst/nantomics/vcf2mafConversion/maf/somatic"
mafdir <- "/N/project/phi_asha_archive/peds_pst/nantomics/vcf2mafConversion/maf/germline"
setwd(mafdir)

#germlinemafs <- list.files(mafdir, pattern = "*\\.germline.maf$", recursive = FALSE)
pattrn <- "*\\.germline.maf"#"indels"#"snvs"# "somaticSV"
mafs <- file.path(list.files(mafdir, pattern = pattrn, recursive = TRUE))
mafs

# Get headers for all maf files into one list. final header = unique elements of list
header <- NULL
for(maf in mafs) {
  header <- c(header, get_header(maf))
}
header <- unique(header)
print(header)

# Combine maf files
combine_maf_files(mafs, header)