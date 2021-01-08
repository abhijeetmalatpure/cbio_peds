library(dplyr)
library(biomaRt)
library(tidyr)
library(stringr)
library(data.table)

combine_data <- function (filenames) {
  merge_df <- data.table()

  for (i in filenames) {
    # Set df_name to the vendor name
    df <- read.csv(i, sep = ",", header = TRUE, na.strings = c("N/A", "", "unavailable"))

    if(nrow(df)) {
      df$sample <- substr(basename(i), 1, regexpr("\\.", basename(i))-1)
      df$path <- i
      if (nrow(merge_df) == 0) {
          merge_df <- df
      }
      else {
        if(identical(colnames(merge_df), colnames(df))) {
        merge_df <- rbind(merge_df, df)
        }
        else print("Headers don't match")
      }# Check for headers match
    }
  }
  return(merge_df)
}

ashion <- "c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS/data/ashion"
foundation <- "c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS/data/foundation"
setwd(ashion)

ensembl <- useMart(host = 'grch37.ensembl.org',
                   biomart = 'ENSEMBL_MART_ENSEMBL',
                   dataset = 'hsapiens_gene_ensembl')

copynumber_files <- list.files(pattern = '*.copynumber.csv', recursive = TRUE)

cn_df <- combine_data(copynumber_files)

# Remove chr prefix from chromosomes
cn_df$chromosome <- gsub("^chr", replacement = "", x= cn_df$chromosome, ignore.case = FALSE)

for(row in row.names(cn_df)) {
  cn_df$log2[as.integer(row)] <- (fromJSON(gsub('\'', '"', cn_df$attributes[as.integer(row)]))$LOG2FC)
}

cn_transpose <- data.frame(matrix(ncol = length(unique(cn_df$sample))+1,
                                          nrow = length(unique(cn_df$gene))))
colnames(cn_transpose) <- c("Hugo_Symbol", unique(cn_df$sample))

cn_transpose$Hugo_Symbol <- unique(cn_df$gene)

for (row in seq_len(nrow(cn_df))) {
  expvalue <- cn_df[row, "copy_number"]
  gene <- cn_df[row, "gene"]
  sample <- cn_df[row, "sample"]
  cn_transpose[cn_transpose$Hugo_Symbol == gene, sample] <- expvalue
}

Entrez_Gene_Ids <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = 'hgnc_symbol',
                         values=cn_transpose$Hugo_Symbol, ensembl)
colnames(Entrez_Gene_Ids) <- c("Hugo_Symbol", "Entrez_Gene_Id")

cn_final <- left_join(cn_transpose, Entrez_Gene_Ids, by = 'Hugo_Symbol')

cn_final <- cn_final[, (colnames(cn_final) %in% c("Hugo_Symbol") == FALSE)]


code <- "c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS"
setwd(code)

cnFile <- "data_copynumber_continuous.txt"

write.table(cn_final %>% dplyr::select(c('Entrez_Gene_Id', everything())), cnFile, sep="\t", col.names = TRUE, row.names = FALSE,
            quote = FALSE, append = FALSE, na = "NA")

cnCL <- ("case_lists/cases_cna.txt")

f <- file(cnCL)
writeLines(c(
 "cancer_study_identifier: PST_PEDS_2020",
 "stable_id: PST_PEDS_2020_cna",
 "case_list_name: Copy Number Alterations",
 "case_list_description: Copy Number Alterations  [Continuous]",
 "case_list_category: all_cases_with_cna_data",
 paste("case_list_ids: ", paste(unique(cn_df$sample), collapse = '\t'))
), f
)
close(f)

cnMetaFile <- "meta_copynumber_continuous.txt"

f <- file(cnMetaFile)
writeLines(c(
  "cancer_study_identifier: PST_PEDS_2020",
  "genetic_alteration_type: COPY_NUMBER_ALTERATION",
  "datatype: CONTINUOUS",
  "stable_id: linear_CNA",
  "show_profile_in_analysis_tab: false",
  "profile_name: Copy Number Alterations",
  "profile_description: Copy Number Alterations",
  paste("data_filename: ", cnFile)
), f
)
close(f)
print("CNA metafile completed")
