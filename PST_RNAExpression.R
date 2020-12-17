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
      df$filename <- basename(i)
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

expression_files <- list.files(pattern = '*.rgel', recursive = TRUE)

expression_df <- combine_data(expression_files)
expression_df$sample <- substr(expression_df$filename, 1, regexpr("\\.", expression_df$filename)-1)

expression_transpose <- data.frame(matrix(ncol = length(unique(expression_df$sample))+2,
                                          nrow = length(unique(expression_df$gene_name))))
colnames(expression_transpose) <- c("Hugo_Symbol", "Entrez_Gene_Id", unique(expression_df$sample))

expression_transpose$Hugo_Symbol <- unique(expression_df$gene_name)

# For each row in lt_transpose, move values into tests dataframe
for (row in seq_len(nrow(expression_df))) {
  expvalue <- expression_df[row, "expression"]
  gene <- expression_df[row, "gene_name"]
  sample <- expression_df[row, "sample"]
  expression_transpose[expression_transpose$Hugo_Symbol == gene, sample] <- expvalue
}

Entrez_Gene_Ids <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"), filters = 'hgnc_symbol',
                         values=expression_transpose$Hugo_Symbol, ensembl)
colnames(Entrez_Gene_Ids)[2] <- "Entrez_Gene_Id"

expression_final <- left_join(expression_transpose, Entrez_Gene_Ids, by = c('Hugo_Symbol' = 'hgnc_symbol'))

code <- "c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS"
setwd(code)
samples <- read.csv("data_clinical_sample_formatted.txt", sep = '\t', header = FALSE) %>% dplyr::select('V2')

expression_exists <- unique(expression_final[(expression_final$Tumor_Sample_Barcode %in% samples$V2), ])
missing <- unique(expression_final[!(expression_final$Tumor_Sample_Barcode %in% samples$V2), "Tumor_Sample_Barcode" ])

expressionFile <- "data_rna_expression.txt"

write.table(expression_final, expressionFile, sep="\t", col.names = TRUE, row.names = FALSE,
            quote = FALSE, append = FALSE, na = "NA")

fusionCL <- ("case_lists/cases_sv.txt")

f <- file(fusionCL)
writeLines(c(
 "cancer_study_identifier: PST_PEDS_2020",
 "stable_id: PST_PEDS_2020_sv",
 "case_list_name: RNA Fusion",
 "case_list_description: RNA Fusion",
 paste("case_list_ids: ", paste(unique(fusion$Tumor_Sample_Barcode), collapse = '\t'))
), f
)
close(f)

FusionMetaFile <- "meta_fusion.txt"

f <- file(FusionMetaFile)
writeLines(c(
  "cancer_study_identifier: PST_PEDS_2020",
  "genetic_alteration_type: FUSION",
  "datatype: FUSION",
  "stable_id: FUSION",
  "show_profile_in_analysis_tab: true",
  "profile_name: RNA Fusion",
  "profile_description: RNA Fusion",
  paste("data_filename: ", fusionFile)
), f
)
close(f)
print("RNA Fusion metafile completed")


copynumber_files <- list.files(pattern = '*.copynumber.csv', recursive = TRUE)
somatic_files <- list.files(pattern = '*.somatic.vcf$', recursive = TRUE)
germline_files <- list.files(pattern = '*.germline.vcf$', recursive = TRUE)

