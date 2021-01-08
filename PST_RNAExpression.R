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

expression_transpose <- data.frame(matrix(ncol = length(unique(expression_df$sample))+1,
                                          nrow = length(unique(expression_df$gene_name))))
colnames(expression_transpose) <- c("Hugo_Symbol", unique(expression_df$sample))

# Remove sample ids that do not occur in the main samples file
# code <- "c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS"
# setwd(code)
# samples <- read.csv("data_clinical_sample_formatted.txt", sep = '\t', header = FALSE) %>% dplyr::select('V2')
#
# expression_exists <- unique(expression_df[(expression_df$sample %in% samples$V2), ])
# missing <- unique(expression_df[!(expression_df$sample %in% samples$V2), "sample" ])

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
colnames(Entrez_Gene_Ids) <- c("Hugo_Symbol", "Entrez_Gene_Id")

# expression_transpose <- expression_transpose[!names(expression_transpose) %in% c("Entrez_Gene_Id",
#                                                                                  "C051_0043_035165_T1_K1ID2_ps20201129101455",
#                                                                                  "C051_0042_034525_T1_K1ID2_ps20201121222158")]
expression_final <- left_join(expression_transpose, Entrez_Gene_Ids, by = 'Hugo_Symbol')

newexp <- expression_final[!is.na(expression_final$Entrez_Gene_Id), ] %>% dplyr::select(-c("Hugo_Symbol"))
dataset <- "c:/Users/abhmalat/OneDrive - Indiana University/cBio_PEDS"
setwd(dataset)

expressionFile <- "data_rna_expression.txt"

write.table(newexp %>% dplyr::select(Entrez_Gene_Id, everything()), expressionFile, sep="\t",
            col.names = TRUE, row.names = FALSE,
            quote = FALSE, append = FALSE, na = "NA")


expressionCL <- ("case_lists/cases_rna_seq_mrna.txt")

f <- file(expressionCL)
writeLines(c(
 "cancer_study_identifier: PST_PEDS_2020",
 "stable_id: PST_PEDS_2020_rna_seq_mrna",
 "case_list_name: RNA Expression samples",
 "case_list_description: RNA expression samples [Continuous]",
 paste("case_list_ids: ", paste(unique(expression_df$sample), collapse = '\t'))
), f
)
close(f)

expressionMetafile <- "meta_rna_expression.txt"

f <- file(expressionMetafile)
writeLines(c(
 "cancer_study_identifier: PST_PEDS_2020",
 "genetic_alteration_type: MRNA_EXPRESSION",
 "datatype: CONTINUOUS",
 "stable_id: rna_seq_mrna",
 "show_profile_in_analysis_tab: false",
 "profile_name: RNA Expression",
 "profile_description: RNA Expression [Continuous]",
  paste("data_filename: ", expressionFile)
), f
)
close(f)
print("RNA Expression metafile completed")


copynumber_files <- list.files(pattern = '*.copynumber.csv', recursive = TRUE)
somatic_files <- list.files(pattern = '*.somatic.vcf$', recursive = TRUE)
germline_files <- list.files(pattern = '*.germline.vcf$', recursive = TRUE)

