# Load necessary packages
library(dplyr)
library(readr)

directory = '/re_gecip/hearing_and_sight/dgreen/albinism_digenic/hgmd_variants'


# TYR function
merge_files_for_tyr <- function(gene_name) {
  
  # Define file paths
  file1 <- paste0(directory, "/CNV_outputs/CNV_SV_", gene_name, "_counts.csv")
  file2 <- paste0(directory, "/outputs/", gene_name, "_plink_out_processed")
  file3 <- paste0(directory, "/status-by-ID-full-genotypes.txt")  # Replace with actual path

  # Read csv files
  df1 <- read.csv(file1, sep = '\t')
  df2 <- read.csv(file2, sep = '\t')
  df3 <- read.csv(file3, sep = '\t')

  # Merge the dataframes based on 'platekey' column
  merged_df <- df1 %>%
    full_join(df2, by = "plate_key") %>%
    full_join(df3, by = "plate_key")  # merge the third file
  
  # Replace NA values with 0
  merged_df[is.na(merged_df)] <- 0
  
  # Sum across the 'het_count' and 'hom_count' columns
  merged_df <- merged_df %>%
    rowwise() %>%
    mutate(het_count = sum(c_across(starts_with('het_count'))),
           hom_count = sum(c_across(starts_with('hom_count')))) %>%
    ungroup()

  # Convert 'snp_402', 'snp_192', and 'promotor' columns
  for (snp_col in c('snp_402', 'snp_192', 'promotor')) {
    merged_df <- merged_df %>%
      mutate("{snp_col}" := str_count(get(snp_col), "1"))
  }
  
  write_csv(merged_df, paste0(directory, '/outputs/', gene_name, "_processed_corrected.csv"))

  return(merged_df)
}



# OCA2 function
merge_files_for_gene <- function(gene_name) {
  
  # Define file paths
  file1 <- paste0(directory, "/CNV_outputs/CNV_SV_", gene_name, "_counts.csv")
  file2 <- paste0(directory, "/outputs/", gene_name, "_plink_out_processed")

  # Read csv files
  df1 <- read.csv(file1, sep = '\t')
  df2 <- read.csv(file2, sep = '\t')

  # Merge the dataframes based on 'platekey' column
  merged_df <- df1 %>%
    full_join(df2, by = "plate_key")
  
  # Replace NA values with 0
  merged_df[is.na(merged_df)] <- 0
  
  # Sum across the 'het_count' and 'hom_count' columns
  merged_df <- merged_df %>%
    rowwise() %>%
    mutate(het_count = sum(c_across(starts_with('het_count'))),
           hom_count = sum(c_across(starts_with('hom_count')))) %>%
    ungroup()
  
  write_csv(merged_df, paste0(directory, '/outputs/', gene_name, "_processed_corrected.csv"))

  return(merged_df)
}



# Call the function for each gene
oca2_df <- merge_files_for_gene("OCA2")
tyr_df <- merge_files_for_tyr("TYR")
