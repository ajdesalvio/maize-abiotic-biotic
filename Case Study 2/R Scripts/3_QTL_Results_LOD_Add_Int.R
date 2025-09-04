library(tidyverse)
library(data.table)
library(qtl2)
library(ggplot2)
library(jsonlite)
library(purrr)
library(stringr)

# Path to data (replace with your own path)
data_path <- 'C:/Users/aaron.desalvio/Downloads/Abiotic_Data/'

# Obtain .rds file paths for each QTL analysis that was run
qtl_names <- c('Additive', 'Interactive')
qtl_rds <- c('Additive.RCC.QTL.Outputs.rds', 'Interactive.RCC.QTL.Outputs.rds')

# Loop through scan1 results and extract LOD scores and peaks with intervals
scan1_results_list <- list()
lod_peaks <- list()
permutation_list <- list()

for (i in 1:length(qtl_rds)) {
  # Loop variables
  qtl_name_i <- qtl_names[i]
  rds_file_i <- qtl_rds[i]
  
  ## LOD scores ##
  # Read the rds file
  message('Reading .rds file: ', rds_file_i)
  qtl_results <- read_rds(paste0(data_path, rds_file_i))
  # Extract the LOD score data frame
  LOD_df <- qtl_results[[1]]$scan.result %>% as.data.frame()
  LOD_df$QTL.Analysis.Name <- qtl_name_i
  # Make columns for chromosome and position
  LOD_df$marker <- rownames(LOD_df)
  LOD_df <- LOD_df %>%
    separate(marker, into = c('chr', 'pos'), sep = '_') %>%
    mutate(chr = str_remove(chr, 'S'), pos = as.numeric(pos))
  # Transform wide to long
  pheno_cols <- colnames(qtl_results[[1]]$pheno.loop.final)
  LOD_df_long <- LOD_df %>%
    pivot_longer(cols = all_of(pheno_cols), names_to = 'Phenotype', values_to = 'LOD')
  # Save results in list
  scan1_results_list[[qtl_name_i]] <- LOD_df_long
  
  ## Permutation test results ##
  perm_result <- as.data.frame(summary(qtl_results[[1]]$permutation.result))
  perm_result$QTL.Analysis.Name <- qtl_name_i
  perm_result_long <- perm_result %>%
    pivot_longer(cols = all_of(pheno_cols), names_to = 'Phenotype', values_to = 'Threshold')
  # Save results in list
  permutation_list[[qtl_name_i]] <- perm_result_long
  
  ## LOD peaks with intervals ##
  LOD_peaks_df <- qtl_results[[1]]$LOD.peaks %>% as.data.frame()
  if (nrow(LOD_peaks_df) == 0) {
    message('LOD peaks data frame had 0 rows for ', qtl_name_i)
    next
  }
  LOD_peaks_df$QTL.Analysis.Name <- qtl_name_i
  # Save results in list
  lod_peaks[[qtl_name_i]] <- LOD_peaks_df
}

# Stack the LOD scores
Lod_combined_long <- rbindlist(scan1_results_list, fill = TRUE)
# Stack the permutation test results
perm_combined_long <- rbindlist(permutation_list, fill = TRUE)
# Stack the LOD peaks
lod_peaks_long <- rbindlist(lod_peaks, fill = TRUE)

# Export the LOD scores, LOD peaks with intervals, and permutation cutoffs
#fwrite(Lod_combined_long, paste0(data_path, 'QTL.Results.LOD.Add.Int.csv'))
#fwrite(perm_combined_long, paste0(data_path, 'QTL.Results.Perm.Add.Int.csv'))
#fwrite(lod_peaks_long, paste0(data_path, 'QTL.Results.LODPeaks.Add.Int.csv'))