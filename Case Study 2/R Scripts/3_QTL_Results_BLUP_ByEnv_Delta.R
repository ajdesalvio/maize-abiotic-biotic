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
qtl_names <- c('Delta', 'TXH2.2021', 'TXH3.2021')
qtl_rds <- c('Delta.RCC.QTL.Outputs.rds', 'TXH2.2021.QTL.Outputs.rds', 'TXH3.2021.QTL.Outputs.rds')

for (i in 1:length(qtl_rds)) {
  # Loop variables
  qtl_name_i <- qtl_names[i]
  rds_file_i <- qtl_rds[i]
  
  # Read the .rds file
  results <- read_rds(paste0(data_path, rds_file_i))
  my_res <- results[[1]]

  # Extract the genome scan results
  LOD_peaks <- my_res$LOD.peaks
  if (nrow(LOD_peaks) == 0) {
    next  # Skip this iteration if LOD_peaks has zero rows
  }

  scan1blup_list <- list()
  # Loop to run all scan1blup instances for this Env/Tester combination
  for (j in 1:nrow(LOD_peaks)) {
    # Extract the chromosome number and trait name
    chr <- LOD_peaks$chr[j]
    trait <- LOD_peaks$lodcolumn[j]
    # Run the scan1blup() function for this chr/trait combination
    message(paste0('scan1blup() being run for ', qtl_name_i, ', Chr ', chr, ', trait ', trait))
    coef <- scan1blup(genoprobs = my_res$probs.final[,chr],
                      pheno = my_res$pheno.loop.final[,trait, drop = FALSE],
                      kinship = my_res$kinship[[chr]],
                      cores = 0)
    # Save the scan1blup result as a dataframe
    scan1blup_results <- as.data.frame(coef)
    scan1blup_results$QTL.Analysis.Name <- qtl_name_i
    scan1blup_results$Trait <- trait
    scan1blup_results$marker <- rownames(scan1blup_results)
    # Split the marker name (e.g., "S1_39898") into chr and pos
    scan1blup_results <- scan1blup_results %>%
      separate(marker, into = c('chr', 'pos'), sep = '_') %>%
      mutate(chr = str_remove(chr, 'S'), pos = as.numeric(pos))
    # Save in list
    scan1blup_list[[paste(qtl_name_i, 'Chr', chr, trait, sep = '.')]] <- scan1blup_results
  }
}

# Aggregate results from scan1 and scan1blup
scan1blup_results_wide <- rbindlist(scan1blup_list, fill = TRUE)
scan1blup_results_long <- scan1blup_results_wide %>%
  pivot_longer(cols = AA:intercept, names_to = 'DH6.Parent', values_to = 'BLUP')

# Export CSV files
fwrite(scan1blup_results_wide, paste0(data_path, 'QTL.Results.BLUPs.ByEnv.Delta.wide.csv'))
fwrite(scan1blup_results_long, paste0(data_path, 'QTL.Results.BLUPs.ByEnv.Delta.long.csv'))
