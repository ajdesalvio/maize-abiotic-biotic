library(tidyverse)
library(data.table)
library(qtl2)
library(ggplot2)
library(jsonlite)
library(purrr)
library(stringr)

# Steps used to prepare JSON files, subset genomic data files, and subset response value files (FPCs)
# Path to data (replace with your own path)
data_path <- 'C:/Users/aaron.desalvio/Downloads/Abiotic_Data/'
# Download the supplementary data from Michel et al., 2022 (https://doi.org/10.1093/genetics/iyac063)
# The CSV files you need all start with "cm.mbp.ss.100k.sites"
# Copy all of these CSV files into the data_path folder!

# First, we need to create a phenotype file that will replace ss.blue.phenotypes.csv from the Michel supplementary data. In our
# case, we need a column for Inbred, then the remaining columns will be FPC#_VI.

# Load in the genomic data from the supplementary data so we can get a vector of unique inbred names
supp.genomic <- fread(paste0(data_path, 'cm.mbp.ss.100k.sites.SSpopulation_geno.csv')) %>% as.data.frame()
supp.inbreds <- supp.genomic$lines

# Load the FPCs and put them into wide format
VI <- fread(paste0(data_path, 'TXH2.TXH3.2021.RCC.FPCA.Scores.csv')) %>% as.data.frame()
# Sort by Environment, Pedigree
VI <- VI %>% arrange(Env, Pedigree)
# Although we ran FPCA including the commercial check hybrids, we only want the WI SS MAGIC inbreds
WI.SS <- unique(VI$Pedigree[grep('W10004', VI$Pedigree)])
# Remove pedigrees not contained in the overlap vector
VI <- VI %>% filter(Pedigree %in% WI.SS)
# Unique environments
envs <- sort(unique(VI$Env))

# Transform the VI FPCs to wide format (only retaining FPCs needed to explain 95% of FVE)
VI.wide <- VI %>%
  pivot_wider(
    id_cols = c(Pedigree.Env, Env, Pedigree),
    names_from = Vegetation.Index,
    names_glue = '{.value}_{Vegetation.Index}',
    values_from = c(FPC1:FPC6),
    values_fill = NA
  ) %>%
  as.data.frame()
# Create the Inbred column
VI.wide$Inbred <- lapply(strsplit(as.character(VI.wide$Pedigree), '/'), '[', 1) %>% unlist()
VI.wide <- VI.wide %>%
  relocate(Pedigree.Env:Pedigree, Inbred)

# To create the cross files, we need to subset the genomic data file for each cross. We will subset in a loop,
# stopping the loop if the number of inbreds in the subsetted genotypic data file does not equal the number
# of inbreds obtained from the specific Env the loop is on
for (env in envs) {
  df.loop <- VI.wide %>% filter(Env == env)
  genomic.loop <- supp.genomic %>% filter(lines %in% df.loop$Inbred)
  
  n.inbred  <- length(df.loop$Inbred) # How many rows we expect
  n.genomic <- length(unique(genomic.loop$lines)) # How many unique lines we obtained
  
  if (n.genomic != n.inbred) {
    stop(
      sprintf("Mismatch for Env '%s': %d genomic lines vs %d inbred rows",
              et, n.genomic, n.inbred)
    )
  }
  
  # Clean up the df.loop file to only contain a column called "Inbred" followed by the FPCs
  df.loop.clean <- df.loop %>% select(Inbred, FPC1_RCC:FPC6_RCC)
  
  fwrite(genomic.loop, paste0(data_path, env, '.subset.genotype.data.csv'))
  fwrite(df.loop.clean, paste0(data_path, env, '.subset.phenotype.AGDD.FPC.data.csv'))
}

# Now we will read in the supplementary data file JSON and change the "pheno" and "file" fields
template <- read_json(paste0(data_path, 'cm.mbp.ss.100k.sites.json'), simplifyVector = F)
geno_and_pheno_files <- list.files(path = data_path, pattern = '\\.csv$')
geno_files <- geno_and_pheno_files[grepl('genotype', geno_and_pheno_files)]
pheno_files <- geno_and_pheno_files[grepl('phenotype', geno_and_pheno_files)]
geno_pheno_file_df <- data.frame(Geno = geno_files, Pheno = pheno_files)
geno_pheno_file_df$Env.Geno <- lapply(strsplit(as.character(geno_pheno_file_df$Geno), '\\.'), function(x) paste(x[1], x[2], sep = '.')) %>% unlist()
geno_pheno_file_df$Env.Pheno <- lapply(strsplit(as.character(geno_pheno_file_df$Pheno), '\\.'), function(x) paste(x[1], x[2], sep = '.')) %>% unlist()
# Test that the data frame was constructed correctly
identical(geno_pheno_file_df$Env.Geno, geno_pheno_file_df$Env.Pheno)

# Edit both the "pheno" and "file" fields
for (env in envs) {
  
  new_json <- template # Make a new copy of the template within each loop iteration
  geno_pheno_file_df_loop <- geno_pheno_file_df %>% filter(Env.Geno == env)
  new_json$geno <- geno_pheno_file_df_loop$Geno
  new_json$pheno <- geno_pheno_file_df_loop$Pheno
  
  out_name <- paste0(env, '.AGDD.cross.setup.json')
  out_path <- file.path(data_path, out_name)
  
  write_json(new_json,
             path       = out_path,
             auto_unbox = TRUE,
             pretty     = TRUE)
  
}

# Running QTL analysis for each environment
for (env in envs) {
  # Import cross data for this loop iteration's specific environment/tester combination
  cross.loop <- read_cross2(paste0(data_path, env, '.AGDD.cross.setup.json'))
  print(paste('Checking cross with check_cross2:', check_cross2(cross.loop)))
  print(summary(cross.loop))
  # Identify potential sample duplicates (any with >95% match are removed)
  compare_geno_result <- compare_geno(cross.loop)
  compare_geno_df <- as.data.frame(summary(compare_geno_result))
  sample_duplicates <- compare_geno_df %>%
    filter(prop_match > 0.95) %>%
    select(ind1, ind2) %>%
    unlist(use.names = FALSE) %>%
    unique()
  genotypes_retain <- !grepl(paste(sample_duplicates, collapse = '|'), ind_ids(cross.loop))
  # Subset the cross object to exclude these individuals
  cross.loop.subset1 <- cross.loop[genotypes_retain,]
  
  # Extract components from the cross
  pheno.loop <- cross.loop.subset1$pheno
  gmap <- cross.loop.subset1$gmap
  pmap <- cross.loop.subset1$pmap
  # Calculate genotype probabilities (following error probability in Michel et al., 2022)
  message('Calculating genotype probabilities for ', env)
  probs <- calc_genoprob(cross.loop.subset1, error_prob = 0.01, cores = 0)
  plot_genoprob(probs, gmap, ind = 'W10004_0007', chr = 1)
  # Determine number of crossover events per individual
  message('Determining maximum marginal probabilities for ', env)
  maxmarg_result <- maxmarg(probs, cores = 0)
  plot_onegeno(maxmarg_result, gmap, ind = 'W10004_0007')
  # Locate crossover events
  xo_locations <- locate_xo(maxmarg_result, map = pmap, cores = 0)
  # Count crossovers in Subset A (W10004_0001 to W10004_04xx) and Subset B (W10004_0500 to W1004_xxxx)
  # Remove any lines in A with > 150 crossovers and any in B with > 250 crossovers
  # Return the genotype names with ind_ids(cross.loop.subset1)
  xo_counts <- map_int(
    ind_ids(cross.loop.subset1),                # Iterate over every genotype name
    function(g) {                  # g is the current genotype name
      # For this genotype, sum the vector lengths across all chromosomes
      sum(
        map_int(
          xo_locations,            # Iterate over every chromosome list
          function(chr) {          # chr is one chromosome's list
            vec <- chr[[g]]        # Numeric vector for genotype g (may be NULL)
            if (is.null(vec)) 0L else length(vec)
          }
        )
      )
    }
  )
  # Data frame of crossover tallies
  xo_df <- data.frame(Genotype = ind_ids(cross.loop.subset1), TotalXO = xo_counts)
  # Identify Subset A or B and flag those in A with > 150 or B with > 250
  xo_df <- xo_df %>% mutate(
    suffix  = as.integer(str_extract(Genotype, "(?<=_)\\d+")),
    Subset  = if_else(suffix < 500, "A", "B"),
    Flag    = case_when(
      Subset == "A" & TotalXO > 150 ~ TRUE,
      Subset == "B" & TotalXO > 250 ~ TRUE,
      TRUE                          ~ FALSE)
  ) %>%
    select(-suffix)
  # Count the Subset A and B individuals that will be removed
  print(table(xo_df$Subset, xo_df$Flag))
  # Save the final vector of genotypes that need to be retained - this will be used for subsetting and 
  # sorting everything downstream!
  genotypes_retain_final <- xo_df[xo_df$Flag == 'FALSE',]$Genotype
  
  # Subset the phenotypic data and genotype probabilities
  pheno.loop.final <- pheno.loop[genotypes_retain_final,]
  probs.final <- probs[genotypes_retain_final,]
  # Ensure identical ordering
  identical(get_common_ids(pheno.loop.final, probs.final), genotypes_retain_final)
  
  # Calculate kinship matrix
  message('Calculating kinship matrix for ', env)
  kinship <- calc_kinship(probs = probs.final, type = 'loco', cores = 0)
  # Perform a genome scan
  message('Performing genome scan for ', env)
  scan_result <- scan1(genoprobs = probs.final, pheno = pheno.loop.final, kinship = kinship, cores = 0)
  plot_scan1(x = scan_result, map = pmap, lodcolumn = 'FPC1_RCC')
  # Establish statistical significance with a permutation test
  message('Performing permutation test for ', env)
  perm_result <- scan1perm(genoprobs = probs.final, pheno = pheno.loop.final, n_perm = 1000, cores = 0)
  summary(perm_result)
  # Find LOD peaks (calculate Bayesian credible intervals)
  LOD_peaks <- find_peaks(scan1_output = scan_result,
                          map = pmap,
                          threshold = summary(perm_result),
                          drop = 5)
  
  # Gather all of the useful outputs and save them into a list that will be saved as a .rds
  QTL_list <- list()
  QTL_list[[env]] <- list('sample.duplicates' = sample_duplicates,
                         'xo.locations' = xo_locations,
                         'xo.dataframe' = xo_df,
                         'max.marg.prob' = maxmarg_result,
                         'genotypes.retain.final' = genotypes_retain_final,
                         'pheno.loop.final' = pheno.loop.final,
                         'probs.final' = probs.final,
                         'gmap' = gmap,
                         'pmap' = pmap,
                         'kinship' = kinship,
                         'scan.result' = scan_result,
                         'permutation.result' = perm_result,
                         'LOD.peaks' = LOD_peaks)
  message('Saving .rds for ', env)
  saveRDS(QTL_list, file = paste0(data_path, env, 'QTL.Outputs.rds'))
  
}
