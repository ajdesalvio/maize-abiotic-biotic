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
# Arrange by Env, Pedigree
VI.wide <- VI.wide %>%
  relocate(Pedigree.Env:Pedigree, Inbred) %>%
  arrange(Env, Pedigree)

# To create the cross file, we need to subset the genomic data
genomic.subset <- supp.genomic %>% filter(lines %in% VI.wide$Inbred)

# Clean up the combined phenotype file to only contain a column called "Inbred", then Env, followed by the FPCs
VI.wide.clean <- VI.wide %>% select(Inbred, Env, FPC1_RCC:FPC6_RCC)

fwrite(genomic.subset, paste0(data_path, 'Combined.subset.genotype.data.csv'))
fwrite(VI.wide.clean, paste0(data_path, 'Combined.subset.phenotype.AGDD.FPC.data.csv'))

# Now we will read in the supplementary data file JSON and change the "pheno" and "file" fields
template <- read_json(paste0(data_path, 'cm.mbp.ss.100k.sites.json'), simplifyVector = F)

# read_cross2 expects unique genotype identifiers, but by definition we have duplicate IDs since the same genotypes
# were planted in two environments. We will read in a previous cross and then expand the genotype probabilities
# after importing the file.

# Import cross data from previous analysis
cross.loop <- read_cross2(paste0(data_path, 'Delta.AGDD.cross.setup.json'))
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
probs <- calc_genoprob(cross.loop.subset1, error_prob = 0.01, cores = 0)
plot_genoprob(probs, gmap, ind = 'W10004_0007', chr = 1)
# Determine number of crossover events per individual
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
          vec <- chr[[g]]        # Numeric vector for genotype g
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
genotypes_retain_final <- xo_df %>% filter(!Flag) %>% pull(Genotype)

# Final genotype subset
cross.loop.final <- cross.loop.subset1[genotypes_retain_final,]
probs.final <- calc_genoprob(cross.loop.final, error_prob = 0.01, cores = 0)
kinship.final <- calc_kinship(probs.final, type = 'loco', cores = 0)

# Phenotypic data with combined data from both environments
pheno.combined <- VI.wide.clean %>%
  distinct(Inbred, Env, .keep_all = TRUE) %>%
  mutate(
    env_fac = factor(Env, levels = c('TXH2.2021', 'TXH3.2021')),
    id_env = paste(Inbred, Env, sep = '.')
  ) %>%
  relocate(Inbred, Env, env_fac, id_env) %>%
  semi_join(data.frame(Inbred = genotypes_retain_final), by = 'Inbred')

# Phenotype data frame with only FPCs as columns (rownames become Inbred x Env)
pheno.mat <- pheno.combined %>%
  select(starts_with('FPC'))
rownames(pheno.mat) <- pheno.combined$id_env

# Expand arrays to Inbred x Env
idx <- match(pheno.combined$Inbred, ind_ids(cross.loop.final))

# Expand probs and kinship in the pheno.combined order
probs.env <- probs.final
for (chr in names(probs.final)) {
  A <- probs.final[[chr]][idx, , , drop = FALSE]
  dimnames(A)[[1]] <- pheno.combined$id_env
  probs.env[[chr]] <- A
}

K.env <- kinship.final
for (chr in names(kinship.final)) {
  K <- kinship.final[[chr]][idx, idx, drop = FALSE]
  rownames(K) <- colnames(K) <- pheno.combined$id_env
  K.env[[chr]] <- K
}

# Covariates
addcovar <- model.matrix(~ env_fac, pheno.combined)[, -1, drop = FALSE]
rownames(addcovar) <- pheno.combined$id_env
intcovar <- addcovar

# Alignment checks - compare rownames of phenotype data with 1st chromosome within probs and kinship
identical(rownames(pheno.mat), dimnames(probs.env[[1]])[[1]])
identical(rownames(pheno.mat), rownames(K.env[[1]]))

# Run the additive model scan
scan_add <- scan1(genoprobs = probs.env, pheno = pheno.mat, kinship = K.env,
                  addcovar = addcovar, cores = 0)
# Run the interactive model scan
scan_int <- scan1(genoprobs = probs.env, pheno = pheno.mat, kinship = K.env,
                  addcovar = addcovar,
                  intcovar = intcovar, cores = 0)

# Establish statistical significance with a permutation test (using the interaction model)
perm_result <- scan1perm(genoprobs = probs.env, pheno = pheno.mat, kinship = K.env,
                         addcovar = addcovar, intcovar = intcovar,
                         n_perm = 1000, cores = 10)
summary(perm_result)

# Find LOD peaks (calculate Bayesian credible intervals)
LOD_peaks_add <- find_peaks(scan1_output = scan_add,
                        map = pmap,
                        threshold = summary(perm_result),
                        drop = 5)
LOD_peaks_int <- find_peaks(scan1_output = scan_int,
                            map = pmap,
                            threshold = summary(perm_result),
                            drop = 5)

# Gather all of the useful outputs and save them into a list that will be saved as a .rds
QTL_list_add <- list()
QTL_list_add[['Additive']] <- list('sample.duplicates' = sample_duplicates,
                            'xo.locations' = xo_locations,
                            'xo.dataframe' = xo_df,
                            'max.marg.prob' = maxmarg_result,
                            'genotypes.retain.final' = genotypes_retain_final,
                            'pheno.loop.final' = pheno.mat,
                            'probs.final' = probs.env,
                            'gmap' = gmap,
                            'pmap' = pmap,
                            'kinship' = K.env,
                            'scan.result' = scan_add,
                            'permutation.result' = perm_result,
                            'LOD.peaks' = LOD_peaks_add)
saveRDS(QTL_list_add, file = paste0(data_path, 'Additive.RCC.QTL.Outputs.rds'))

QTL_list_int <- list()
QTL_list_int[['Interactive']] <- list('sample.duplicates' = sample_duplicates,
                                'xo.locations' = xo_locations,
                                'xo.dataframe' = xo_df,
                                'max.marg.prob' = maxmarg_result,
                                'genotypes.retain.final' = genotypes_retain_final,
                                'pheno.loop.final' = pheno.mat,
                                'probs.final' = probs.env,
                                'gmap' = gmap,
                                'pmap' = pmap,
                                'kinship' = K.env,
                                'scan.result' = scan_int,
                                'permutation.result' = perm_result,
                                'LOD.peaks' = LOD_peaks_int)
saveRDS(QTL_list_int, file = paste0(data_path, 'Interactive.RCC.QTL.Outputs.rds'))

