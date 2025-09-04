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

# Separate FPC data frames for each environment
TXH2 <- VI.wide %>% filter(Env == 'TXH2.2021') %>% arrange(Pedigree)
TXH3 <- VI.wide %>% filter(Env == 'TXH3.2021') %>% arrange(Pedigree)
identical(TXH2$Pedigree, TXH3$Pedigree)

# Create the "delta" phenotype
delta <- TXH2[,5:10] - TXH3[,5:10]
colnames(delta) <- paste0('Delta.', colnames(delta))
delta <- cbind(TXH2[,c(1,3,4)], delta)

# To create the cross file, we need to subset the genomic data
genomic.subset <- supp.genomic %>% filter(lines %in% delta$Inbred)
  
# Clean up the delta file to only contain a column called "Inbred" followed by the FPCs
delta.clean <- delta %>% select(Inbred, Delta.FPC1_RCC:Delta.FPC6_RCC)
  
fwrite(genomic.subset, paste0(data_path, 'Delta.subset.genotype.data.csv'))
fwrite(delta.clean, paste0(data_path, 'Delta.subset.phenotype.AGDD.FPC.data.csv'))

# Now we will read in the supplementary data file JSON and change the "pheno" and "file" fields
template <- read_json(paste0(data_path, 'cm.mbp.ss.100k.sites.json'), simplifyVector = F)

# Edit both the "pheno" and "file" fields
new_json <- template
new_json$geno <- 'Delta.subset.genotype.data.csv'
new_json$pheno <- 'Delta.subset.phenotype.AGDD.FPC.data.csv'
  
out_name <- 'Delta.AGDD.cross.setup.json'
out_path <- file.path(data_path, out_name)
  
write_json(new_json,
            path       = out_path,
            auto_unbox = TRUE,
            pretty     = TRUE)

#### Running QTL analysis for delta traits ####
# Import cross data
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
  ind_ids(cross.loop.subset1),                # iterate over every genotype name
  function(g) {                  # g is the current genotype name
    # For this genotype, sum the vector lengths across all chromosomes
    sum(
      map_int(
        xo_locations,            # iterate over every chromosome list
        function(chr) {          # chr is one chromosome's list
          vec <- chr[[g]]        # numeric vector for genotype g (may be NULL)
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
kinship <- calc_kinship(probs = probs.final, type = 'loco', cores = 0)
# Perform a genome scan
scan_result <- scan1(genoprobs = probs.final, pheno = pheno.loop.final, kinship = kinship, cores = 0)
plot_scan1(x = scan_result, map = pmap, lodcolumn = 'Delta.FPC1_RCC')
# Establish statistical significance with a permutation test
perm_result <- scan1perm(genoprobs = probs.final, pheno = pheno.loop.final, n_perm = 1000, cores = 0)
summary(perm_result)
# Find LOD peaks (calculate Bayesian credible intervals)
LOD_peaks <- find_peaks(scan1_output = scan_result,
                        map = pmap,
                        threshold = summary(perm_result),
                        drop = 5)
  
# Gather all of the useful outputs and save them into a list that will be saved as a .rds
QTL_list <- list()
QTL_list[['Delta']] <- list('sample.duplicates' = sample_duplicates,
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
saveRDS(QTL_list, file = paste0(data_path, 'Delta.RCC.QTL.Outputs.rds'))

