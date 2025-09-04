library(nlme)
library(sjstats)
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(MuMIn)
library(data.table)

# Path to data (replace with your own path)
data_path <- 'C:/Users/aaron.desalvio/Downloads/Rust_Data/'

# Import VIs
df <- fread(paste0(data_path, '2021_Rust_VIs_FieldInfo.csv')) %>% as.data.frame()

# Removing outlier plots based on visual inspection during shapefile alignment (e.g., heavy weed pressure / bad germination)
df <- df %>% filter(!Range_Row %in% c('2_3', '15_19'))

# Create a new column to use for the loop by concatenating DAP and Experiment to allow running ANOVA within each Exp/DAP
df$Experiment.DAP <- paste(df$Experiment, df$DAP, sep = '.')

# Initialize empty lists to store results
G_list <- list()
VC_list <- list()

# Vector to loop through containing Experiment.DAPs
exp.dap.loop <- unique(df$Experiment.DAP)

# Vector for inner loop containing VI names
VI.Names <- df %>%
  select(BCC_mean:VEG_median) %>%
  names()

# Loop over each experiment
for (exp.dap in exp.dap.loop) {
  
  # Subset data for the current experiment
  df.subset <- df[df$Experiment.DAP == exp.dap, ]
  
  # Convert variables to appropriate data types
  df.subset$Range <- as.factor(df.subset$Range)
  df.subset$Row <- as.factor(df.subset$Row)
  df.subset$Replicate <- as.factor(df.subset$Replicate)
  df.subset$Pedigree <- as.factor(df.subset$Pedigree)
  
  for (vi.name in VI.Names) {
    
    print(paste0('Current ANOVA model: ', vi.name, ' for ', exp.dap))
    
    # Fit the linear mixed-effects model after setting "Inf" and NaN as NA
    response <- df.subset[,vi.name]
    response[is.na(response) | response == 'Inf' | is.nan(response) | response == '-Inf'] <- NA
    fit <- tryCatch(
      {
        lmer(response ~ (1 | Pedigree) + (1 | Range) + (1 | Row) + (1 | Replicate),
             data = df.subset)
      },
      error = function(e) {
        message("lmer() failed: ", conditionMessage(e), " — skipping this iteration")
        NULL
      }
    )
    # Skip if model didn't fit
    if (is.null(fit) || !inherits(fit, "merMod")) next
    
    # Calculate RMSE
    rmse_value <- sqrt(mean(residuals(fit)^2))
    
    # Calculate R-squared
    R <- MuMIn::r.squaredGLMM(fit)[, "R2c"]  # Conditional R-squared
    
    # Extract variance components of random effects
    VC <- as.data.frame(VarCorr(fit))
    VC$Percent <- round(VC$vcov / sum(VC$vcov) * 100, 2)
    
    # Calculate heritability
    nRep <- length(unique(df.subset$Replicate))
    Vg <- VC[VC$grp == "Pedigree", "vcov"]
    Ve <- VC[VC$grp == "Residual", "vcov"]
    heritability <- Vg / (Vg + (Ve / nRep))
    heritability <- round(heritability, 3)
    
    # Add additional information to VC data frame
    VC$Rmse <- rmse_value
    VC$R_squared <- R
    VC$Experiment.DAP <- exp.dap
    VC$Vegetation.Index <- vi.name
    VC$Heritability <- heritability
    
    # Store VC in the list
    VC_list[[paste(exp.dap, vi.name, sep = '.')]] <- VC
    
    # Extract BLUPs for Pedigree
    G <- coef(fit)$Pedigree
    G$Pedigree <- rownames(G)
    G$Experiment.DAP <- exp.dap
    G$Vegetation.Index <- vi.name

    # Store G in the list
    G_list[[paste(exp.dap, vi.name, sep = '.')]] <- G
  }
}

# Combine the lists into data frames
G_df <- rbindlist(G_list, fill = TRUE) %>% as.data.frame()
VC_df <- rbindlist(VC_list, fill = TRUE) %>% as.data.frame()

# View the first few rows of the results
head(G_df)
head(VC_df)

# Data wrangling before exporting
VC_df <- VC_df %>% select(!c('var1', 'var2'))
VC_df$DAP <- lapply(strsplit(as.character(VC_df$Experiment.DAP), "\\."), "[", 2) %>% unlist()
VC_df$Experiment <- lapply(strsplit(as.character(VC_df$Experiment.DAP), "\\."), "[", 1) %>% unlist()

G_df$DAP <- lapply(strsplit(as.character(G_df$Experiment.DAP), '\\.'), "[", 2) %>% unlist()
G_df$Experiment <- lapply(strsplit(as.character(G_df$Experiment.DAP), '\\.'), "[", 1) %>% unlist()
colnames(G_df)[1] <- 'VI.BLUP'

#fwrite(VC_df, paste0(data_path, '2021_Rust_VarComp_BLUPs.csv'))
#fwrite(G_df, paste0(data_path, '2021_Rust_VI_BLUPs.csv'))
