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

# Create a new column to use for the loop by concatenating DAP and Experiment to allow running ANOVA within each Experiment/DAP
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
        lmer(response ~ Pedigree + (1 | Range) + (1 | Row) + (1 | Replicate),
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
    
    # Add additional information to VC data frame
    VC$Rmse <- rmse_value
    VC$R_squared <- R
    VC$Experiment.DAP <- exp.dap
    VC$Vegetation.Index <- vi.name
    
    # Store VC in the list
    VC_list[[paste(exp.dap, vi.name, sep = '.')]] <- VC
    
    # Extract fixed effects for 'Pedigree'
    G <- coef(summary(fit)) %>% as.data.frame()
    
    # The baseline (reference) level for a factor in R is the first element
    # returned by the levels() function for that factor. Use this to retrieve the
    # "missing" Pedigree - the one that was set as the intercept
    G$Pedigree <- gsub('Pedigree', '', rownames(G))
    # Retrieve the Intercept pedigree
    G$Pedigree[1] <- levels(df.subset$Pedigree)[1]
    
    # To retrieve BLUEs, we need to add the reference level (intercept) to all
    # other estimate values that were returned
    G$VI.BLUE <- NA
    G$VI.BLUE[1] <- G$Estimate[1]
    G$VI.BLUE[2:nrow(G)] <- G$VI.BLUE[1] + G$Estimate[2:nrow(G)]
    
    # Useful variables for data wrangling
    G$Experiment.DAP <- exp.dap
    G$Vegetation.Index <- vi.name
    
    # Store G in the list
    G_list[[paste(exp.dap, vi.name, sep = '.')]] <- G
  }
}

# Combine the lists into data frames
G_df <- rbindlist(G_list, fill = TRUE)
VC_df <- rbindlist(VC_list, fill = TRUE)

# View the first few rows of the results
head(G_df)
head(VC_df)

# Data wrangling before exporting
VC_df <- VC_df %>% select(!c('var1', 'var2'))
VC_df$DAP <- lapply(strsplit(as.character(VC_df$Experiment.DAP), "\\."), "[", 2) %>% unlist()
VC_df$Experiment <- lapply(strsplit(as.character(VC_df$Experiment.DAP), "\\."), "[", 1) %>% unlist()

G_df$DAP <- lapply(strsplit(as.character(G_df$Experiment.DAP), '\\.'), "[", 2) %>% unlist()
G_df$Experiment <- lapply(strsplit(as.character(G_df$Experiment.DAP), '\\.'), "[", 1) %>% unlist()

#fwrite(VC_df, paste0(data_path, '2021_Rust_VarComp_BLUEs.csv'))
#fwrite(G_df, paste0(data_path, '2021_Rust_VI_BLUEs.csv'))
