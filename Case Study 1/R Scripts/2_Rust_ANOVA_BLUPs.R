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
df <- fread(paste0(data_path, '2021_Rust_Scores_Range_Row_V2.csv')) %>% as.data.frame() %>%
  rename(Experiment = TEST, Replicate = Rep, Plot = PLOT, Pedigree = PEDIGREE)

# Removing outlier plots based on visual inspection during shapefile alignment (e.g., heavy weed pressure / bad germination)
df <- df %>% filter(!Range_Row %in% c('2_3', '15_19'))

# Initialize empty lists to store results
G_list <- list()
VC_list <- list()

# Vector to loop through containing each experiment's name
exp.loop <- unique(df$Experiment)

# Vector to loop through containing rust score column names
rust.scores <- c('Rust.1', 'Rust.2')

# Loop over each experiment
for (exp in exp.loop) {
  
  # Subset data for the current experiment
  df.subset <- df %>% filter(Experiment == exp)
  
  # Convert variables to appropriate data types
  df.subset$Range <- as.factor(df.subset$Range)
  df.subset$Row <- as.factor(df.subset$Row)
  df.subset$Replicate <- as.factor(df.subset$Replicate)
  df.subset$Pedigree <- as.factor(df.subset$Pedigree)
  
  for (rust in rust.scores) {
    
    print(paste0('Current ANOVA model: ', rust, ' for ', exp))
    
    # Fit the linear mixed-effects model after setting "Inf" and NaN as NA
    response <- df.subset[,rust]
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
    VC$Experiment <- exp
    VC$Rust.Score <- rust
    VC$Heritability <- heritability
    
    # Store VC in the list
    VC_list[[paste(exp, rust, sep = '.')]] <- VC
    
    # Extract BLUPs for Pedigree
    G <- coef(fit)$Pedigree
    G$Pedigree <- rownames(G)
    G$Experiment <- exp
    G$Rust.Score <- rust

    # Store G in the list
    G_list[[paste(exp, rust, sep = '.')]] <- G
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

colnames(G_df)[1] <- 'Rust.BLUP'

#fwrite(VC_df, paste0(data_path, '2021_Rust_VarComp_RustScore_BLUPs.csv'))
#fwrite(G_df, paste0(data_path, '2021_Rust_RustScore_BLUPs.csv'))
