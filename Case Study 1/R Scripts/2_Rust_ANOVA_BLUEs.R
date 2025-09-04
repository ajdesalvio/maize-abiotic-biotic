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
  VC$Experiment <- exp
  VC$Rust.Score <- rust
  
  # Store VC in the list
  VC_list[[paste(exp, rust, sep = '.')]] <- VC
  
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
  G$Rust.BLUE <- NA
  G$Rust.BLUE[1] <- G$Estimate[1]
  G$Rust.BLUE[2:nrow(G)] <- G$Rust.BLUE[1] + G$Estimate[2:nrow(G)]
  
  # Useful variables for data wrangling
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

#fwrite(VC_df, paste0(data_path, '2021_Rust_VarComp_RustScore_BLUEs.csv'))
#fwrite(G_df, paste0(data_path, '2021_Rust_RustScore_BLUEs.csv'))
