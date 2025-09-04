library(fdapace)
library(dplyr)
library(nlme)
library(sjstats)
library(lme4)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(data.table)

# Path to data (replace with your own path)
data_path <- 'C:/Users/aaron.desalvio/Downloads/Rust_Data/'

# VI BLUEs
VI <- fread(paste0(data_path, '2021_Rust_VI_BLUEs.csv')) %>% as.data.frame()

# FPCA
fpc_out <- list()
predcurves_out <- list()

# VIs to loop through
veg.indices <- sort(unique(VI$Vegetation.Index))

# Set up vector to determine where the temporal data get truncated (ending at peak of NGRDI)
end.flight <- paste0('VI.BLUE.', c(unique(VI$DAP)[5:15]))

for (vi.name in veg.indices) {
  
  for (end in end.flight) {
    print(paste0('Current FPCA: ', vi.name, ' ending on ', end))

    # Filter tall data frame and sort in order of Experiment, DAP, Pedigree
    df <- VI %>% filter(Vegetation.Index == vi.name) %>% arrange(Experiment, DAP, Pedigree)

    wide <- df %>%
      filter(Vegetation.Index == vi.name) %>%
      pivot_wider(
        id_cols   = c(Pedigree, Experiment),
        names_from  = DAP,
        values_from = 'VI.BLUE',
        names_glue  = "VI.BLUE.{DAP}",
        values_fill = NA
      ) %>%
      as.data.frame()

    # Ensure proper ordering of DAPs for the MakeFPCAInputs function
    df.dap <- sort(unique(df$DAP))
    df.dap <- paste0('VI.BLUE.', df.dap)
    wide <- wide[,c(colnames(wide)[1:2], df.dap)]
    rownames(wide) <- paste(wide$Pedigree, wide$Experiment, sep = '.')
    
    # Filtering based on end point
    wide <- wide %>% select(Pedigree:all_of(end))
    
    # Simulate missing data (replace 20% of values with NA)
    set.seed(42)
    missing <- 0.2 # Missing data percentage
    wide <- wide %>%
      mutate(across(where(is.numeric), \(x) {
        ok <- which(!is.na(x))
        n  <- floor(missing * length(ok))
        if (n > 0)
          x[sample(ok, n)] <- NA
        x
      }))
    
    # Extract DAPs from the wide data frame
    DAPs <- as.numeric(gsub('VI.BLUE.', '', colnames(wide)[3:ncol(wide)]))

    # FPCA input preparation
    L3 <- MakeFPCAInputs(IDs = rep(rownames(wide), each = length(sort(DAPs))), 
                        tVec = rep(sort(DAPs), length(rownames(wide))), 
                        t(wide[,-c(1:2)]))

    # Perform FPCA
    FPCAsparse <- FPCA(L3$Ly, L3$Lt,list(dataType = 'Sparse', plot = T,
                                        methodMuCovEst = 'smooth',
                                        methodBwCov = 'GCV',
                                        methodBwMu = 'GCV'))

    # Save functional principal components (as many as possible)
    fpca <- as.data.frame(FPCAsparse$xiEst)
    nExtract <- ncol(fpca)
    fpca <- as.data.frame(fpca[,c(1:nExtract)])
    fpca$Pedigree.Experiment <- rownames(wide)
  
    # Functional variation explained
    # FPCAsparse$cumFVE contains the cumulative variation explained
    cumFVE <- FPCAsparse$cumFVE
    nFVE <- min(length(cumFVE), nExtract)
  
    # Store FVE in a vector
    fpca_variation <- numeric(nFVE)
    if(nFVE >= 1) fpca_variation[1] <- cumFVE[1]
    if(nFVE > 1){
      for(i in 2:nFVE){
        fpca_variation[i] <- cumFVE[i] - cumFVE[i-1]
      }
    }
  
    # Save FVE as new columns
    for(i in 1:nFVE){
      fpca[[paste0("FPC", i, "_FVE")]] <- fpca_variation[i]
    }
  
    # Save column names for FPC scores just as FPC#; FVE is saved in separate columns
    colnames(fpca)[1:nExtract] <- paste0(rep('FPC', nExtract), 1:nExtract)

    # Predicted values
    predicted_scores <- predict(FPCAsparse, newLy = L3$Ly, newLt = L3$Lt)
    predscores <- as.data.frame(predicted_scores$scores)
    rownames(predscores) <- rownames(wide)
    predcurves <- as.data.frame(predicted_scores$predCurves)
    rownames(predcurves) <- rownames(wide)
    predgrid <- as.data.frame(predicted_scores$predGrid)
    predgrid.vec <- predgrid$`predicted_scores$predGrid`

    # Example plot of predicted values
    plot(predgrid.vec, predcurves[1,])

    # Mean of predicted scores across all pedigrees
    mean_pred_curves <- colMeans(predcurves)
    # Plot the mean predicted curve
    plot(predgrid.vec, mean_pred_curves, xlab = 'DAP')
    # Save data frame of mean predicted values
    pred.save <- data.frame(DAP = predgrid.vec, Mean.Pred.Vals = mean_pred_curves, Vegetation.Index = vi.name)

    # Save useful variables for splitting the FPC score data later
    fpca$Pedigree <- lapply(strsplit(as.character(fpca$Pedigree.Experiment), '\\.'), '[', 1) %>% unlist()
    fpca$Experiment <- lapply(strsplit(as.character(fpca$Pedigree.Experiment), '\\.'), '[', 2) %>% unlist()
    fpca$Experiment.VI <- paste(fpca$Experiment, vi.name, sep = '.')
    fpca$Vegetation.Index <- vi.name
    fpca$Ending.Flight <- end
    head(fpca)

    # Save useful variables for splitting the predcurves data later
    colnames(predcurves) <- predgrid.vec
    predcurves$Pedigree.Experiment <- rownames(predcurves)
    predcurves$Pedigree <- lapply(strsplit(as.character(predcurves$Pedigree.Experiment), '\\.'), '[', 1) %>% unlist()
    predcurves$Experiment <- lapply(strsplit(as.character(predcurves$Pedigree.Experiment), '\\.'), '[', 2) %>% unlist()
    predcurves$Vegetation.Index <- vi.name
    predcurves$Ending.Flight <- end
    predcurves <- predcurves[, c(52:56, 1:51)]
  
    # Save outputs into lists
    fpc_out[[paste0(vi.name, '.ending.on.', end)]] <- fpca
    predcurves_out[[paste0(vi.name, '.ending.on.', end)]] <- predcurves
  
  }
}

fpc_out_df <- rbindlist(fpc_out, fill = T)
fpc_out_df <- fpc_out_df %>% relocate(any_of(contains('_FVE'))) %>%
  relocate(Pedigree.Experiment:Ending.Flight)
#fwrite(fpc_out_df, paste0(data_path, '2021_Rust_FPC_Scores_VI_BLUEs_MissingData.csv'))

# Tall format for exporting
stack_out <- list()
for (vi.name in veg.indices) {
  for (end in end.flight) {
    temp <- predcurves_out[[paste0(vi.name, '.ending.on.', end)]]
    long <- temp %>% pivot_longer(
      cols = colnames(temp)[6:ncol(temp)],
      names_to = 'DAP',
      values_to = 'Predicted.Value')
    stack_out[[paste0(vi.name, '.ending.on.', end)]] <- long
  }
}

stack_out_df <- rbindlist(stack_out, fill = TRUE)
#fwrite(stack_out_df, paste0(data_path, '2021_Rust_FPCA_Predicted_Vals_MissingData.csv'))

#save.image(paste0(data_path, 'RData/20250814_FPCA_VI_BLUEs_MissingData_V1.RData'))