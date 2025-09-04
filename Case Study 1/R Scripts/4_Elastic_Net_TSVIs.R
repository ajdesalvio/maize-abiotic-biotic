library(caret)
library(glmnet)
library(tidyverse)
library(ggplot2)
library(data.table)
library(doParallel)
library(purrr)

# Path to data (replace with your own path)
data_path <- 'C:/Users/aaron.desalvio/Downloads/Rust_Data/'

# Import time series VIs
tsvi <- fread(paste0(data_path, '2021_Rust_VI_BLUEs.csv')) %>% as.data.frame()

# Work with the median values
tsvi <- tsvi %>% filter(str_detect(Vegetation.Index, 'median'))
tsvi$Pedigree.Experiment <- paste(tsvi$Pedigree, tsvi$Experiment, sep = '.')

# Retain only necessary columns
tsvi <- tsvi %>% select(Pedigree.Experiment, Pedigree, Experiment, DAP, Vegetation.Index, VI.BLUE)

# Unique experiments
exps <- sort(unique(tsvi$Experiment))

# Import rust BLUEs
rust <- read.csv(paste0(data_path, '2021_Rust_RustScore_BLUEs.csv'))
rust$Pedigree.Experiment <- paste(rust$Pedigree, rust$Experiment, sep = '.')

# Rust score names
rust.scores <- c('Rust.1', 'Rust.2')

# Empty lists to store results
varimp_results <- list()
elastic_net_models <- list()

# Ending flights to use for innermost loop
end.flight <- c('VI.BLUE.59','VI.BLUE.62','VI.BLUE.65','VI.BLUE.69','VI.BLUE.78',
                'VI.BLUE.84','VI.BLUE.94','VI.BLUE.103','VI.BLUE.109','VI.BLUE.116',
                'VI.BLUE.120')

# Unique VI names
vi.names <- sort(unique(tsvi$Vegetation.Index))

# Flight dates
flights <- sort(unique(tsvi$DAP))

# We need to create separate data frames for each FPCA ending type
tsvi.list.grouped <- list()
for (end in end.flight) {
  # Determine the end DAP for this loop iteration
  end_DAP <- as.numeric(strsplit(end, '\\.')[[1]][3])
  flights_keep <- flights[1:which(flights == end_DAP)]
  # Build the data frame of time series VIs for the specified time range
  tsvi.list.grouped_i <- tsvi %>%
    filter(DAP %in% flights_keep) %>%
    pivot_wider(
      id_cols = c(Pedigree.Experiment, Pedigree, Experiment),
      names_from = c(DAP, Vegetation.Index),
      values_from = VI.BLUE,
      names_glue = '{Vegetation.Index}.{DAP}.BLUE',
      values_fill = NA)
  # We'd like to order the columns by VI then DAP
  VI_col_order <- paste(rep(vi.names, each = length(flights_keep)), rep(flights_keep, times = length(vi.names)), 'BLUE', sep = '.')
  tsvi.list.grouped_i <- tsvi.list.grouped_i[,c('Pedigree.Experiment', 'Pedigree', 'Experiment', VI_col_order)]
  tsvi.list.grouped[[end]] <- tsvi.list.grouped_i
}

# Leave-one-experiment-out prediction
elastic_net_models <- list()
varimp_results <- list()
loeo_preds <- list()
loeo_metrics <- list()

for (heldout.exp in exps) {
  
  for (rust.score.num in rust.scores) {
    
    # Filter rust for rust score 1 or 2 and retain only necessary columns
    rust.filt <- rust %>% filter(Rust.Score == rust.score.num) %>%
      select(Pedigree.Experiment, Rust.Score, Rust.BLUE)
    
    # Loop through each FPCA ending flight date
    for (end in end.flight) {
      
      message('Elastic net for heldout experiment ', heldout.exp, ', ', rust.score.num, ', end flight date ', end)
      loop_iter_key <- paste(heldout.exp, rust.score.num, end, sep = '_')
      
      # Join the rust scores with the VI predictors
      df0 <- rust.filt %>% inner_join(tsvi.list.grouped[[end]], by = 'Pedigree.Experiment')
      
      # Split LOEO
      df_train <- df0 %>% filter(Experiment != heldout.exp)
      df_test <- df0 %>% filter(Experiment == heldout.exp)
      
      # Obtain predictor columns (making sure to exclude the rust score column)
      pred_cols <- df_train %>%
        select(where(is.numeric)) %>%
        select(-Rust.BLUE) %>%
        names()
      
      Xtr <- df_train[, pred_cols, drop = FALSE]
      ytr <- df_train$Rust.BLUE
      
      # Drop columns with any non-finite values
      bad_cols <- vapply(Xtr, function(z) any(!is.finite(z)), logical(1))
      if (any(bad_cols)) Xtr <- Xtr[, !bad_cols, drop = FALSE]
      
      # Drop columns with near-zero variance
      nzv <- nearZeroVar(Xtr)
      if (length(nzv)) Xtr <- Xtr[, -nzv, drop = FALSE]
      
      # Row-wise complete cases
      keep_tr <- complete.cases(Xtr) & is.finite(ytr)
      Xtr <- Xtr[keep_tr, , drop = FALSE]
      ytr <- ytr[keep_tr]
      exps_here <- df_train$Experiment[keep_tr]
      
      # Build group folds from experiments included in the training data
      foldsOut <- split(seq_along(exps_here), exps_here)
      indexOut <- foldsOut
      index    <- lapply(foldsOut, function(idx) setdiff(seq_along(exps_here), idx))
      
      # Set up k-fold cross-validation
      ctrl <- trainControl(
        method   = "cv",
        index   = index,
        indexOut = indexOut,
        verboseIter = TRUE,
        allowParallel = TRUE,
        returnResamp = "final",
        savePredictions = "final"
      )
      
      # Start parallel after preparing data; ensure cleanup
      cl <- parallel::makePSOCKcluster(max(1, parallel::detectCores() - 4))
      doParallel::registerDoParallel(cl)
      
      # Fit the model
      set.seed(123)
      elastic_fit <- tryCatch(
        train(
          x = Xtr,
          y = ytr,
          method = "glmnet",
          trControl = ctrl,
          preProcess = c("zv","center","scale"),
          tuneLength = 25
        ),
        error = function(e) { message("Elastic net failed: ", conditionMessage(e)); NULL },
        finally = {
          stopCluster(cl)
          foreach::registerDoSEQ()
        }
      )
      if (is.null(elastic_fit) || all(is.na(elastic_fit$results$RMSE))) next
      
      best <- elastic_fit$bestTune
      final_model <- elastic_fit$finalModel
      elastic_net_models[[loop_iter_key]] <- final_model
      
      # Variable importance
      vip <- caret::varImp(elastic_fit)$importance
      if (all(is.nan(vip$Overall))) {
        cf <- as.matrix(stats::coef(final_model, s = best$lambda))[-1, 1]
        vip <- data.frame(Overall = abs(cf))
        rownames(vip) <- rownames(stats::coef(final_model))[-1]
      }
      vip <- vip[order(-vip$Overall), , drop = FALSE]
      
      # Variable importance data wrangling
      vip$Vegetation.Index <- lapply(strsplit(as.character(rownames(vip)), '\\.'), '[', 1) %>% unlist()
      vip$DAP <- lapply(strsplit(as.character(rownames(vip)), '\\.'), '[', 2) %>% unlist()
      vip$Heldout.Experiment <- heldout.exp
      vip$Rust.Score <- rust.score.num
      vip$End.Flight.Date <- end
      
      # Add results to list
      varimp_results[[loop_iter_key]] <- vip
      
      # Predict on the held-out experiment
      if (nrow(df_test) > 0) {
        # Ensure test has the same predictors and order as training
        # (if any training columns are missing in df_test, this will produce an error)
        Xte <- df_test[, colnames(Xtr), drop = FALSE]
        
        # Predict with the best elastic net
        pred_vec <- as.numeric(predict(elastic_fit, newdata = Xte))
        obs_vec  <- df_test$Rust.BLUE
        
        # Evaluate on rows where both sides are finite
        keep  <- is.finite(pred_vec) & is.finite(obs_vec)
        n_ok  <- sum(keep)
        px <- pred_vec[keep]
        ox <- obs_vec[keep]
        n_ok <- length(px)
        sd_pred <- sd(px)
        sd_obs  <- sd(ox)
        
        pearson  <- if (n_ok >= 3 && sd_pred > 0 && sd_obs > 0) cor(px, ox, method = "pearson")  else NA_real_
        spearman <- if (n_ok >= 3 && sd_pred > 0 && sd_obs > 0) cor(px, ox, method = "spearman") else NA_real_
        rmse    <- if (n_ok >= 1) sqrt(mean((px - ox)^2)) else NA_real_
        mae      <- if (n_ok >= 1) mean(abs(px - ox))      else NA_real_
        
        # Determine how sparse the fitted model is given the best lambda
        nz_coefs <- {
          cf <- as.matrix(stats::coef(elastic_fit$finalModel, s = elastic_fit$bestTune$lambda))[-1, , drop = FALSE]
          sum(abs(cf[,1]) > 0)
        }
        
        # Store row-level predictions aligned to df_test
        loeo_preds[[loop_iter_key]] <- dplyr::mutate(
          df_test,
          Pred      = pred_vec,
          Residual  = Rust.BLUE - pred_vec,
          Heldout   = heldout.exp,
          RustScore = rust.score.num,
          EndFlight = end
        )
        
        # Store summary metrics to report
        loeo_metrics[[loop_iter_key]] <- tibble::tibble(
          Heldout   = heldout.exp,
          RustScore = rust.score.num,
          EndFlight = end,
          n_test    = n_ok,
          Pearson   = pearson,
          Spearman  = spearman,
          RMSE      = rmse,
          MAE       = mae,
          SD_Pred   = sd_pred,
          SD_Obs    = sd_obs,
          NZ_Coefs  = nz_coefs
        )
      }
    }
  }
}

# Compile all results and save as an RDS file
elastic_all_results <- list(Final_Models = elastic_net_models,
                            VarImp_Results = varimp_results,
                            LOEO_Preds = loeo_preds,
                            LOEO_Metrics = loeo_metrics)
#write_rds(elastic_all_results, file = paste0(data_path, 'RData/20250813_Elastic_Net_TSVIs.rds'))

# Save variable importance and LOEO metrics as data frames
vip_df <- rbindlist(varimp_results, fill = TRUE)
loeo_metrics_df <- rbindlist(loeo_metrics, fill = TRUE)
#fwrite(vip_df, paste0(data_path, '2021_Rust_Elastic_Net_TSVIs_VarImp.csv'))
#fwrite(loeo_metrics_df, paste0(data_path, '2021_Rust_Elastic_Net_TSVIs_Metrics.csv'))
