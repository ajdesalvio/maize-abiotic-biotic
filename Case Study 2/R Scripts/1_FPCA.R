library(tidyverse)
library(data.table)
library(fdapace)

# Path to data (replace with your own path)
data_path <- 'C:/Users/aaron.desalvio/Downloads/Abiotic_Data/'

# Import yield BLUEs
yield <- read.csv(paste0(data_path, 'TXH2.TXH3.2021.Yield.BLUEs.csv'))

# Pedigree overlap between environments
tx2 <- yield %>% filter(Env == 'TXH2.2021') %>% pull(Pedigree) %>% unique()
tx3 <- yield %>% filter(Env == 'TXH3.2021') %>% pull(Pedigree) %>% unique()
overlap <- Reduce(intersect, list(tx2, tx3))

# Import vegetation index BLUEs
vi <- fread(paste0(data_path, 'TXH2.TXH3.2021.VI.BLUEs.RCC.csv')) %>% as.data.frame()
# Work with RCC for this example; filter for overlapping Pedigrees
vi.name <- 'RCC'
vi <- vi %>% filter(Vegetation.Index == vi.name, Pedigree %in% overlap)

# Read in EnvRtype weather data to create a DAP / GDD conversion
envrtype <- read.csv(paste0(data_path, 'EnvRtype.Weather.Data.Cleaned.csv'))
dap.gdd <- envrtype %>% select(Env, DAP, GDD, YYYYMMDD)

# Create the accumulated GDD column (AGDD)
dap.gdd <- dap.gdd %>%
  arrange(Env, DAP) %>%
  group_by(Env) %>%
  mutate(AGDD = cumsum(GDD)) %>%
  ungroup()

# Left-join the AGDD data frame to the VI data frame
vi <- vi %>% left_join(dap.gdd, by = c('Env', 'DAP'))

# Visualize the mean temporal trajectories using DAP vs. AGDD
vi %>% group_by(Env, AGDD) %>%
  summarize(MeanVI = mean(VI.BLUE)) %>%
  ggplot() +
  geom_line(aes(x = AGDD, y = MeanVI, group = Env, color = Env), linewidth=2) +
  ylab('Mean NGRDI values per AGDD')

vi %>% group_by(Env, DAP) %>%
  summarize(MeanVI = mean(VI.BLUE)) %>%
  ggplot() +
  geom_line(aes(x = DAP, y = MeanVI, group = Env, color = Env), linewidth=2) +
  ylab('Mean NGRDI values per DAP')

# FPCA
wide <- vi %>%
  pivot_wider(
    id_cols   = c(Pedigree, Env),
    names_from  = AGDD,
    values_from = 'VI.BLUE',
    names_glue  = "VI.BLUE.{AGDD}",
    values_fill = NA
  ) %>%
  as.data.frame()

# Ensure proper ordering of AGDDs for the MakeFPCAInputs function
df.agdd <- sort(unique(vi$AGDD))
df.agdd <- paste0('VI.BLUE.', df.agdd)
wide <- wide[,c(colnames(wide)[1:2], df.agdd)]
rownames(wide) <- paste(wide$Pedigree, wide$Env, sep = '.')

# Extract AGDDs from the wide data frame
AGDDs <- as.numeric(gsub('VI.BLUE.', '', colnames(wide)[3:ncol(wide)]))

# FPCA input preparation
L3 <- MakeFPCAInputs(IDs = rep(rownames(wide), each = length(sort(AGDDs))), 
                     tVec = rep(sort(AGDDs), length(rownames(wide))), 
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
fpca$Pedigree.Env <- rownames(wide)

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

# Save useful variables
fpca$Pedigree <- lapply(strsplit(as.character(fpca$Pedigree.Env), '\\.'), '[', 1) %>% unlist()
fpca$Env <- lapply(strsplit(as.character(fpca$Pedigree.Env), '\\.'), function(x) paste(x[2], x[3], sep = '.')) %>% unlist()
fpca$Env.VI <- paste(fpca$Env, vi.name, sep = '.')
fpca$Vegetation.Index <- vi.name
head(fpca)

# Save useful variables
colnames(predcurves) <- predgrid.vec
predcurves$Pedigree.Env <- rownames(predcurves)
predcurves$Pedigree <- lapply(strsplit(as.character(predcurves$Pedigree.Env), '\\.'), '[', 1) %>% unlist()
predcurves$Env <- lapply(strsplit(as.character(predcurves$Pedigree.Env), '\\.'), function(x) paste(x[2], x[3], sep = '.')) %>% unlist()
predcurves$Vegetation.Index <- vi.name
predcurves <- predcurves[, c(52:55, 1:51)]

# Convert predcurves to long format
long <- predcurves %>% pivot_longer(
  cols = colnames(predcurves)[5:ncol(predcurves)],
  names_to = 'AGDD',
  values_to = 'Predicted.Value')

# Save results
#write_rds(FPCAsparse, paste0(data_path, 'TXH2.TXH3.2021.RCC.FPCA.rds'))
#fwrite(fpca, paste0(data_path, 'TXH2.TXH3.2021.RCC.FPCA.Scores.csv'))
#fwrite(predcurves, paste0(data_path, 'TXH2.TXH3.2021.RCC.FPCA.Predicted.Vals.csv'))
#fwrite(long, paste0(data_path, 'TXH2.TXH3.2021.RCC.FPCA.Predicted.Vals.Long.csv'))
