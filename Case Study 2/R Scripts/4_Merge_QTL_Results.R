library(tidyverse)
library(data.table)

# Path to data (replace with your own path)
data_path <- 'C:/Users/aaron.desalvio/Downloads/Abiotic_Data/'

# BLUPs
blup.add.int <- read.csv(paste0(data_path, 'QTL.Results.BLUPs.Add.Int.long.csv'))
blup.byenv.delta <- read.csv(paste0(data_path, 'QTL.Results.BLUPs.ByEnv.Delta.long.csv'))
blup <- rbind(blup.add.int, blup.byenv.delta)
fwrite(blup, paste0(data_path, 'QTL.Results.BLUPs.Combined.long.csv'))

# Gene Hits
ghits.add.int <- read.csv(paste0(data_path, 'QTL.Results.Gene.Hits.Add.Int.csv'))
ghits.byenv.delta <- read.csv(paste0(data_path, 'QTL.Results.Gene.Hits.ByEnv.Delta.csv'))
ghits <- rbind(ghits.add.int, ghits.byenv.delta)
fwrite(ghits, paste0(data_path, 'QTL.Results.Gene.Hits.Combined.csv'))

# LOD
lod.add.int <- fread(paste0(data_path, 'QTL.Results.LOD.Add.Int.csv')) %>% as.data.frame()
lod.byenv.delta <- fread(paste0(data_path, 'QTL.Results.LOD.ByEnv.Delta.csv')) %>% as.data.frame()
lod <- rbind(lod.add.int, lod.byenv.delta)
fwrite(lod, paste0(data_path, 'QTL.Results.LOD.Combined.csv'))

# LOD peaks (including the close-to-but-not-significant peak on Chr 6 FPC1 RCC)
lodpeaks.add.int <- read.csv(paste0(data_path, 'QTL.Results.LODPeaks.Add.Int.csv'))
lodpeaks.byenv.delta <- read.csv(paste0(data_path, 'QTL.Results.LODPeaks.ByEnv.Delta.NonSigEx.csv'))
lodpeaks <- rbind(lodpeaks.add.int, lodpeaks.byenv.delta)
fwrite(lodpeaks, paste0(data_path, 'QTL.Results.LODPeaks.Combined.csv'))

# Permutation results
perm.add.int <- read.csv(paste0(data_path, 'QTL.Results.Perm.Add.Int.csv'))
perm.byenv.delta <- read.csv(paste0(data_path, 'QTL.Results.Perm.ByEnv.Delta.csv'))
perm <- rbind(perm.add.int, perm.byenv.delta)
fwrite(perm, paste0(data_path, 'QTL.Results.Perm.Combined.csv'))
