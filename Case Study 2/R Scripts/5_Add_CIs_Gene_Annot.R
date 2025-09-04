library(tidyverse)
library(data.table)

# Add confidence intervals to annotated gene hits

# Path to data (replace with your own path)
data_path <- 'C:/Users/aaron.desalvio/Downloads/Abiotic_Data/'

annot <- read.csv(paste0(data_path, 'QTL.Results.Gene.Hits.Annot.csv'))
ci <- read.csv(paste0(data_path, 'QTL.Results.LODPeaks.Combined.csv'))

# Column reorganization before joining
ci2 <- ci %>% rename(peak_pos = pos)
ci2 <- ci2 %>% select(QTL.Analysis.Name, ci_lo, ci_hi, peak_pos)

# left join
annot <- annot %>% left_join(ci2, by = c('QTL.Analysis.Name', 'peak_pos'))

# Reorganize column order
annot <- annot %>% relocate(chr, start, end, gene_id, peak_pos, ci_lo, ci_hi)
write.csv(annot, paste0(data_path, 'QTL.Results.Gene.Hits.Annot.CIs.csv'))
