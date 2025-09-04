library(tidyverse)
library(data.table)

# Path to data (replace with your own path)
data_path <- 'C:/Users/aaron.desalvio/Downloads/Rust_Data/'

# Field map with range/row/pedigree information
map <- read.csv(paste0(data_path, '2021_Rust_Scores_Range_Row_V2.csv'))

# Raw VIs
VI <- fread(paste0(data_path, '2021_Rust_VIs.csv')) %>% as.data.frame()

# Make columns in the VI file for Flight_Date, Range, Row
VI$Flight_Date <- lapply(strsplit(as.character(VI$file_name), '_'), '[', 1) %>% unlist()
VI$Range <- lapply(strsplit(as.character(VI$file_name), '_'), '[', 2) %>% unlist()
VI$Row <- lapply(strsplit(as.character(VI$file_name), '_'), '[', 3) %>% unlist()
VI$Range_Row <- paste(VI$Range, VI$Row, sep = '_')

# Confirm number of unique plots
length(unique(VI$Range_Row))

# Add the data from the field map
VI <- VI %>% left_join(map, by = 'Range_Row')

# Create DAP column
VI$Planting_Date <- '20210329'
VI$Flight_Date <- as.Date(as.character(VI$Flight_Date), format = '%Y%m%d')
VI$Planting_Date <- as.Date(as.character(VI$Planting_Date), format = '%Y%m%d')
VI$DAP <- as.numeric(difftime(VI$Flight_Date, VI$Planting_Date, units = 'days'))

# Relocate useful columns to the front of the data frame
VI <- VI %>% relocate(Flight_Date:DAP)

# Clean up duplicate columns and tidy up names
VI <- VI %>%
  select(!c(Range.y, Row.y)) %>%
  rename(Range = Range.x, Row = Row.x, Pedigree = PEDIGREE, Plot = PLOT, Replicate = Rep, Experiment = TEST)

# Export the VI data frame
#fwrite(VI, paste0(data_path, '2021_Rust_VIs_FieldInfo.csv'))
