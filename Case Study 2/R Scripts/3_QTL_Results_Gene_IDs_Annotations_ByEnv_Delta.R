library(readr)
library(dplyr)
library(patchwork)

# Path to data (replace with your own path)
data_path <- 'C:/Users/aaron.desalvio/Downloads/Abiotic_Data/'

#Use the B73 v5 reference to obtain gene IDs
# Read the GFF3 file (skip comment lines)
gff <- read_tsv(paste0(data_path, 'Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3'),
                comment = "#", col_names = FALSE)

# Assign column names (standard GFF3)
colnames(gff)[1:9] <- c("seqid", "source", "type", "start", "end",
                        "score", "strand", "phase", "attributes")

# Keep only gene rows
genes <- gff %>% filter(type == "gene")

genes <- genes %>%
  mutate(gene_id = sub(".*ID=gene:([^;]+);.*", "\\1", attributes)) %>%
  select(chr = seqid, start, end, gene_id)

# These lines will import the matching chromosomes from your qtl hits & the reference genome
# NOTE: the "Delta" result was (barely) not significant, but since our analysis primarily concerns plasticity
# between genotypes planted in different environmental conditions, we will analyze it anyway
qtl <- read.csv(paste0(data_path, 'QTL.Results.LODPeaks.ByEnv.Delta.NonSigEx.csv'))  # columns: chr, pos, model, LOD

genes <- genes %>%
  filter(grepl("^chr", chr))
qtl <- qtl %>%
  mutate(chr = paste0("chr", chr))

search_window <- 200000

candidates <- list()

for (i in 1:nrow(qtl)) {
  chr_i <- qtl$chr[i]
  pos_i <- qtl$pos[i]

  hits <- genes %>%
    filter(chr == chr_i,
           start <= (pos_i + search_window),
           end >= (pos_i - search_window)) %>%
    mutate(peak_pos = pos_i,
           LOD = qtl$lod[i],
           Phenotype = qtl$lodcolumn[i],
           QTL.Analysis.Name = qtl$QTL.Analysis.Name[i])
  
  candidates[[i]] <- hits
}

candidate_df <- bind_rows(candidates)

# View and save
head(candidate_df)

gene_hits <- as.data.frame(candidate_df$gene_id)

# Extract the gene ID
gene_hits_ids <- sub("ID=([^;]+);.*", "\\1", gene_hits$`candidate_df$gene_id`)
candidate_df$gene_id <- sub("ID=([^;]+);.*", "\\1", candidate_df$gene_id)

#data.table::fwrite(candidate_df, paste0(data_path, 'QTL.Results.Gene.Hits.ByEnv.Delta.csv'))
