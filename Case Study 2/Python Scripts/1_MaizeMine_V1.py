from intermine.webservice import Service
import pandas as pd
import os

# Path to data (replace with your own path)
data_path = 'C:/Users/aaron.desalvio/Downloads/Abiotic_Data'

# Import QTL peaks and their associated Gene IDs (e.g., Zm00001eb140300)
qtl = pd.read_csv(os.path.join(data_path, 'QTL.Results.Gene.Hits.Combined.csv'))

# Connect to MaizeMine
svc = Service("https://maizemine-v15.rnet.missouri.edu/maizemine/service")

# Attributes to retrieve
views = [
    "primaryIdentifier",
    "secondaryIdentifier",
    "symbol", "name", "description",
    "organism.name",
    "chromosome.primaryIdentifier",
    "chromosomeLocation.start", "chromosomeLocation.end",
    "goAnnotation.ontologyTerm.identifier"
]

# Prepare ID list and split into manageable chunks
gene_ids = qtl["gene_id"].unique().tolist()
chunks   = [gene_ids[i:i + 200] for i in range(0, len(gene_ids), 200)]

# Run the query chunk-by-chunk
rows = []
for chunk in chunks:
    q = svc.new_query("Gene")
    q.add_view(*views)
    q.add_constraint("primaryIdentifier", "ONE OF", chunk, code="A")
    rows.extend(dict(r.items()) for r in q.rows())

# Assemble into a DataFrame and merge back to your QTL table
annot = pd.DataFrame(rows)
annot.columns = annot.columns.str.replace(r"^Gene\.", "", regex=True)

qtl_enriched = qtl.merge(
    annot,
    left_on="gene_id",
    right_on="primaryIdentifier",
    how="left"
)

qtl_enriched.to_csv(os.path.join(data_path, 'QTL.Results.Gene.Hits.Annot.csv'), index=False)
