###MTG index and feature-table extraction from phyloseq###

# Extract mtg index based on metadata information after loading metadata
# Create a SampleID2 column in metadata and give the same ID to oral and fecal samples from the same person

metadata <- read.csv("./yonsei.metadata.csv", header=TRUE)