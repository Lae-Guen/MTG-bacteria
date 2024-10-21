###MTG index and feature-table extraction from phyloseq###

# Extract mtg index based on metadata information after loading metadata
# Create a SampleID2 column in metadata and give the same ID to oral and fecal samples from the same person

metadata <- read.csv("./yonsei.metadata.csv", header=TRUE)

# Initialize an empty list to store the results
results_list <- list()
# Loop through samples signature1 to signature507
for (i in 1:416) {
  sample_id <- paste("sig", i, sep = "")
  
  tryCatch({
    # Find the rows in metadata with matching SampleID2
    matching_rows <- metadata[metadata$SampleID2 == sample_id, ]
    
    if (nrow(matching_rows) == 0) {
      # Handle the case where no matching metadata is found
      cat("No matching metadata found for sample ID:", sample_id, "\n")
      next
    }
    
    # Use the first row of matching metadata
    matching_metadata <- matching_rows[1, ]
    
    # Extract the information from matching_metadata
    sex_info <- matching_metadata$Sex
    age_info <- matching_metadata$DOB
    group_info <- matching_metadata$Group
    
    # Subset the phyloseq object for the current sample
    subset_d <- subset_samples(d, SampleID2 == sample_id)
    # Perform the filtering and ASV/taxa extraction
    asv_table <- as.data.frame(otu_table(subset_d))
    filtered_asv_table <- asv_table[rowSums(asv_table >= 1) == 2, ]
    taxa_for_sample <- as.data.frame(tax_table(subset_d))
    
    # Find the corresponding oral and feces columns based on SampleID
    oral_col <- as.character(metadata[metadata$SampleID2 == sample_id & metadata$Type == "Mouth", "SampleID"])
    feces_col <- as.character(metadata[metadata$SampleID2 == sample_id & metadata$Type == "Feces", "SampleID"])
    
    if (length(oral_col) == 0 || length(feces_col) == 0) {
      cat("Oral or fecal sample not found for SampleID:", sample_id, "\n")
      next
    }
    
    oral_sample_id <- ifelse(length(oral_col) > 0, oral_col, NA)
    fecal_sample_id <- ifelse(length(feces_col) > 0, feces_col, NA)
    
    # Compile the data only if there are filtered ASVs
    if (nrow(filtered_asv_table) >= 1) {
      asvs_with_taxa <- cbind("ASV" = rownames(filtered_asv_table), filtered_asv_table[, c(oral_col, feces_col)], taxa_for_sample[rownames(filtered_asv_table), ])
      asvs_with_taxa$SampleID2 <- sample_id
      asvs_with_taxa$Oral_sample <- oral_sample_id
      asvs_with_taxa$Feces_sample <- fecal_sample_id
      asvs_with_taxa$Sex <- sex_info
      asvs_with_taxa$DOB <- age_info
      asvs_with_taxa$Group <- group_info
      
      # Rename the columns
      colnames(asvs_with_taxa)[2] <- "Readcount in mouth"
      colnames(asvs_with_taxa)[3] <- "Readcount in feces"
      
      # Store the result in the list
      results_list[[i]] <- asvs_with_taxa
      
      # Check if the result has 1 or more rows
      if (nrow(asvs_with_taxa) >= 1) {
        # Define the filename for the txt file
        filename <- paste("yonsei_mtg_", i, "_results.txt", sep = "")
        
        # Save the result to a txt file
        write.table(asvs_with_taxa, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
      }
    }
  }, error = function(e) {
    cat("Error with sample ID:", sample_id, ". Error message: ", e$message, "\n")
  })
}

# Combine the resulting files created for each person into one file.
file_list <- list.files(pattern = "*results.txt")
combined_data <- lapply(file_list, function(file) {
  read.csv(file, header = TRUE)
}) %>%
  bind_rows()

# Save as CSV file format
write.csv(combined_data, "yonsei.mtg.csv")