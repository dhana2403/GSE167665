# Load necessary libraries
library(edgeR)

# Define the path to the filtered expression data
filtered_data_path <- "data/processed/readcounts/"

# List of specific tissue types (for example)
tissue_types <- c("cerebral_cortex", "liver", "lung", "skeletal_muscle")  # Update with your actual tissue types

# Loop over each tissue type to load, normalize, and save data
for (tissue in tissue_types) {
  # Construct the path for the specific tissue
  tissue_path <- paste0(filtered_data_path, tissue, ".rds")  # Assuming each tissue has a .rds file
  
  # Load the filtered expression data for the specific tissue
  filtered_data <- readRDS(tissue_path)
  
  # Convert to a DGEList object for edgeR
  dge <- DGEList(counts = as.matrix(filtered_data))
  
  # Perform TMM normalization
  dge <- calcNormFactors(dge)
  
  # Get normalized counts (counts per million)
  normalized_counts <- cpm(dge)
  
  # Define the path to save the normalized counts
  normalized_data_path <- paste0("data/processed/readcounts_tmm/", tissue, ".rds")
  
  # Create the directory to save the normalized counts if it doesn't exist
  dir.create(dirname(normalized_data_path), recursive = TRUE, showWarnings = FALSE)
  
  # Save the normalized counts as an RDS file
  saveRDS(normalized_counts, file = normalized_data_path)
  
  # Optionally, print the first few rows of the normalized counts for the tissue
  print(paste("Normalized counts for", tissue))
  print(head(normalized_counts))
}
