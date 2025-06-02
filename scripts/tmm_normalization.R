# Install edgeR package if not already installed
if (!requireNamespace("edgeR", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("edgeR")
}

# Load necessary libraries
library(edgeR)

# Define the path to the filtered expression data
filtered_data_path <- "data/processed/filtered_combined_exp_all.rds"

# Load the filtered expression data
filtered_data <- readRDS(filtered_data_path)

# Convert to a DGEList object for edgeR
dge <- DGEList(counts = as.matrix(filtered_data))

# Perform TMM normalization
dge <- calcNormFactors(dge)

# Get normalized counts
normalized_counts <- cpm(dge)  # cpm() will give you counts per million

# Define the path to save the normalized counts
normalized_data_path <- "data/processed/normalized_counts_tmm.rds"
dir.create(dirname(normalized_data_path), recursive = TRUE, showWarnings = FALSE)

# Save the normalized counts as an RDS file
saveRDS(normalized_counts, file = normalized_data_path)

# Optionally, print the first few rows of the normalized counts
print(head(normalized_counts))
