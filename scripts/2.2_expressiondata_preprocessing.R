

#Gene expression data

# Load necessary libraries
library(dplyr)
library(tibble)

# Define a function to read and combine data
combine_gene_expression_data <- function(ma_id) {
  # Define the path to the folder containing the text files
  folder_path <- paste0("./data/raw count/", ma_id)
  
  # List all .txt files in the folder
  file_list <- list.files(path = folder_path, pattern = "*.txt", full.names = TRUE)
  
  # Initialize an empty list to store DataFrames
  df_list <- list()
  
  # Loop through each file and read the data
  for (file in file_list) {
    # Read the gene expression data from the file
    temp_df <- read.table(file, header = TRUE)  # Assuming the first row is the header
    
    # Extract the sample name from the file name
    sample_name <- gsub("Sample_(.*?)_count\\.txt", "Sample_\\1", basename(file))
    
    # Rename the columns using the full sample name
    colnames(temp_df) <- c("Gene", sample_name)  # Assuming the first column is 'Gene'
    
    # Append the DataFrame to the list
    df_list[[sample_name]] <- temp_df
  }
  
  # Combine all DataFrames into one
  result_df <- Reduce(function(x, y) {
    merge(x, y, by = "Gene", all = TRUE)  # Merge on the 'Gene' column
  }, df_list)
  
  # Define the path to save the combined DataFrame
  data_path <- "data/raw count/combined/"
  
  # Save the combined DataFrame as an RDS file
  saveRDS(result_df, file = paste0(data_path, "combined_exp_", ma_id, ".rds"))
  
  # Return the combined DataFrame (optional)
  return(result_df)
}

# Combine data for MA1, MA2, and MA3
result_ma1 <- combine_gene_expression_data("MA1")
result_ma2 <- combine_gene_expression_data("MA2")
result_ma3 <- combine_gene_expression_data("MA3")

# Remove the first five rows from result_ma1
result_ma1 <- result_ma1[-(1:5), ]
result_ma2 <- result_ma2[-(1:5), ]
result_ma3 <- result_ma3[-(1:5), ]

final_result <- Reduce(function(x, y) {
  merge(x, y, by = "Gene", all = TRUE)  # Merge on the 'Gene' column
}, list(result_ma1, result_ma2, result_ma3))

final_result <- column_to_rownames(final_result, var = "Gene")

# Define the path to save the final combined DataFrame
final_data_path <- "data/raw count/combined/"
dir.create(final_data_path, recursive = TRUE, showWarnings = FALSE)

# Save the final combined DataFrame as an RDS file
saveRDS(final_result, file = paste0(final_data_path, "combined_exp_all.rds"))

# Define genes of interest
genes_of_interest <- c('ENSMUSG00000022982', 'ENSMUSG00000002985', 'ENSMUSG00000022551', 'ENSMUSG00000048756', 'ENSMUSG00000028991', 'ENSMUSG00000025486')

# Filter the final_result for the genes of interest
filtered_final_result <- final_result[rownames(final_result) %in% genes_of_interest, ]

# Define the path to save the filtered DataFrame
filtered_data_path <- "data/raw count/combined/"
dir.create(filtered_data_path, recursive = TRUE, showWarnings = FALSE)

# Save the filtered DataFrame as an RDS file
saveRDS(filtered_final_result, file = paste0(filtered_data_path, "filtered_combined_exp_all.rds"))

# Optionally print the filtered result
print(filtered_final_result)

# Load the filtered combined expression data
filtered_data_path <- "data/processed/filtered_combined_exp_all.rds"
filtered_result <- readRDS(filtered_data_path)

samplesx = intersect(metadata$sample_no, colnames(filtered_combined_exp_all))

# Apply the filtering for each tissue type
dat_filtered_tissues = lapply(samples_by_tissues, function(samps) {
  samps = intersect(samps, samplesx)
  filtered_combined_exp_all[, samps, drop = FALSE]
})

# Create directory if it does not exist
dir.create('data/processed/readcounts/', recursive = TRUE, showWarnings = FALSE)

# Save the filtered data
sapply(names(dat_filtered_tissues), function(nm) {
  saveRDS(dat_filtered_tissues[[nm]], file.path('data/processed/readcounts/', paste0(nm, '.rds')))
})


