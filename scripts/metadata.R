# Load necessary libraries
library(dplyr)

# Define the path
path <- "Github/RDS files/"

# Read the RDS files
tissue_ids <- readRDS(paste0(path, "tissue_ids.rds"))
individual_ids <- readRDS(paste0(path, "individual_ids.rds"))
ages <- readRDS(paste0(path, "ages.rds"))

# Combine the data frames (assumes they can be combined by columns)
combined_data <- cbind(tissue_ids, individual_ids, ages)

# Save the combined data to a new RDS file
saveRDS(combined_data, file = paste0(path, "combined_data.rds"))

# Combine into a DataFrame
df1 <- data.frame(age = ctx_age, id = ctx_id, sample_no = ctx_samples, name= ctx_name, tissue = 'cerebral_cortex')

df2 <- data.frame(age = lng_age, id = lng_id, sample_no = lng_samples, name = lng_name, tissue = 'lung')

df3 <- data.frame(age = lv_age, id = lv_id, sample_no = lv_samples, name = lv_name, tissue = 'liver')

df4 <- data.frame(age = ms_age, id = ms_id, sample_no = ms_samples, name = ms_name, tissue = 'skeletal_muscle')

combine <- rbind(df1,df2,df3,df4)

rownames(combine) <- 1:nrow(combine)



# Change the order of the tissue column
combine$tissue <- factor(combine$tissue, 
                         levels = c("cerebral_cortex", "liver", "lung", "skeletal_muscle"))

# Sort the DataFrame based on the factor levels
combine <- combine %>%
  arrange(tissue)

result_df1 <- cbind(combine, batch = metadata$batch)

data_path <- "data/metadata/"

saveRDS(result_df1, file = paste0(data_path, "combined_metadata_dat.rds"))

rm(list = ls())

#expression data

# Load necessary libraries
library(dplyr)
library(readr)  # For read_csv or read_table functions

# Define the path to the folder containing the text files
folder_path <- "./data/raw count/MA1"

# List all .txt files in the folder
file_list <- list.files(path = folder_path, pattern = "*.txt", full.names = TRUE)

# Initialize an empty list to store DataFrames
df_list <- list()

# Loop through each file and read the data
for (file in file_list) {
  # Read the gene expression data from the file
  temp_df <- read.table(file, header = TRUE)  # Adjust as necessary
  
  # Extract the full sample name from the file name
  sample_name <- gsub("Sample_(.*?)-count\\.txt", "Sample_\\1", basename(file))
  
  # Rename the columns using the full sample name
  colnames(temp_df) <- c("Gene", sample_name)  # Assuming the first column is 'Gene'
  
  # Append the DataFrame to the list
  df_list[[sample_name]] <- temp_df
}

# Combine all DataFrames into one
result_df <- Reduce(function(x, y) {
  merge(x, y, by = "Gene", all = TRUE)  # Merge on the 'Gene' column
}, df_list)

# Display the resulting DataFrame
print(result_df)
