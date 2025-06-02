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

saveRDS(result_df1, file = paste0(data_path, "combined_dat.rds"))




rm(list = ls())