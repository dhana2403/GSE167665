library(dplyr)

# Combine individual tissue details into a DataFrame
cerebral_cortex <- data.frame(age = ctx_age, id = ctx_id, sample_no = ctx_samples, name= ctx_name, tissue = 'cerebral_cortex')

lung <- data.frame(age = lng_age, id = lng_id, sample_no = lng_samples, name = lng_name, tissue = 'lung')

liver <- data.frame(age = lv_age, id = lv_id, sample_no = lv_samples, name = lv_name, tissue = 'liver')

skeletal_muscle <- data.frame(age = ms_age, id = ms_id, sample_no = ms_samples, name = ms_name, tissue = 'skeletal_muscle')

combined_tissues <- rbind(cerebral_cortex,lung,liver,skeletal_muscle)

rownames(combined_tissues) <- 1:nrow(combined_tissues)

# Change the order of the tissue column to combine it to metadata for batch
combined_tissues$tissue <- factor(combined_tissues$tissue, 
                         levels = c("cerebral_cortex", "liver", "lung", "skeletal_muscle"))

# Sort the DataFrame based on the factor levels
combined_tissues <- combined_tissues %>%
arrange(tissue)

rm(list=ls())

path <- "./data/metadata/"

metadata <- readRDS(paste0(path, "combined_metadata_dat.rds"))

# Rename columns using colnames()
colnames(metadata)[colnames(metadata) == "sample_no"] <- "mouse_no"
colnames(metadata)[colnames(metadata) == "name"] <- "sample_no"

samples_by_tissues = tapply(as.character(metadata$sample_no), INDEX = metadata$tissue, FUN = function(x) unique(c(x)))

#Gene expression data




# Apply the filtering for each tissue type
dat_filtered_tissues = lapply(samples_by_tissues, function(samps) {
  samps = intersect(samps, samplesx)
  dat_filtered[, samps, drop = FALSE]
})























rm(list = ls())