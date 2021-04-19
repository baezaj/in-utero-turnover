
# Document Setup ----------------------------------------------------------


library(tidyverse) 
library(rio) # data import
library(janitor) # for cleaning data

# Data Import -------------------------------------------------------------


data <- import("data/20190205_fetal_Dev_8_16_24_48_72-(2)_PeptideGroups.txt", setclass = "tibble")
input_files <- import("data/20190205_fetal_Dev_8_16_24_48_72-(2)_InputFiles.txt", setclass = "tibble")


# Formatting data ---------------------------------------------------------


# formatting names and selecting columns for use
data <- data %>% 
  clean_names() %>% 
  select(-contains("found_in_"), 
         -contains("ratio"),
         -contains("abundances_"))

# Input file setup
input_files <- input_files %>% 
  clean_names() %>% 
  mutate(file_id = tolower(file_id)) %>% 
  select(file_id, file_name)


# Formatting raw file names
input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
  unlist(strsplit(x, split = "\\", fixed = TRUE))[8]
}))

input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
  unlist(strsplit(x, split = ".", fixed = TRUE))[1]
}))

input_files <- input_files %>% 
  separate(file_name, c("temp1", "temp2", "temp3", "temp4"), remove = FALSE, convert = TRUE) %>% 
  select(-contains("temp")) %>% 
  mutate(biorep = 1:nrow(.))

# Tidy data ---------------------------------------------------------------


data_tidy <- data %>% 
  filter(contaminant == FALSE) %>% 
  select(sequence, modifications, number_of_missed_cleavages, theo_m_hplus_in_da,
         master_protein_accessions, master_protein_descriptions, 
         contains("abundance")) %>% 
  gather(temp, value, contains("abundance")) %>% 
  separate(temp, c("temp2", "file_id", "isotope", "temp3", "time", "tissue")) %>% 
  right_join(input_files, .) %>%
  mutate(value = log2(value)) %>% 
  select(-contains("temp")) %>%
  mutate(time = as.numeric(str_replace_all(time, "hr", ""))) %>% 
  unite(id, tissue, time, biorep, remove = FALSE)


# Normalization -----------------------------------------------------------


# log2 median normalization
global_median <- median(data_tidy$value, na.rm = TRUE)
group_median <- tapply(data_tidy$value, data_tidy$id, median, na.rm = TRUE)

# Log2 median normalization
data_tidy <- data_tidy %>%
  mutate(abundance_norm = value - group_median[.$id] + global_median)

# removing dataframe
rm(global_median, group_median);gc()

# Formatting file names
# Filtering for 8, 24 and 48 hour runs
data_tidy <- data_tidy %>% 
  mutate(id = str_replace_all(id, "liver_8", "liver_08"),
         id = str_replace_all(id, "heart_8", "heart_08")) %>% 
  filter(time %in% c(8, 24, 48))
         

# Data Export -------------------------------------------------------------

export(data_tidy, file = "intermediate_data/SuppFig02_Peptide_abundance_tidy.csv")


