
# Document Setup ----------------------------------------------------------


library(tidyverse) 
library(rio) # data import
library(janitor) # data cleaning


# Data Import -------------------------------------------------------------


input_files <- import("data/20170526_Fetal_Tissue_16hr_injection_InputFiles.txt", setclass = "tibble")
data <- import("data/20170526_Fetal_Tissue_16hr_injection_PeptideGroups.txt", setclass = "tibble")


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
  separate(file_name, c("temp1", "temp2", "tissue", "biorep"), remove = FALSE, convert = TRUE) %>% 
  select(-contains("temp"))


# Tidy data ---------------------------------------------------------------

data_tidy <- data %>% 
  filter(contaminant == FALSE) %>% 
  select(sequence, modifications, number_of_missed_cleavages, theo_m_hplus_in_da,
         master_protein_accessions, master_protein_descriptions, 
         contains("abundance")) %>% 
  gather(temp, value, contains("abundance")) %>% 
  separate(temp, c("temp2", "file_id", "isotope", "temp3")) %>% 
  right_join(input_files, .) %>% 
  mutate(value = log2(value)) %>% 
  select(-contains("temp")) %>% 
  unite(id, tissue, biorep, remove = FALSE)


# Normalization -----------------------------------------------------------


# log2 median normalization
global_median <- median(data_tidy$value, na.rm = TRUE)
group_median <- tapply(data_tidy$value, data_tidy$file_name, median, na.rm = TRUE)

# Log2 median normalization
data_tidy <- data_tidy %>%
  mutate(abundance_norm = value - group_median[.$file_name] + global_median)

# removing dataframe
rm(global_median, group_median);gc()


# Data Export -------------------------------------------------------------


export(data_tidy, file = "intermediate_data/Fig02_Peptide_abundance_tidy.csv")

