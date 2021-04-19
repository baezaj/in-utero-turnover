
# Document Setup ----------------------------------------------------------

library(tidyverse)
library(rio)


# Data Import -------------------------------------------------------------

data <- import("data/20171221_Baeza_DevAtlas_Brain_Liver_timecourse_PeptideGroups.txt", setclass = "tibble")
input_files <- import("data/20171221_Baeza_DevAtlas_Brain_Liver_timecourse_InputFiles.txt", setclass = "tibble")
proteins <- import("data/20171221_Baeza_DevAtlas_Brain_Liver_timecourse_Proteins.txt", setclass = "tibble")


# Formatting data ---------------------------------------------------------


# Formatting proteins df
# removing variables that won't be needed
proteins <- proteins %>% 
  rename_all(tolower) %>% 
  rename_all(~str_replace_all(., "\\W", "_")) %>% 
  rename_all(~str_replace_all(., "_{2,}", "_")) %>% 
  select(-contains("found_in"),
         -contains("abundances_"),
         -contains("abundance_"),
         -unique_sequence_id,
         -marked_as)

# formatting names and selecting columns for use
data <- data %>% 
  rename_all(tolower) %>% 
  rename_all(~str_replace_all(., "\\W", "_")) %>% 
  rename_all(~str_replace_all(., "_{2,}", "_")) %>% 
  select(-contains("found_in_"), 
         -contains("ratio"),
         -contains("abundances_"))

# Cleaning up the protein description
data$master_protein_descriptions <- gsub(" \\[OS=Mus musculus]", "", data$master_protein_descriptions)

# Formatting the mida file
mida_median <- mida_median %>% 
  select(file_id, RIA_median)

# Input file setup
input_files <- input_files %>% 
  rename_all(tolower) %>% 
  rename_all(~str_replace_all(., "\\W", "_")) %>% 
  rename_all(~str_replace_all(., "_{2,}", "_")) %>% 
  # rename_all(~str_replace_all(., "__", "_")) %>% 
  mutate(file_id = tolower(file_id),
         biorep = 1:nrow(.)) %>% 
  select(file_id, file_name, biorep) #%>% 
  # full_join(., mida_median)


# Formatting raw file names
input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
  unlist(strsplit(x, split = "\\", fixed = TRUE))[8]
}))

input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
  unlist(strsplit(x, split = ".", fixed = TRUE))[1]
}))

input_files <- input_files %>% 
  separate(file_name, c("temp1", "temp2", "temp3", "tissue", "fetus_id"), sep = "_", remove = FALSE) %>% 
  separate(fetus_id, c("time", "rep"), sep = "-", convert = TRUE, remove = FALSE) %>% 
  select(-contains("temp"))


# Tidy data ---------------------------------------------------------------


data_tidy <- data %>% 
  filter(contaminant == FALSE) %>% 
  select(sequence, modifications, master_protein_accessions, master_protein_descriptions, theo_mhplus_in_da, 
         contains("abundance")) %>% 
  gather(temp, value, contains("abundance")) %>%
  separate(temp, c("temp2", "file_id", "isotope", "temp3")) %>% 
  right_join(input_files, .) %>% 
  mutate(value = log2(value)) %>% 
  select(-contains("temp"))


data <- data %>% 
  filter(contaminant == FALSE) %>% 
  select(-contains("abundance"))


# Generating model input data ---------------------------------------------


# Calculating fraction Heavy / Total
data_ria <- data_tidy %>% 
  filter(!is.na(value)) %>%
  select(-file_id, -file_name) %>% 
  spread(isotope, value) %>% 
  mutate(total = ifelse(is.na(heavy) & !is.na(light), light, log2(2^heavy + 2^light)),
         fraction = ifelse(is.na(heavy), NA, (2^heavy / 2^total)),
         # synthesized = log2((2^total) * (fraction / RIA_median)),
         has_fraction = ifelse(is.na(heavy), FALSE, TRUE)) %>% 
  mutate(heavy = 2^heavy,
         # synthesized = 2^synthesized,
         light = 2^light,
         total = 2^total)

# Adding a time point of zero
data_temp <- data_ria %>% 
  select(tissue, time, rep, biorep) %>% 
  mutate(time = 0) %>% 
  distinct(tissue, time, rep, .keep_all = TRUE) %>% 
  mutate(biorep = 21:26) %>% 
  unite(fetus_id, time, rep, sep = "-", remove = FALSE)

# Dataframe with zero time point for all peptides
data_zero_tp <- data_ria %>% 
  arrange(time) %>% 
  select(-fetus_id, -biorep, -rep#, 
         # -RIA_median
         ) %>% 
  distinct(tissue, sequence, modifications, master_protein_accessions, theo_mhplus_in_da, .keep_all = TRUE) %>% 
  mutate(time = 0,
         # RIA_median = 1,
         heavy = 0,
         light = NA,
         total = NA,
         fraction = 0,
         # synthesized = 0,
         has_fraction = FALSE) %>% 
  right_join(data_temp, .)

all(names(data_zero_tp) %in% names(data_ria))

# Merging the zero time point with the data
data_lm <- bind_rows(data_ria, data_zero_tp)

# removing tables
rm(data_temp, data_zero_tp);gc()


# Data export -------------------------------------------------------------

export(data_tidy, file = "intermediate_data/02_out_Peptide_tidy_data.csv")
export(data, file = "intermediate_data/02_out_Peptide_metadata.csv")
export(data_lm, file = "intermediate_data/02_out_Peptide_modeling_input_data.csv")
export(proteins, file = "intermediate_data/02_out_Protein_metadata.csv")
