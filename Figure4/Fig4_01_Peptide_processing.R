
# libraries ---------------------------------------------------------------


library(tidyverse)
library(rio)
library(ggthemes)
library(broom)


# data import -------------------------------------------------------------


data <- import("data/20180811_Mouse_Fetal_Development_Liver-(1)_PeptideGroups.txt")
input_files <- import("data/20180811_Mouse_Fetal_Development_Liver-(1)_InputFiles.txt")
sample_list <- import("data/20180717_Sample_prep.xlsx")
# mida_median <- import("processed_data/20181120_MIDA_analysis_summary.csv")


# Formatting Input file ---------------------------------------------------


# Input file names
input_files <- input_files %>%
  rename_all(tolower) %>% 
  rename_all(~str_replace_all(., "\\W", "_"))


# converting to NA
input_files[][input_files[] == ""] <- NA

# Subsetting columns and rows
input_files <- input_files %>% 
  filter(!is.na(study_file_id)) %>% 
  select(study_file_id, file_name) %>% 
  mutate(study_file_id = tolower(study_file_id)) %>% 
  rename(file_id = "study_file_id")

# Formatting rawfile name
input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
  unlist(strsplit(x, split = "\\", fixed = TRUE))[9]
}))
input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
  unlist(strsplit(x, split = ".", fixed = TRUE))[1]
}))

input_files <- input_files %>% 
  separate(file_name, c("temp1", "temp2", "temp3", "temp4", "fetus_id"), remove = FALSE, convert = TRUE) %>% 
  select(-contains("temp")) 

input_files$fetus_id[input_files$file_name == "20180811_E13_Liver_2_126"] <- 39
input_files$fetus_id[input_files$file_name == "20180811_E13_Liver_4_130"] <- 46
input_files$fetus_id[input_files$file_name == "20180811_E13_Liver_6_138"] <- 54
input_files$fetus_id[input_files$file_name == "20180811_E14_Liver_2_39"] <- 68
input_files$fetus_id[input_files$file_name == "20180811_E14_Liver_4_46"] <- 77
input_files$fetus_id[input_files$file_name == "20180811_E14_Liver_6_54"] <- 86
input_files$fetus_id[input_files$file_name == "20180811_E16_Liver_2_68"] <- 126
input_files$fetus_id[input_files$file_name == "20180811_E16_Liver_4_77"] <- 130
input_files$fetus_id[input_files$file_name == "20180811_E16_Liver_6_86"] <- 138


# Formatting sample list --------------------------------------------------


sample_list <- sample_list %>% 
  select(fetus_id, sample_name, batch)


sample_list$fetus_id[sample_list$sample_name == "Batch_1"] <- 1
sample_list$fetus_id[sample_list$sample_name == "Batch_2"] <- 2
sample_list$fetus_id[sample_list$sample_name == "Batch_3"] <- 3
sample_list$fetus_id[sample_list$sample_name == "Batch_4"] <- 4
sample_list$fetus_id[sample_list$sample_name == "Batch_5"] <- 5
sample_list$fetus_id[sample_list$sample_name == "Batch_6"] <- 6



# Merging with input files
input_files <- right_join(sample_list, input_files) %>% 
  select(-file_name)

input_files$sample_name[input_files$sample_name == "Batch_1"] <- "Batch_Liver_0hr_1"
input_files$sample_name[input_files$sample_name == "Batch_2"] <- "Batch_Liver_0hr_2"
input_files$sample_name[input_files$sample_name == "Batch_3"] <- "Batch_Liver_0hr_3"
input_files$sample_name[input_files$sample_name == "Batch_4"] <- "Batch_Liver_0hr_4"


# Cleaning input file names
input_files <- input_files %>% 
  separate(sample_name, c("dev_stage", "tissue", "time", "temp"), remove = FALSE, convert = TRUE) %>% 
  mutate(time = as.numeric(str_sub(time, start = 1L, end = -3L)),
         sample_name = str_replace_all(sample_name, "hr", "")) %>% 
  select(-contains("temp"))

# removing sample list
rm(sample_list);gc()


# Peptide abundance -------------------------------------------------------

# formatting names and selecting columns for use
data <- data %>% 
  rename_all(tolower) %>% 
  rename_all(~str_replace_all(., "\\W", "_")) %>% 
  rename_all(~str_replace_all(., "___", "_")) %>% 
  rename_all(~str_replace_all(., "__", "_")) %>% 
  select(-contains("found_in_"), 
         -contains("ratio"), 
         -contains("abundances_"))

data$master_protein_descriptions <- gsub(" \\[OS=Mus musculus]", "", data$master_protein_descriptions)

# data wrangling
data_tidy <- data %>% 
  filter(contaminant == FALSE) %>% 
  select(sequence, modifications, theo_mhplus_in_da,
         master_protein_accessions, master_protein_descriptions, 
         contains("abundance_")) %>% 
  gather(temp, value, contains("abundance_")) %>% 
  separate(temp, c("temp1", "file_id", "isotope", "temp4", "temp5", "temp6", "temp7", "temp8"), sep = "_") %>% 
  right_join(input_files, .) %>%
  mutate(value = log2(value)) %>% 
  select(-contains("temp"), -file_id)

# Metadata
data <- data %>% 
  filter(contaminant == "FALSE") %>% 
  select(-contains("abundance_"))

# Removing tables
rm(input_files);gc()


# Generating model input data ---------------------------------------------


# calculating fraction labeled
data_ria <- data_tidy %>% 
  filter(!is.na(value)) %>% 
  spread(isotope, value) %>% 
  mutate(total = ifelse(is.na(heavy) & !is.na(light), light, log2(2^heavy + 2^light)),
         fraction = ifelse(is.na(heavy), NA, (2^heavy / 2^total)),
         has_fraction = ifelse(is.na(heavy), FALSE, TRUE)) %>% 
  mutate(heavy = 2^heavy,
         light = 2^light,
         total = 2^total)


# Adding a time point of zero
data_temp <- data_ria %>% 
  select(dev_stage, fetus_id, tissue, time, batch) %>% 
  mutate(time = 0) %>% 
  distinct(dev_stage, tissue, time, batch, .keep_all = TRUE) %>% 
  mutate(fetus_id = batch) %>%
  unite(sample_name, dev_stage, tissue, time, fetus_id, sep = "_", remove = FALSE) %>% 
  filter(dev_stage != "Batch")

# creating a dataframe for the zero timepoint
data_zero_tp <- data_ria %>% 
  arrange(time) %>% 
  select(-dev_stage, -fetus_id, -batch, -sample_name) %>% 
  distinct(tissue, sequence, modifications, theo_mhplus_in_da, 
           master_protein_accessions, master_protein_descriptions, .keep_all = TRUE) %>% 
  mutate(time = 0,
         heavy = 0,
         light = NA,
         total = NA,
         fraction = 0,
         has_fraction = TRUE) %>% 
  right_join(data_temp, .)

all(names(data_ria) %in% names(data_zero_tp))

# combining the zero time point with the data
data_lm <- bind_rows(data_ria, data_zero_tp) %>% 
  filter(!is.na(fraction))

# removing tables
rm(data_temp, data_zero_tp);gc()

# Data export -------------------------------------------------------------

export(data_tidy, file = "intermediate_data/02_out_Liver_Peptide_tidy_data.csv")
export(data, file = "intermediate_data/02_out_Liver_Peptide_metadata.csv")
export(data_lm, file = "intermediate_data/02_out_Liver_Peptide_modeling_input_data.csv")

