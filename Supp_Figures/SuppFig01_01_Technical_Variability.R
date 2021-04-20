
# libraries ---------------------------------------------------------------


library(tidyverse)
library(rio)
library(janitor)


# data import -------------------------------------------------------------

data <- import("data/20180828_Mouse_liver_standard_CV_PeptideGroups.txt")
input_files <- import("data/20180828_Mouse_liver_standard_CV_InputFiles.txt")


# Formatting data ---------------------------------------------------------


# formatting names
# Removing columns that won't be needed
data <- data %>% 
  clean_names() %>%
  select(-contains("abundances_"), 
         -contains("found_in_"),
         -contains("abundance_ratio"))


# Formatting Input file ---------------------------------------------------


# Input file names
input_files <- input_files %>% 
  clean_names()

# converting to NA
input_files[][input_files[] == ""] <- NA

# Subsetting columns and rows
input_files <- input_files %>% 
  filter(!is.na(study_file_id)) %>% 
  select(study_file_id, file_name) %>% 
  mutate(study_file_id = tolower(study_file_id)) %>% 
  rename(condition = "study_file_id")

# Formatting rawfile name
input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
  unlist(strsplit(x, split = "\\", fixed = TRUE))[5]
}))
input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
  unlist(strsplit(x, split = ".", fixed = TRUE))[1]
}))

input_files <- input_files %>% 
  separate(file_name, c("temp1", "temp2", "temp3", "temp4", "temp5", "TechRep"), remove = FALSE, convert = TRUE) %>% 
  select(-contains("temp")) 



# Subetting Protein data --------------------------------------------------


# Cleaning up data
data_tidy <- data %>% 
  filter(contaminant == FALSE) %>% 
  select(sequence, modifications, number_of_missed_cleavages, theo_m_hplus_in_da,
         master_protein_accessions, master_protein_descriptions, 
         contains("abundance_")) %>% 
  gather(temp, abundance, contains("abundance_")) %>% 
  separate(temp, c("temp1", "condition", "isotope", "temp4", "temp5"), sep = "_") %>% 
  select(-contains("temp")) %>% 
  full_join(input_files, .) %>% 
  mutate(abundance = log2(abundance))

data <- data %>% 
  select(-contains("abundance_")) %>% 
  filter(contaminant == FALSE)


# Normalization -----------------------------------------------------------



group_median <- tapply(data_tidy$abundance, data_tidy$TechRep, median, na.rm = TRUE)
global_median <- median(data_tidy$abundance, na.rm = TRUE)

data_tidy <- data_tidy %>% 
  mutate(abundance_norm = abundance - group_median[.$TechRep] + global_median)




# Peptide Summary ---------------------------------------------------------


data_cv <- data_tidy %>% 
  filter(abundance > 0) %>% 
  group_by(sequence, modifications, number_of_missed_cleavages, theo_m_hplus_in_da,
           master_protein_accessions, master_protein_descriptions,
           isotope) %>% 
  summarise(mean = mean(2^abundance_norm, na.rm = TRUE),
            sd = sd(2^abundance_norm, na.rm = TRUE),
            count = n()) %>% 
  ungroup() %>% 
  mutate(cv = sd / mean * 100)

median(data_cv$cv[data_cv$isotope == "light"], na.rm = TRUE)
median(data_cv$cv[data_cv$isotope == "heavy"], na.rm = TRUE)


# Plots -------------------------------------------------------------------


ggplot(data_cv %>% filter(count == 6)) +
  geom_density(aes(x = cv, color = isotope)) +
  theme_classic(base_size = 14) +
  scale_color_brewer(palette = "Set1") +
  annotate("text", x = 150, y = 0.04,
           label = paste0("light = ", round(median(data_cv$cv[data_cv$isotope == "light"], na.rm = TRUE), digits = 2), "%")) +
  annotate("text", x = 150, y = 0.03,
           label = paste0("heavy = ", round(median(data_cv$cv[data_cv$isotope == "heavy"], na.rm = TRUE), digits = 2), "%")) +
  labs(x = "Coefficient of Variation",
       y = "Density")



