
# Document setup ----------------------------------------------------------


library(tidyverse)
library(rio)
library(janitor)
library(broom)

# Data Import -------------------------------------------------------------

data_pd <- import("data/20171221_Baeza_DevAtlas_Brain_Liver_timecourse_PSMs.txt")
input_files <- import("data/20171221_Baeza_DevAtlas_Brain_Liver_timecourse_InputFiles.txt")
data_mq <- import("data/maxquant/evidence.txt")


# Formatting Input file ---------------------------------------------------


# Input file names
input_files <- input_files %>% 
  clean_names()

# Formatting rawfile name
input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
  unlist(strsplit(x, split = "\\", fixed = TRUE))[8]
}))
input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
  unlist(strsplit(x, split = ".", fixed = TRUE))[1]
}))

# Extracting the tissue and fetus id
input_files <- input_files %>% 
  select(file_id, file_name) %>% 
  separate(file_name, c("temp1", "temp2", "temp3", "tissue", "fetus_id"), sep = "_", remove = FALSE) %>% 
  select(-contains("temp"))


# Formatting PD data ------------------------------------------------------

# Cleaning up the names and merging with the input files table
data_pd <- data_pd %>% 
  clean_names() %>% 
  full_join(input_files, .)


# Formatting MaxQuant data ------------------------------------------------

# Cleaning up the names and merging with the input files table
data_mq <- data_mq %>% 
  clean_names() %>% 
  rename(file_name = "raw_file") %>% 
  full_join(input_files, .) %>% 
  mutate(intensity = as.numeric(intensity))

# Preparing MIDA input ----------------------------------------------------


# labeling peptides
pd_input <- data_pd %>% 
  filter(!is.na(precursor_abundance)) %>% 
  select(file_id, tissue, fetus_id, sequence, charge, modifications,
         master_protein_accessions, master_protein_descriptions, 
         theo_m_hplus_in_da, precursor_abundance) %>% 
  mutate(k_count = str_count(sequence, "K"),
         r_count = str_count(sequence, "R")) %>% 
  filter(k_count == 2 & r_count == 0 | k_count == 0 & r_count == 2)

# Selecting PSM with highest abundance
pd_input <- pd_input %>% 
  group_by(file_id, tissue, fetus_id, sequence, charge, 
           modifications, theo_m_hplus_in_da,
           master_protein_accessions, master_protein_descriptions) %>% 
  summarize(precursor_abundance = max(precursor_abundance, na.rm = TRUE)) %>% 
  ungroup()

# Annotating peptides with modification types
pd_input <- pd_input %>% 
  mutate(c_count = str_count(sequence, "C"),
         carb_count = str_count(modifications, "\\(Carbamidomethyl"),
         m_count = str_count(sequence, "M"),
         ox_count = str_count(modifications, "Oxidation"),
         pyro_glu = str_count(modifications, "pyro-Glu"),
         n_pyro = str_count(modifications, "N-Term\\(Pyro-carbamidomethyl"),
         n_ac = str_count(modifications, "\\(Acetyl"),
         label_count = str_count(modifications, "Label"),
         isotope = if_else(label_count == 0, "LL", 
                           if_else(label_count == 1, "HL", "HH")))

# Removing duplicate entries. Possibly due to modification positional isomers
pd_input <- pd_input %>% 
  distinct(file_id, tissue, fetus_id, sequence, charge, theo_m_hplus_in_da,
           master_protein_accessions, master_protein_descriptions,
           precursor_abundance,
           .keep_all = TRUE) 


# PD MIDA -----------------------------------------------------------------


pd_mida <- pd_input %>% 
  select(-label_count, -modifications, -theo_m_hplus_in_da) %>% 
  spread(isotope, precursor_abundance) %>% 
  mutate(RIA = (2*(HH/HL)) / (1 + (2*(HH/HL))),
         quant = "PD") %>% 
  filter(!is.na(RIA)) %>% 
  separate(fetus_id, c("time", "biorep"), remove = FALSE, convert = TRUE) %>% 
  select(file_id, tissue, fetus_id, time, biorep,
         master_protein_accessions, master_protein_descriptions,
         sequence, charge, RIA, quant) %>% 
  rename(proteins = "master_protein_accessions",
         protein_names = "master_protein_descriptions")


# Preparing MaxQuant MIDA -------------------------------------------------


# Filtering for peptides with 2K or 2R
mq_input <- data_mq %>%
  filter(intensity > 0) %>% 
  select(file_id, tissue, fetus_id, proteins, protein_names, 
         sequence, charge, modifications, intensity) %>% 
  mutate(k_count = str_count(sequence, "K"),
         r_count = str_count(sequence, "R")) %>% 
  filter(k_count == 2 & r_count == 0 | k_count == 0 & r_count == 2)

# Summarizing peptide intensities
mq_input <- mq_input %>%   
  group_by(file_id, tissue, fetus_id, proteins, protein_names, 
           sequence, charge, modifications) %>% 
  summarize(intensity_max = max(intensity, na.rm = TRUE)) %>% 
  ungroup()

# Labeling peptide modifications
mq_input <- mq_input %>% 
  mutate(nterm_ac = str_detect(modifications, "Acetyl"),
         oxidation = str_detect(modifications, "Oxidation"),
         silac = str_detect(modifications, "SILAC"),
         label_type = ifelse(silac == FALSE, "LL", "HL"))

mq_input$label_type[grep("2 Lys8_SILAC", mq_input$modifications)] <- "HH"
mq_input$label_type[grep("2 Arg10_SILAC", mq_input$modifications)] <- "HH"


# MaxQuant MIDA -----------------------------------------------------------


# Equation = (2 * HH/HL) / (1 + (2 * HH/HL))
mq_mida <- mq_input %>% 
  select(file_id, tissue, fetus_id, proteins, protein_names, sequence, charge,
         nterm_ac, oxidation, intensity_max, label_type) %>% 
  spread(label_type, intensity_max) %>% 
  mutate(ratio = HH / HL,
         RIA = (2 * ratio) / (1 + (2 * ratio)),
         quant = "MaxQuant") %>% 
  filter(!is.na(RIA)) %>% 
  separate(fetus_id, c("time", "biorep"), remove = FALSE, convert = TRUE) %>% 
  select(file_id, tissue, fetus_id, time, biorep, proteins, protein_names, sequence, charge, 
         RIA, quant)


# Joining MaxQuant and PD data --------------------------------------------


mida <- bind_rows(mq_mida, pd_mida) %>% 
  mutate(file_id = tolower(file_id))



# MIDA summary ------------------------------------------------------------


# Summarizing the precursor RIA
mida_median <- mida %>%
  group_by(file_id, tissue, time, biorep, fetus_id) %>% 
  summarise(RIA_median = median(RIA, na.rm = TRUE),
            RIA_mad = mad(RIA, na.rm = TRUE),
            count = n()) %>% 
  ungroup()

# Summarizing the precursor RIA
mida_brain <- mida_median %>%
  filter(tissue == "Brain")

# Summarizing the precursor RIA
mida_liver <- mida_median %>%
  filter(tissue == "Liver")


summary(lm(RIA_median ~ time, data = mida_brain))
summary(lm(RIA_median ~ time, data = mida_liver))




# Supplemental Figure 3C --------------------------------------------------


# Supplemental Figure 3C
ggplot(mida, aes(x = fetus_id, y = RIA, color = factor(time))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(shape = 22, width = 0.2) + 
  theme_bw(base_size = 15) +
  facet_wrap(~tissue) +
  scale_color_brewer(palette = "Dark2") +
  expand_limits(y = c(0,1)) +
  guides(color = FALSE) +
  labs(y = "RIA",
       x = "Fetus ID",
       color = "Time")


# Supplemental Figure 3D --------------------------------------------------




# non-linear model
mod_nls_liver <- nls(RIA_median ~ exp(-k*time),
                     data = mida_liver,
                     start = list(
                       k = coef(lm(RIA_median ~ time, data = mida_liver))[[2]]
                     )
)

tidy(mod_nls_liver)
log(2)/.124


# non-linear model
mod_nls_brain <- nls(RIA_median ~ exp(-k*time),
                     data = mida_brain,
                     start = list(
                       k = coef(lm(RIA_median ~ time, data = mida_brain))[[2]]
                     )
)

tidy(mod_nls_brain)
log(2)/.113

ggplot(mida_median, 
       aes(x = as.numeric(time), y = RIA_median, color = tissue)) +
  geom_jitter(shape = 22, size = 2, width = 0.2) +
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c(15, 17, 18)) +
  facet_wrap(~tissue) +
  geom_smooth(method = "nls", 
              formula = y ~ exp(-k*x),
              method.args = list(start= c(k = -0.04)),
              se = FALSE, 
              fullrange = TRUE) +
  theme_bw(base_size = 15) +
  expand_limits(x = c(0,7), y = c(0,1)) +
  guides(color = FALSE) +
  labs(y = "Fraction",
       x = "Time (hours)") 


