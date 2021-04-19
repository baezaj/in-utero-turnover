
# Setup -------------------------------------------------------------------

library(tidyverse)
library(rio) # data import
library(ggridges) # ridge plot
library(ggthemes) # extra ggplot themes

# data import -------------------------------------------------------------


data_tidy <- import("intermediate_data/SuppFig02_Peptide_abundance_tidy.csv")


# Formatting --------------------------------------------------------------


# Dataframe with total and fraction
data_fraction <- data_tidy %>% 
  filter(!is.na(abundance_norm)) %>% 
  select(id, tissue, biorep, sequence, modifications, number_of_missed_cleavages, theo_m_hplus_in_da,
         master_protein_accessions, master_protein_descriptions, 
         isotope, abundance_norm) %>% 
  spread(isotope, abundance_norm) %>% 
  mutate(total = ifelse(is.na(heavy) & !is.na(light), light, log2(2^heavy + 2^light)),
         fraction = ifelse(is.na(heavy), NA, (2^heavy / 2^total) * 100),
         has_fraction = ifelse(is.na(heavy), FALSE, TRUE))


# Formatting data ---------------------------------------------------------

# separating id into it's components
data_fraction <- data_fraction %>% 
  separate(id, c("tissue", "time", "biorep"), sep = "_", remove = FALSE)



# stacked density plot ----------------------------------------------------


# Ridge plot with fraction
ggplot(data_fraction) +
  geom_density_ridges(aes(y = id, x = fraction, fill = time), scale = 2.0) +
  scale_fill_brewer(palette = "Set1") +
  labs(y = NULL,
       x = "% Incorporated") +
  guides(fill = FALSE) +
  theme_pander(base_size = 14)




