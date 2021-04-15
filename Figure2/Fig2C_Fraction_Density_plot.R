
# Setup -------------------------------------------------------------------


library(tidyverse)
library(rio) # data import
library(ggridges) # ridge plot
library(ggthemes) # extra ggplot themes

# data import -------------------------------------------------------------

# Importing tidy data
data_tidy <- import("intermediate_data/01_out_Peptide_abundance_tidy.csv")


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

data_fraction <- data_fraction %>% 
  mutate(id = str_replace_all(id, "_", " - "))

sd(c(
  mean(data_fraction$fraction[data_fraction$tissue == "Heart"], na.rm = TRUE),
  mean(data_fraction$fraction[data_fraction$tissue == "Liver"], na.rm = TRUE),
  mean(data_fraction$fraction[data_fraction$tissue == "Lung"], na.rm = TRUE)
))


# Plots -------------------------------------------------------------------

# Ridge plot with fraction
# Figure 2C
ggplot(data_fraction) +
  geom_density_ridges(aes(y = id, x = fraction, fill = tissue), scale = 2.0) +
  scale_fill_brewer(palette = "Set1") +
  labs(y = NULL,
       x = "% Incorporated") +
  guides(fill = FALSE) +
  theme_pander(base_size = 14)







