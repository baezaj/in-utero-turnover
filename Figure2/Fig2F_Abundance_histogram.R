
# Document Setup ----------------------------------------------------------


library(tidyverse) 
library(rio) # data import


# data import -------------------------------------------------------------


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

# changing id names for aesthetics
data_fraction <- data_fraction %>% 
  mutate(id = str_replace_all(id, "_", " - "))

# Summarizing peptide abundances.
# Selecting the max peptide abundance as a proxy for protein abundance
data_prot <- data_fraction %>% 
  group_by(id, tissue, biorep, master_protein_accessions, master_protein_descriptions) %>% 
  summarize(total = max(total, na.rm = TRUE),
            fraction = median(fraction, na.rm = TRUE)) %>% 
  mutate(has_fraction = ifelse(is.na(fraction), FALSE, TRUE)) %>% 
  filter(total != "-Inf")




# Plots -------------------------------------------------------------------

# Histogram of total protein abundance
ggplot(data_prot) +
  geom_histogram(data = data_prot,
                 aes(x = total, fill = "Total"), 
                 color = "black",
                 bins = 25) +
  geom_histogram(data = data_prot %>% filter(has_fraction == TRUE),
                 aes(x = total, fill = "Fraction"), 
                 color = "black",
                 bins = 25) +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~tissue, nrow = 3) +
  theme_bw(base_size = 12) +
  scale_x_continuous(limits = c(15, 35))


