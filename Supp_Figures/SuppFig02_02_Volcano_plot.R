
# Setup -------------------------------------------------------------------

library(tidyverse)
library(rio)
library(purrr)
library(broom)


# Data Import -------------------------------------------------------------


data_tidy <- import("intermediate_data/SuppFig02_Peptide_abundance_tidy.csv", setclass = "tibble")


# Functions ---------------------------------------------------------------

# Function to run an anova across the dataframe
anova_nest_function <- function(data){
  lm(abundance_norm ~ 1 + time, 
     data = data,
     na.action = na.exclude)
}


# Stats -------------------------------------------------------------------


# Preparing the data by nesting
data_nest <- data_tidy %>% 
  select(-modifications) %>% 
  filter(isotope == "heavy",
         biorep < 33,
         !is.na(abundance_norm)) %>% 
  mutate(time = factor(time, levels = c("48", "24", "8"))) %>% 
  group_by(tissue, sequence, theo_m_hplus_in_da, 
           master_protein_accessions, master_protein_descriptions) %>% 
  nest()

# Performing the linear model (Anova)
data_nest <- data_nest %>% 
  mutate(temp = map_dbl(data, nrow)) %>% 
  filter(tissue == "liver" & temp == 9 |
           tissue == "heart" & temp == 12) %>% 
  mutate(data_anova = map(data, anova_nest_function))

# Extracting the model information
data_anova <- data_nest %>% 
  mutate(anova_results = map(data_anova, tidy)) %>% 
  unnest(anova_results) %>% 
  filter(term != "(Intercept)") %>% 
  select(-data, -data_anova) 

# p.value correction using the Benjamini-Hocherg method
data_anova <- data_anova %>% 
  group_by(tissue, term) %>% 
  mutate(adj.p.value = p.adjust(p.value, method = "fdr")) %>% 
  ungroup()



# Plots -------------------------------------------------------------------

# Volcano Plot
ggplot(data_anova) +
  geom_point(aes(x = estimate, y = -log10(adj.p.value)), alpha = 0.1) +
  facet_wrap(~term + tissue, nrow = 4) +
  geom_hline(yintercept = -log10(0.01)) +
  labs(x = "log2fc (time - 48hr)") +
  theme_bw()


# Data Export -------------------------------------------------------------


export(data_anova, file = "intermediate_data/SuppFig02_Anova_Liver_Heart_8hr_vs_24hr_48hr.csv")

