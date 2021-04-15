
# Libraries ---------------------------------------------------------------


library(tidyverse)
library(rio)


# data import -------------------------------------------------------------


model_tidy <- import("intermediate_data/05_out_model_nls_peptide_fraction_broom_tidy.csv")

# Filtering ---------------------------------------------------------------


# Filtering data
data_filtered <- model_tidy %>% 
  filter(p.value < 0.05,
         std.error < 0.025) %>% 
  select(tissue, sequence, modifications, master_protein_accessions, master_protein_descriptions, 
         theo_mhplus_in_da, estimate) %>% 
  spread(tissue, estimate)

# Annotating by protein groups
data_filtered <- data_filtered %>% 
  mutate(group = NA) %>% 
  mutate(group = ifelse(str_detect(master_protein_descriptions, "(?i)histone (?i)h"), "histone", group),
         group = ifelse(str_detect(master_protein_descriptions, "(?i)ribosomal pr"), "ribosome", group),
         group = ifelse(str_detect(master_protein_descriptions, "(?i)hemoglobin"), "hemoglobin", group))


# Plot --------------------------------------------------------------------


ggplot(data_filtered) +
  geom_point(aes(x = Liver, y = Brain, color = group), alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  scale_color_brewer(palette = "Set1", na.value = "grey80") +
  scale_y_log10() +
  scale_x_log10() +
  theme_light(base_size = 14) +
  annotate("text", y = .0001, x = .1, label = paste0("rho == ", round(cor(data_filtered$Liver, data_filtered$Brain, use = "pairwise.complete.obs", method = "spearman"), digits = 4)), parse = TRUE) +
  labs(color = NULL)

