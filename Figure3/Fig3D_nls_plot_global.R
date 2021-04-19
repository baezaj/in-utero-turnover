
# Setup -------------------------------------------------------------------

library(tidyverse)
library(rio)

# Data Import -------------------------------------------------------------

bfl <- import("intermediate_data/Fig03_nls_model_BestFitLine_data.csv")
model_tidy <- import("intermediate_data/Fig03_nls_model_peptide_fraction_broom_tidy.csv")


# Turnover calculation ----------------------------------------------------


# Calculating turnover rates
model_tidy <- model_tidy %>% 
  filter(p.value < 0.05, 
         std.error < 0.025) %>% 
  mutate(t_half = log(2)/estimate)


# joining model summary data with bfl
data <- inner_join(model_tidy, bfl) %>% 
  unite(id, tissue, sequence, modifications, theo_m_hplus_in_da, sep = "_", remove = FALSE)


# plots -------------------------------------------------------------------


ggplot(data) +
  geom_line(aes(x = time, y = pred, color = tissue, group = id), alpha = 0.02) +
  scale_color_brewer(palette = "Set1") +
  facet_wrap(~tissue) +
  theme_bw(base_size = 14) +
  labs(y = "Fraction",
       x = "Time (hours)")





