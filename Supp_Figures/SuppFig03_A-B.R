
# setup -------------------------------------------------------------------

library(tidyverse)
library(rio)
library(ggridges)
# library(modelr)
# library(broom)
# library(patchwork)
# library(ggthemes)
# library(ggforce)
# library(ggrepel)


# data import -------------------------------------------------------------

data_nls <- import("../Figure3/intermediate_data/02_out_Peptide_modeling_input_data.csv")


# fractional distribution -------------------------------------------------


dat_w <- data_nls %>% 
  filter(time > 0) %>% 
  unite(id, tissue, time, sep = "-", remove = FALSE)

# Supplemental Figure 3A
# Density plot of fractional abundance
# Entire range
ggplot(dat_w) +
  geom_density_ridges(aes(x = fraction, y = tissue, color = factor(time)), alpha = 0) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 14) +
  labs(color = "Time",
       x = "% Incorporated",
       y = NULL)


# Supplemental Figure 3A (Inset)
# Density plot of fractional abundance
# zoomed in range
ggplot(dat_w) +
  geom_density_ridges(aes(x = fraction, y = tissue, color = factor(time)), alpha = 0) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 14) +
  coord_cartesian(xlim = c(-0.05,0.25)) +
  labs(color = "Time",
       x = "% Incorporated",
       y = NULL)

# Supplemental Figure 3B
# Pairwise comparison 
pairwise.wilcox.test(dat_w$fraction, dat_w$id, p.adjust.method = "BH")



