
# Setup -------------------------------------------------------------------

library(tidyverse)
library(rio)
library(tidytext) # using with PCA
library(RColorBrewer) # color palettes
library(irlba) # for sparsed matrix PCA


# data import -------------------------------------------------------------

data <- import("intermediate_data/02_out_Peptide_tidy_data.csv")


# formatting --------------------------------------------------------------

# Filtering for heavy peptides
# Selecting only the needed variables
pca_data_irlba <- data %>% 
  filter(isotope == "heavy",
         !is.na(value)) %>% 
  select(tissue, time, rep, sequence, theo_mhplus_in_da, master_protein_accessions, value) %>% 
  unite(id, tissue, time, rep, sep = "_") %>% 
  unite(feature, master_protein_accessions, sequence, theo_mhplus_in_da, sep = "_")

# Setting seed for reproducibility
set.seed(10) 

# Cast Sparse data
sparse_data <- pca_data_irlba %>%
  cast_sparse(row = id, column = feature, value = value)

# Performing the PCA using the sparsed data
pca_irlba <- prcomp_irlba(sparse_data, n = 8, retx = TRUE, center = TRUE, scale. = TRUE)

# extracting data for PCA plotting
pca_plot_irlba <- data.frame(id = sparse_data@Dimnames[[1]],
                             pca_irlba$x) %>%
  separate(id, c("tissue", "time", "biorep"), sep = "_", remove = FALSE, convert = TRUE) %>%
  unite(condition, time, tissue, sep = "_", remove = FALSE) 

# PC components variation
percent_variation <- pca_irlba$sdev^2 / sum(pca_irlba$sdev^2)

# PCA plot
ggplot(pca_plot_irlba, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = factor(time), shape = tissue), size = 5) +
  theme_bw(base_size = 12) +
  scale_color_brewer(palette = "Set1", direction = -1) +
  labs(title = "Principle Components Analysis",
       x = paste0("Principal component 1 (",  round(percent_variation[1], digits = 2) * 100, "%)"),
       y = paste0("Principal component 2 (",  round(percent_variation[2], digits = 2) * 100, "%)")) +
  expand_limits(x = c(-60, 60), y = c(-40, 40)) +
  coord_fixed()
