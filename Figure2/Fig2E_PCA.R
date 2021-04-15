
# document setup ----------------------------------------------------------


library(tidyverse)
library(rio) # data import
library(tidytext) # using with PCA
library(RColorBrewer) # color palettes
library(irlba) # for sparsed matrix PCA
library(plotly) # for 3d plotting



# data import -------------------------------------------------------------

data_tidy <- import("intermediate_data/01_out_Peptide_abundance_tidy.csv")


# PCA ---------------------------------------------------------------------

# cleaning up data for PCA
pca_data_irlba <- data_tidy %>%
  filter(!is.na(abundance_norm)) %>%
  select(master_protein_accessions, tissue, isotope, biorep, abundance_norm) %>%
  unite(id, tissue, isotope, biorep, sep = "_")

set.seed(10)
# Cast Sparse data
sparse_data <- pca_data_irlba %>%
  cast_sparse(row = id, column = master_protein_accessions, value = abundance_norm)

# Irlba pca
pca_irlba <- prcomp_irlba(sparse_data, n = 11, retx = TRUE, center = TRUE, scale. = TRUE)

# dataframe for PCA data
pca_plot_irlba <- data.frame(id = sparse_data@Dimnames[[1]],
                             pca_irlba$x) %>%
  separate(id, c("tissue", "isotope", "biorep"), sep = "_", remove = FALSE, convert = TRUE) %>%
  unite(condition, tissue, isotope, sep = "_", remove = FALSE) 

# PC components variation
percent_variation <- pca_irlba$sdev^2 / sum(pca_irlba$sdev^2)

# PCA plot
# Figure 2E
ggplot(pca_plot_irlba, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = condition), size = 5, shape = 21) +
  theme_bw(base_size = 12) +
  scale_fill_brewer(palette = "Paired", direction = -1) +
  labs(x = paste0("Principal component 1 (",  round(percent_variation[1], digits = 2) * 100, "%)"),
       y = paste0("Principal component 2 (",  round(percent_variation[2], digits = 2) * 100, "%)")) +
  expand_limits(x = c(-60, 60), y = c(-40, 40)) +
  coord_fixed()



