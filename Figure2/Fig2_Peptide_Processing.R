
# Document Setup ----------------------------------------------------------


library(tidyverse) 
library(rio) # data import
library(janitor) # data cleaning
# library(ggridges) # ridge plot
# library(GGally) # scatterplot matrix
# library(ggthemes) # extra ggplot themes
# library(tidytext) # using with PCA
# library(RColorBrewer) # color palettes
# library(irlba) # for sparsed matrix PCA
# library(plotly) # for 3d plotting
# library(corrplot) # For correlation matrix
# library(viridisLite) # viridis color palette


# Data Import -------------------------------------------------------------


# data <- import("data/20170526_Fetal_Tissue_16hr_injection_Proteins.txt", setclass = "tibble")
input_files <- import("data/20170526_Fetal_Tissue_16hr_injection_InputFiles.txt", setclass = "tibble")
data <- import("data/20170526_Fetal_Tissue_16hr_injection_PeptideGroups.txt", setclass = "tibble")


# Formatting data ---------------------------------------------------------


# formatting names and selecting columns for use
data <- data %>% 
  clean_names() %>% 
  select(-contains("found_in_"), 
         -contains("ratio"),
         -contains("abundances_"))

# Input file setup
input_files <- input_files %>% 
  clean_names() %>% 
  mutate(file_id = tolower(file_id)) %>% 
  select(file_id, file_name)

# Formatting raw file names
input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
  unlist(strsplit(x, split = "\\", fixed = TRUE))[8]
}))

input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
  unlist(strsplit(x, split = ".", fixed = TRUE))[1]
}))

input_files <- input_files %>% 
  separate(file_name, c("temp1", "temp2", "tissue", "biorep"), remove = FALSE, convert = TRUE) %>% 
  select(-contains("temp"))


# Tidy data ---------------------------------------------------------------

data_tidy <- data %>% 
  filter(contaminant == FALSE) %>% 
  select(sequence, modifications, number_of_missed_cleavages, theo_m_hplus_in_da,
         master_protein_accessions, master_protein_descriptions, 
         contains("abundance")) %>% 
  gather(temp, value, contains("abundance")) %>% 
  separate(temp, c("temp2", "file_id", "isotope", "temp3")) %>% 
  right_join(input_files, .) %>% 
  mutate(value = log2(value)) %>% 
  select(-contains("temp")) %>% 
  unite(id, tissue, biorep, remove = FALSE)


# Normalization -----------------------------------------------------------


# log2 median normalization
global_median <- median(data_tidy$value, na.rm = TRUE)
group_median <- tapply(data_tidy$value, data_tidy$file_name, median, na.rm = TRUE)

# Log2 median normalization
data_tidy <- data_tidy %>%
  mutate(abundance_norm = value - group_median[.$file_name] + global_median)

# removing dataframe
rm(global_median, group_median);gc()


# Data Export -------------------------------------------------------------


export(data_tidy, file = "intermediate_data/01_out_Peptide_abundance_tidy.csv")


# Plotting tables ---------------------------------------------------------

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

data_f_wide <- data_fraction %>% 
  filter(!is.na(fraction)) %>% 
  select(id, sequence, modifications, 
         master_protein_accessions, master_protein_descriptions, fraction) %>% 
  spread(id, fraction)


# wide format of data
data_wide <- data_tidy %>%
  filter(!is.na(abundance_norm)) %>% 
  select(id, sequence, modifications, number_of_missed_cleavages, theo_m_hplus_in_da,
         master_protein_accessions, master_protein_descriptions, isotope, abundance_norm) %>% 
  spread(id, abundance_norm)

data_append <- data_wide %>% 
  slice(1:2)
length(8:ncol(data_append))
data_append[1,8:ncol(data_append)] <- t(rep(10, 12))
data_append[2,8:ncol(data_append)] <- t(rep(38, 12))

data_append[1,1] <- "min"
data_append[2,1] <- "max"

data_wide <- rbind(data_wide, data_append)

rm(data_append);gc()



data_cor <- data_tidy %>% 
  filter(!is.na(abundance_norm)) %>% 
  select(tissue, isotope, biorep, 
         sequence, modifications,
         master_protein_accessions, master_protein_descriptions, abundance_norm) %>% 
  unite(id, tissue, isotope, biorep, sep = "_") %>% 
  spread(id, abundance_norm)

# matrix transformation
mat_abundance <- cor(select_if(data_cor, is.numeric), use = "pairwise.complete.obs", method = "pearson")


# PCA ---------------------------------------------------------------------


# PCA if using light and heavy proteins
pca_data_irlba <- data_tidy %>%
  filter(!is.na(abundance_norm)) %>%
  select(master_protein_accessions, tissue, isotope, biorep, abundance_norm) %>%
  unite(id, tissue, isotope, biorep, sep = "_")

# Cast Sparse data
sparse_data <- pca_data_irlba %>%
  cast_sparse(row = id, column = master_protein_accessions, value = abundance_norm)

# Irlba pca
pca_irlba <- prcomp_irlba(sparse_data, n = 11, retx = TRUE, center = TRUE, scale. = TRUE)

# dataframe for PCA data
pca_plot_irlba <- data.frame(id = sparse_data@Dimnames[[1]],
                             pca_irlba$x) %>%
  separate(id, c("tissue", "isotope", "biorep"), sep = "_", remove = FALSE) %>%
  unite(condition, tissue, isotope, sep = "_", remove = FALSE) %>%
  mutate(biorep = as.numeric(biorep))

# PC components variation
percent_variation <- pca_irlba$sdev^2 / sum(pca_irlba$sdev^2)

# PCA plot
ggplot(pca_plot_irlba, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = condition, shape = factor(tissue)), size = 5) +
  theme_bw(base_size = 12) +
  scale_color_brewer(palette = "Paired", direction = -1) +
  labs(title = "Principle Components Analysis",
       x = paste0("Principal component 1 (",  round(percent_variation[1], digits = 2) * 100, "%)"),
       y = paste0("Principal component 2 (",  round(percent_variation[2], digits = 2) * 100, "%)")) +
  expand_limits(x = c(-60, 60), y = c(-40, 40))

ggsave(filename = "figures/PCA_sparsed_matrix_PC1_PC2.png")
ggsave(filename = "figures/PCA_sparsed_matrix_PC1_PC2.pdf")

# 3D scatter plot
plot_ly(pca_plot_irlba, x = ~PC1, y = ~PC2, z = ~PC3,
        color = ~condition, colors = brewer.pal(6, "Paired")) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'PC1'),
                      yaxis = list(title = 'PC2'),
                      zaxis = list(title = 'PC3')))



# Plots -------------------------------------------------------------------

# Histogram of total protein abundance
ggplot(data_fraction) +
  geom_histogram(data = data_fraction,
                 aes(x = total, fill = "Total"), 
                 color = "black",
                 bins = 25) +
  geom_histogram(data = data_fraction %>% filter(has_fraction == TRUE),
                 aes(x = total, fill = "Fraction"), 
                 color = "black",
                 bins = 25) +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~tissue, nrow = 2) +
  theme_bw(base_size = 12)
ggsave(filename = "figures/histogram_total_vs_fraction.png")
ggsave(filename = "figures/histogram_total_vs_fraction.pdf")

# Ridge plot with fraction
ggplot(data_fraction) +
  geom_density_ridges(aes(y = id, x = fraction, fill = tissue), scale = 2.0) +
  scale_fill_brewer(palette = "Set1") +
  labs(y = NULL,
       x = "% Incorporated") +
  guides(fill = FALSE) +
  theme_pander(base_size = 12)
ggsave(filename = "figures/ridgeplot_fraction_density.png")
ggsave(filename = "figures/ridgeplot_fraction_density.pdf")

median(data_fraction$fraction, na.rm = TRUE)



# ggscatmat(data_wide, columns = 3:ncol(data_wide), color = "isotope") +
#   scale_color_brewer(palette = "Set1") +
#   theme_void(base_size = 5) +
#   labs(y = NULL, x = NULL) +
#   expand_limits(x = c(10,38), y = c(10,40)) +
#   guides(color = FALSE)
# 
# ggsave(filename = "figures/scatterplot_matrix_Void.png", width = 5, height = 5)
# ggsave(filename = "figures/scatterplot_matrix_Void.pdf", width = 5, height = 5)


# Correlation plot matrix
corrplot(mat_abundance, 
         order = "hclust", 
         method = "color", 
         # addrect = 4,
         type = "full",
         tl.pos = "dt", 
         tl.col = "black", 
         tl.srt = 60, 
         tl.cex = 0.65, 
         cl.lim = c(min(mat_abundance),1), 
         # col = viridis(50, option = "D", direction = 1),
         col = colorRampPalette(brewer.pal(9,"RdBu"))(100),
         is.corr = F)

ggplot(data_cor) +
  geom_point(aes(x = Liver_light_2, y = Liver_light_3), shape = 21) +
  expand_limits(x = c(10,40), y = c(10,40)) +
  theme_light(base_size = 12)


# cor(data_cor$Lung_heavy_3, data_cor$Lung_heavy_4, use = "pairwise.complete.obs", method = "pearson")

ggsave(filename = "figures/scatter_plot_liver_light_2_3.pdf", width = 5, height = 5)

ggplot(data_cor) +
  geom_point(aes(x = Lung_heavy_3, y = Lung_heavy_4), shape = 21) +
  expand_limits(x = c(10,40), y = c(10,40)) +
  theme_light(base_size = 12)
ggsave(filename = "figures/scatter_plot_lung_heavy_3_4.pdf", height = 5, width = 5)


# # Scatterplot Matrix using ggpairs
# ggpairs(data = data_wide, 
#         mapping = aes(color = isotope, fill = isotope),
#         columns = 3:ncol(data_wide),
#         upper = list(continuous = wrap("cor", size = 3, alignPercent = 1)),
#         lower = list(continuous = wrap("points", alpha = 0.5, size = 0.1)),
#         diag = list(continuous = wrap("blank", alpha = 0))
#         ) + 
#   theme_light(base_size = 6) 
# 
# ggsave(filename = "figures/scatterplot_matrix_ggpairs.pdf", height = 10, width = 10)
# ggsave(filename = "figures/scatterplot_matrix_ggpairs.png", height = 10, width = 10)
# ggsave(filename = "figures/scatterplot_matrix_ggpairs_Light_lower.pdf", height = 10, width = 10)
# ggsave(filename = "figures/scatterplot_matrix_ggpairs_Light_lower.png", height = 10, width = 10)

