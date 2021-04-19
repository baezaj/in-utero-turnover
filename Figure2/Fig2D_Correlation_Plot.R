
# Setup -------------------------------------------------------------------

library(tidyverse)
library(rio) # data import
library(corrplot) # For correlation matrix
library(RColorBrewer) # color palettes

# data import -------------------------------------------------------------

# Importing tidy data
data_tidy <- import("intermediate_data/Fig02_Peptide_abundance_tidy.csv")


# Formatting --------------------------------------------------------------


# wide format of data
data_wide <- data_tidy %>%
  filter(!is.na(abundance_norm)) %>% 
  select(id, sequence, modifications, number_of_missed_cleavages, theo_m_hplus_in_da,
         master_protein_accessions, master_protein_descriptions, isotope, abundance_norm) %>% 
  spread(id, abundance_norm)


data_append <- data_wide %>% 
  slice(1:2)

data_append[1,3:ncol(data_append)] <- rep(floor(min(data_wide[8:19], na.rm = TRUE)), 17)
data_append[2,3:ncol(data_append)] <- rep(ceiling(max(data_wide[8:19], na.rm = TRUE)), 17)

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



# Plots -------------------------------------------------------------------


# Correlation plot matrix
# Figure 2D
corrplot(mat_abundance, 
         order = "hclust", 
         method = "color", 
         type = "full",
         tl.pos = "dt", 
         tl.col = "black", 
         tl.srt = 60, 
         tl.cex = 0.65, 
         cl.lim = c(min(mat_abundance),1), 
         col = colorRampPalette(brewer.pal(9,"RdBu"))(100),
         is.corr = F)


# Figure 2D Inset (Light)
ggplot(data_cor, aes(x = Liver_light_2, y = Liver_light_3)) +
  geom_point(shape = 22) +
  expand_limits(x = c(10,36), y = c(10,36)) +
  theme_classic(base_size = 14) +
  labs(y = NULL,
       x = NULL) +
  annotate("text", x = 25, y = 12, 
           label = paste("R2 = ", round(cor(data_cor$Liver_light_2, data_cor$Liver_light_3, 
                                            use = "pairwise.complete.obs", method = "pearson"), 
                                        digits = 4)))



# Figure 2D Inset (Heavy)
ggplot(data_cor, aes(x = Heart_heavy_3, y = Heart_heavy_2)) +
  geom_point(shape = 22) +
  expand_limits(x = c(10,36), y = c(10,36)) +
  theme_classic(base_size = 14) +
  labs(y = NULL,
       x = NULL) +
  annotate("text", x = 25, y = 12, 
           label = paste("R2 =", round(cor(data_cor$Heart_heavy_3, data_cor$Heart_heavy_2, 
                                            use = "pairwise.complete.obs", method = "pearson"), 
                                        digits = 4)), sep = " ")




