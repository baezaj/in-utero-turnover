
# libraries ---------------------------------------------------------------

library(tidyverse)
library(rio)
library(umap)


# data --------------------------------------------------------------------

# Importing lung and liver data
data_tidy <- import("intermediate_data/01_out_Liver_Lung_combined_Peptide_tidy_data.csv")


# Single point calibration ------------------------------------------------

# Extracting the reference data
std_data <- data_tidy %>% 
  filter(dev_stage == "Std",
         isotope == "light",
         tissue == "Liver") %>% 
  rename(std = "value") %>% 
  select(-sample_name, -fetus_id, -dev_stage,-isotope, -tissue, -time)

# Merging reference samples with full dataset
data_tidy <- full_join(data_tidy, std_data)

# Single point calibration
data_tidy <- data_tidy %>% 
  mutate(abundance = value - std)


# Normalization -----------------------------------------------------------


# extracting log2 median of each run individually
group_median <- tapply(data_tidy$abundance, data_tidy$sample_name, median, na.rm = TRUE)

# extracting log2 median of entire dataset
global_median <- median(data_tidy$value, na.rm = TRUE)

# performing reference sample normalization
data_tidy <- data_tidy %>% 
  mutate(abundance_norm = abundance - group_median[.$sample_name] + global_median) %>% 
  unite(id, dev_stage, time, fetus_id, sep = "_", remove = FALSE)

rm(global_median, group_median);gc()


# Fig 4B ------------------------------------------------------------------


# Boxplot
ggplot(data_tidy %>% filter(time > 0)) +
  geom_boxplot(aes(x = factor(fetus_id), y = abundance_norm, color = factor(time)), outlier.shape = NA) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 10) +
  scale_y_continuous(limits = c(12,32)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = "Unique fetus",
       y = "Abundance (Log2)",
       color = "Time point") +
  facet_grid(vars(tissue), vars(isotope))








# Fig4C -------------------------------------------------------------------




umap_dat <- data_tidy %>% 
  unite(id, master_protein_accessions, sequence, isotope, theo_m_hplus_in_da, sep = "-") %>% 
  select(id, sample_name, abundance_norm) %>% 
  spread(id, abundance_norm) %>% 
  column_to_rownames("sample_name")

# Removing missing values
umap_dat <- t(na.omit(t(umap_dat)))

# performing UMAP
set.seed(21)
umap_res <- umap(umap_dat)


umap_table <- as.data.frame(umap_res$layout) %>% 
  rownames_to_column("sample_name") %>% 
  rename(UMAP1 = "V1",
         UMAP2 = "V2")

umap_table <- umap_table %>% 
  separate(sample_name, into = c("dev_stage", "tissue", "time", "batch", "fetus_id"), 
           remove = FALSE, convert = TRUE) %>% 
  unite(id, tissue, dev_stage, remove = FALSE) %>% 
  mutate(id = factor(id, levels = c("Liver_E13", "Liver_E14", "Liver_E16", "Liver_P0", 
                                    "Lung_E13", "Lung_E14", "Lung_E16", "Lung_P0", 
                                    "Liver_Std", "Lung_Std")))

ggplot(umap_table) +
  geom_point(aes(x = UMAP1, y = UMAP2, color = factor(id), shape = tissue), size = 4) +
  scale_color_brewer(palette = "Paired") +
  theme_light(base_size = 12) +
  labs(color = NULL,
       shape = NULL) +
  scale_y_continuous(breaks = scales::breaks_pretty()) +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  coord_fixed()

