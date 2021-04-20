
# libraries ---------------------------------------------------------------


library(tidyverse)
library(rio)
library(broom)
library(UpSetR)


# data import -------------------------------------------------------------

model_tidy <- import("intermediate_data/Fig04_nls_model_peptide_fraction_broom_tidy.csv")


# Filtering poor model fits -----------------------------------------------


model_tidy <- model_tidy %>% 
  filter(master_protein_accessions != "") %>% 
  unite(id, tissue, dev_stage, sep = "-", remove = FALSE) %>% 
  filter(p.value < 0.05, std.error < 0.025) %>% 
  mutate(t_half = log(2)/estimate,
         R_count = str_count(sequence, "R"),
         K_count = str_count(sequence, "K"),
         RP_count = str_count(sequence, "RP"),
         KP_count = str_count(sequence, "KP"),
         missed_cleavage = R_count + K_count - RP_count - KP_count - 1) %>% 
  filter(missed_cleavage < 2)


# Functions ---------------------------------------------------------------


# Running a linear model on turnover rate
nested_stats_test <- function(data){
  lm(log10(estimate) ~ dev_stage,
     data = data)
}

count_gestation <- function(x){
  length(unique(x$dev_stage))
}

count_peptides <- function(x){
  length(unique(x$sequence))
}

count_tissue <- function(x){
  length(unique(x$tissue))
}


# Linear model ------------------------------------------------------------


# Data table for stats
data_nest <- model_tidy %>% 
  unite(peptide_id, sequence, theo_m_hplus_in_da, sep = "_", remove = FALSE) %>% 
  group_by(tissue, peptide_id) %>% 
  mutate(id_n = n()) %>% 
  ungroup() %>% 
  mutate(gestation = dev_stage, 
         dev_stage = str_remove_all(dev_stage, "E"),
         dev_stage = str_replace_all(dev_stage, "P0", "20"),
         dev_stage = as.numeric(dev_stage)) %>% 
  group_by(tissue, master_protein_accessions, master_protein_descriptions) %>% 
  nest() %>% 
  mutate(devstage_n = map_dbl(data, count_gestation),
         peptide_n = map_dbl(data, count_peptides))

# Running the model
data_nest <- data_nest %>% 
  filter(devstage_n >= 2) %>%
  mutate(model_data = map(data, nested_stats_test))

# Extracting model coefficients
model_coef <- data_nest %>% 
  mutate(temp = map(model_data, tidy)) %>% 
  unnest(temp) %>% 
  select(-data, -model_data) %>% 
  filter(term != "(Intercept)") %>% 
  arrange(p.value)

# Extracting the data used for the model
data <- data_nest %>% 
  select(-model_data) %>% 
  unnest(data)


# Fig 4D ------------------------------------------------------------------


# data for hemoglobin
hemoglobin <- data %>% 
  mutate(hemoglobin = str_detect(master_protein_descriptions, "(?i)hemoglobin")) %>% 
  filter(hemoglobin == TRUE)

liv_list <- c("P01942", "P02088")

# Boxplot of hemoglobin proteins
ggplot(hemoglobin %>% filter(id_n >= 4, tissue == "Liver", master_protein_accessions %in% liv_list), 
       aes(x = gestation, y = estimate)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(), width = 0.25) +
  guides(color = FALSE) +
  theme_bw(base_size = 10) +
  facet_wrap(~master_protein_descriptions, nrow = 1, scales = "free_x") +
  scale_y_log10()




# Fig4 E ------------------------------------------------------------------


# Filtering for ribosomal proteins
ribosomal <- data %>% 
  mutate(ribosome = str_detect(master_protein_descriptions, "(?i)ribosom")) %>% 
  filter(ribosome == TRUE)


ribo_prot <- ribosomal %>% 
  filter(peptide_n >= 2, !is.na(estimate)) %>%
  ungroup() %>% 
  group_by(id, master_protein_accessions) %>% 
  summarize(estimate = median(estimate, na.rm = TRUE)) %>% 
  ungroup() %>% 
  spread(id, estimate) %>% 
  column_to_rownames("master_protein_accessions")


ribo_prot[!is.na(ribo_prot)] <- 1
ribo_prot[is.na(ribo_prot)] <- 0

ribo_prot <- ribo_prot %>% 
  rownames_to_column("protein")

# Upset plot
upset(ribo_prot, 
      nsets = 8,
      order.by = c("freq"),
      keep.order = TRUE,
      point.size = 2,
      text.scale = 1.5,
      number.angles = 0,
      nintersects = 5,
      sets = c("Liver-E13", "Liver-E14", "Liver-E16", "Liver-P0",
               "Lung-E13", "Lung-E14", "Lung-E16", "Lung-P0")
)



# Fig 4F (top) ------------------------------------------------------------


ribo_liver <- ribo_prot %>% 
  filter(`Liver-E13` == 1,
         `Liver-E14` == 1,
         `Liver-E16` == 1,
         `Liver-P0` == 1,
         `Lung-E13` == 0,
         `Lung-E14` == 0,
         `Lung-E16` == 0,
         `Lung-P0` == 0)

# Boxplot of ribosomal proteins
# Detected in Liver only
ggplot(ribosomal %>% filter(master_protein_accessions %in% ribo_liver$protein, peptide_n >= 2),
       aes(x = gestation, y = estimate, color = tissue)) +
  geom_boxplot(outlier.color = NA) +
  geom_point(position = position_jitterdodge(), width = 0.25) +
  theme_bw(base_size = 10) +
  facet_wrap(~master_protein_descriptions, nrow = 1, labeller = label_wrap_gen(30)) +
  scale_y_log10()



# Fig 4F (bottom) ---------------------------------------------------------


ribo_lung <- ribo_prot %>% 
  filter(`Liver-E13` == 0,
         `Liver-E14` == 0,
         `Liver-E16` == 0,
         `Liver-P0` == 0,
         `Lung-E13` == 1,
         `Lung-E14` == 1,
         `Lung-E16` == 1,
         `Lung-P0` == 1)


# Boxplot of ribosomal proteins
# Detected in lung only
ggplot(ribosomal %>% filter(master_protein_accessions %in% ribo_lung$protein, peptide_n >= 2),
       aes(x = gestation, y = estimate, color = tissue)) +
  geom_boxplot(outlier.color = NA) +
  geom_point(position = position_jitterdodge(), width = 0.25) +
  theme_bw(base_size = 10) +
  facet_wrap(~master_protein_descriptions, nrow = 1, labeller = label_wrap_gen(30)) +
  scale_y_log10()




