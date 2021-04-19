
# libraries ---------------------------------------------------------------


library(tidyverse)
library(rio)
library(broom)
library(patchwork)
library(UpSetR)


# data import -------------------------------------------------------------


model_tidy <- import("intermediate_data/02_out_model_nls_peptide_fraction_broom_tidy.csv")
model_glance <- import("intermediate_data/02_out_model_nls_peptide_fraction_broom_glance.csv")
bfl <- import("intermediate_data/02_out_nls_model_BestFitLine_data.csv")
data_lm <- import("intermediate_data/01_out_Liver_Lung_combined_Peptide_modeling_input_data.csv")


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


temp <- model_tidy %>% 
  select(dev_stage:n_runs)

# removing poor model fits from glance data
model_glance <- inner_join(model_glance, temp) %>% 
  unite(id, tissue, dev_stage, sep = "-", remove = FALSE)

# removing poor model fits from glance data
bfl <- inner_join(bfl, temp)%>% 
  unite(id, tissue, dev_stage, sequence, modifications, theo_m_hplus_in_da, sep = "_", remove = FALSE)


rm(temp);gc()





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
  # filter(id_n >= 3) %>% 
  mutate(gestation = dev_stage, 
         dev_stage = str_remove_all(dev_stage, "E"),
         dev_stage = str_replace_all(dev_stage, "P0", "20"),
         dev_stage = as.numeric(dev_stage)) %>% 
  group_by(tissue, master_protein_accessions, master_protein_descriptions) %>% 
  nest() %>% 
  mutate(devstage_n = map_dbl(data, count_gestation),
         # tissue_n = map_dbl(data, count_tissue),
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



# Supp Fig 5A -------------------------------------------------------------


# Turnover rate plots
ggplot(bfl) +
  geom_line(aes(x = time, y = pred, group = id), alpha = 0.02) +
  facet_grid(vars(tissue), vars(dev_stage)) +
  theme_light(base_size = 14) +
  labs(y = "Fraction",
       x = "Time (hours)")


# Supp Fig 5B -------------------------------------------------------------

liv_dat <- data_lm %>% 
  filter(master_protein_accessions %in% liv_list, tissue == "Liver",
         sequence != "IGGHGAEYGAEALER") %>% 
  mutate(mox = str_detect(modifications, "(?i)oxidation")) %>% 
  filter(mox == FALSE)

# Turnover rate profiles for hemoglobin proteins
ggplot(bfl %>% filter(master_protein_accessions %in% liv_list, tissue == "Liver"),
       aes(x = time, color = sequence)) +
  geom_line(aes(y = pred, group = id)) +
  geom_jitter(data = liv_dat, aes(y = fraction), width = 0.05) +
  facet_grid(vars(master_protein_descriptions), vars(dev_stage), scales = "free_y") +
  guides(color = FALSE) +
  theme_light(base_size = 10)


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




# Supp Fig 5C -------------------------------------------------------------


ribo_liver_lung <- ribo_prot %>% 
  filter(`Liver-E13` == 1,
         `Liver-E14` == 1,
         `Liver-E16` == 1,
         `Liver-P0` == 1,
         `Lung-E13` == 1,
         `Lung-E14` == 1,
         `Lung-E16` == 1,
         `Lung-P0` == 1)

# Boxplot of ribosomal proteins
# Detected in liver and lung
ggplot(ribosomal %>% filter(master_protein_accessions %in% ribo_liver_lung$protein, peptide_n >= 2),
       aes(x = gestation, y = estimate, color = tissue)) +
  geom_boxplot(outlier.color = NA, size = 0.25) +
  geom_point(position = position_jitterdodge(seed = 1), size = 0.4) +
  theme_bw(base_size = 10) +
  facet_wrap(~master_protein_descriptions, nrow = 5, labeller = label_wrap_gen(30)) +
  scale_y_log10()




# pdf plots ---------------------------------------------------------------

data_plot <- data %>% 
  ungroup() %>% 
  select(master_protein_accessions, master_protein_descriptions, peptide_id, 
         id, estimate) %>% 
  spread(id, estimate) %>% 
  gather(id, estimate, 4:11) %>% 
  separate(id, into = c("tissue", "dev_stage"), remove = FALSE)

prot_list <- data_plot %>% 
  select(master_protein_accessions, master_protein_descriptions) %>% 
  distinct() %>% 
  arrange(master_protein_descriptions)

nrow(prot_list)

prot_list$panel <- rep(c("A","B","C","D","E"), nrow(prot_list)/5)
prot_list$int <- rep(1:(nrow(prot_list)/5), each = 5)

data_plot <- full_join(prot_list, data_plot)


pdf(file = "tissue_turnover_rates.pdf", width = 8.5, height = 11)

for(i in 1:max(prot_list$int)){
  
  temp_data <- data_plot %>% filter(int == i)
  
  
  
  A <- ggplot(temp_data %>% filter(panel == "A")) +
    geom_boxplot(aes(x = dev_stage, y = estimate), outlier.shape = NA, size = 0.1) +
    geom_jitter(aes(x = dev_stage, y = estimate), width = 0.1, size = 0.5) +
    facet_grid(~tissue) +
    theme_light(base_size = 8) +
    labs(title = paste(temp_data$master_protein_descriptions[temp_data$panel == "A"][1]),
         subtitle = paste(temp_data$master_protein_accessions[temp_data$panel == "A"][1])) +
    scale_y_log10()
  
  B <- ggplot(temp_data %>% filter(panel == "B")) +
    geom_boxplot(aes(x = dev_stage, y = estimate), outlier.shape = NA, size = 0.1) +
    geom_jitter(aes(x = dev_stage, y = estimate), width = 0.1, size = 0.5) +
    facet_grid(~tissue) +
    theme_light(base_size = 8) +
    labs(title = paste(temp_data$master_protein_descriptions[temp_data$panel == "B"][1]),
         subtitle = paste(temp_data$master_protein_accessions[temp_data$panel == "B"][1])) +
    scale_y_log10()
  
  C <- ggplot(temp_data %>% filter(panel == "C")) +
    geom_boxplot(aes(x = dev_stage, y = estimate), outlier.shape = NA, size = 0.1) +
    geom_jitter(aes(x = dev_stage, y = estimate), width = 0.1, size = 0.5) +
    facet_grid(~tissue) +
    theme_light(base_size = 8) +
    labs(title = paste(temp_data$master_protein_descriptions[temp_data$panel == "C"][1]),
         subtitle = paste(temp_data$master_protein_accessions[temp_data$panel == "C"][1])) +
    scale_y_log10()
  
  D <- ggplot(temp_data %>% filter(panel == "D")) +
    geom_boxplot(aes(x = dev_stage, y = estimate), outlier.shape = NA, size = 0.1) +
    geom_jitter(aes(x = dev_stage, y = estimate), width = 0.1, size = 0.5) +
    facet_grid(~tissue) +
    theme_light(base_size = 8) +
    labs(title = paste(temp_data$master_protein_descriptions[temp_data$panel == "D"][1]),
         subtitle = paste(temp_data$master_protein_accessions[temp_data$panel == "D"][1])) +
    scale_y_log10()
  
  E <- ggplot(temp_data %>% filter(panel == "E")) +
    geom_boxplot(aes(x = dev_stage, y = estimate), outlier.shape = NA, size = 0.1) +
    geom_jitter(aes(x = dev_stage, y = estimate), width = 0.1, size = 0.5) +
    facet_grid(~tissue) +
    theme_light(base_size = 8) +
    labs(title = paste(temp_data$master_protein_descriptions[temp_data$panel == "E"][1]),
         subtitle = paste(temp_data$master_protein_accessions[temp_data$panel == "E"][1])) +
    scale_y_log10()
  
  
  
  print(A + B + C + D + E +
          plot_layout(ncol = 1))
  
  
}

dev.off()




