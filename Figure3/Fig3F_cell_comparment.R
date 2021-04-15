
# setup -------------------------------------------------------------------


library(tidyverse)
library(rio)
library(janitor)
library(broom)
library(purrr)
library(modelr)
library(ggthemes)
library(forcats)
library(viridisLite)

# data import -------------------------------------------------------------


model_tidy <- import("intermediate_data/05_out_model_nls_peptide_fraction_broom_tidy.csv")
proteins <- import("data/20171221_Baeza_DevAtlas_Brain_Liver_timecourse_Proteins.txt")
input_files <- import("data/20171221_Baeza_DevAtlas_Brain_Liver_timecourse_InputFiles.txt")

# formatting --------------------------------------------------------------


# Selecting peptides for protein quant
data_filtered <- model_tidy %>% 
  filter(p.value < 0.05,
         std.error < 0.025,
         n_runs >= 10) %>%
  select(tissue, sequence, modifications, theo_mhplus_in_da,
         master_protein_accessions, master_protein_descriptions,
         n_runs, n_tp, term, estimate, p.value) %>% 
  distinct() %>% 
  mutate(t_half = log(2)/estimate) %>% 
  arrange(t_half)


input_files <- input_files %>% 
  clean_names() %>% 
  select(file_id, file_name) %>% 
  mutate(file_id = tolower(file_id))

input_files$file_name <- unlist(lapply(input_files$file_name, function(x){
  unlist(strsplit(x, split = "\\", fixed = TRUE))[8]
}))

input_files$tissue <- unlist(lapply(input_files$file_name, function(x){
  unlist(strsplit(x, split = "_", fixed = TRUE))[4]
}))

input_files <- input_files %>% 
  select(file_id, tissue)

prot_long <- proteins %>% 
  clean_names() %>% 
  select(-contains("found_"),
         -contains("abundances_"),
         -contains("abundance_r"),
         -contains("search_engine"),
         # -contains("number_of"),
         -contains("sequence"),
         -contains("pathway")) %>% 
  gather(temp, value, contains("abundance")) %>% 
  separate(temp, into = c("temp1", "file_id", "isotope", "temp2"), sep = "_") %>% 
  select(-contains("temp")) %>% 
  left_join(., input_files) %>% 
  spread(isotope, value) %>% 
  mutate(abundance = ifelse(!is.na(heavy) & !is.na(light), heavy + light, light))

prot_long <- prot_long %>% 
  filter(!is.na(heavy) | !is.na(light))

prot_long <- prot_long %>% 
  group_by(accession, tissue, mw_in_k_da, coverage_in_percent, number_of_peptides, number_of_ps_ms,
           biological_process, molecular_function, cellular_component) %>% 
  summarize(abundance_median = median(abundance, na.rm = TRUE)) %>% 
  ungroup()

prot_long <- prot_long %>% 
  rename(master_protein_accessions = "accession")

prot_go <- prot_long %>% 
  gather(go, name, biological_process, molecular_function, cellular_component) %>% 
  filter(name != "") %>% 
  mutate(name = strsplit(name, split = ";")) %>% 
  unnest(name)


data_filtered$master_protein_accessions <- unlist(lapply(data_filtered$master_protein_accessions, function(x){
  unlist(strsplit(x, split = ";"))[1]
}))




data_filtered <- left_join(data_filtered, prot_go) %>% 
  filter(!is.na(abundance_median))

# Adding a total group which is all compartments
all <- data_filtered %>% 
  select(-name) %>% 
  distinct() %>% 
  mutate(name = "total")

data_filtered <- bind_rows(data_filtered, all)



# plots -------------------------------------------------------------------




ggplot(data_filtered %>% filter(go == "cellular_component", 
                                name != "vacuole"),
       aes(x = name, y = log10(estimate))) +
  # aes(x = fct_reorder(name, name, na.rm = TRUE), y = log10(estimate))) +
  geom_jitter(width = 0.2, alpha = 0.01) +
  geom_boxplot(alpha = 0) +
  theme_bw(base_size = 14) +
  facet_wrap(~tissue) +
  coord_flip() +
  labs(y = "synthesis rate Log10",
       x = "Cellular Compartment")





# Stats -------------------------------------------------------------------

cc_liver <- data_filtered %>% 
  filter(go == "cellular_component",
         tissue == "Liver")

cc_brain <- data_filtered %>% 
  filter(go == "cellular_component",
         tissue == "Brain")



stats_liver <- tidy(TukeyHSD(aov(log10(estimate) ~ name,
                                data = cc_liver,
                                na.action = na.exclude)))

stats_brain <- tidy(TukeyHSD(aov(log10(estimate) ~ name,
                                 data = cc_brain,
                                 na.action = na.exclude)))


liv_plot <- stats_liver %>% 
  separate(contrast, into = c("cc1", "cc2"), sep = "-") %>% 
  mutate(adj.p.value = ifelse(adj.p.value > 0.05, NA, adj.p.value))

ggplot(liv_plot) +
  geom_tile(aes(x = cc1, y = cc2, fill = -log10(adj.p.value))) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(colors = viridis(100, option = "D", direction = 1)) +
  labs(title = "Liver",
       fill = "-Log10(q.value)")
# ggsave(filename = "figures/20201027_liver_stats.pdf")

br_plot <- stats_brain %>% 
  separate(contrast, into = c("cc1", "cc2"), sep = "-") %>% 
  mutate(adj.p.value = ifelse(adj.p.value > 0.05, NA, adj.p.value))

ggplot(br_plot) +
  geom_tile(aes(x = cc1, y = cc2, fill = -log10(adj.p.value))) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_gradientn(colors = viridis(100, option = "D", direction = 1)) +
  labs(title = "Brain",
       fill = "-Log10(q.value)") 

# ggsave(filename = "figures/20201027_brain_stats.pdf")
