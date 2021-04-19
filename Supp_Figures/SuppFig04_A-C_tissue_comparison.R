
# setup -------------------------------------------------------------------

library(tidyverse)
library(rio)
library(janitor)
library(broom)
library(ggrepel)


# data import -------------------------------------------------------------

model_tidy <- import("../Figure3/intermediate_data/05_out_model_nls_peptide_fraction_broom_tidy.csv")
proteins <- import("../Figure3/data/20171221_Baeza_DevAtlas_Brain_Liver_timecourse_Proteins.txt")
input_files <- import("../Figure3/data/20171221_Baeza_DevAtlas_Brain_Liver_timecourse_InputFiles.txt")


# Formatting --------------------------------------------------------------


# Selecting peptides for protein quant
model_tidy <- model_tidy %>% 
  filter(p.value < 0.05,
         std.error < 0.025,
         n_runs >= 10) %>%
  select(tissue, sequence, modifications, theo_mhplus_in_da,
         master_protein_accessions, master_protein_descriptions,
         n_runs, n_tp, term, estimate, p.value) %>% 
  distinct() %>% 
  mutate(t_half = log(2)/estimate)


# Cleaning input file names
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

# Formatting protein names
prot_long <- proteins %>% 
  clean_names() %>% 
  select(-contains("found_"),
         -contains("abundances_"),
         -contains("abundance_r"),
         -contains("search_engine"),
         -contains("number_of"),
         -contains("sequence"),
         -contains("pathway")) %>% 
  gather(temp, value, contains("abundance")) %>% 
  separate(temp, into = c("temp1", "file_id", "isotope", "temp2"), sep = "_") %>% 
  select(-contains("temp")) %>% 
  left_join(., input_files) %>% 
  spread(isotope, value) %>% 
  mutate(abundance = ifelse(!is.na(heavy) & !is.na(light), heavy + light, light)) %>% 
  filter(!is.na(heavy) | !is.na(light))

# Summarizing protein abundance 
# taking the median peptide abundance value as a proxy for protein abundance
prot_long <- prot_long %>% 
  group_by(accession, tissue, mw_in_k_da, 
           coverage_in_percent) %>% 
  summarize(abundance_median = median(abundance, na.rm = TRUE)) %>% 
  ungroup() %>% 
  rename(master_protein_accessions = "accession")

# Data for turnover vs protein abundance correlation
data_cor <- model_tidy %>% 
  filter(master_protein_accessions != "")

data_cor$master_protein_accessions <- unlist(lapply(data_cor$master_protein_accessions, function(x){
  unlist(strsplit(x, split = ";"))[1]
}))


# Merging with protein
# Filtering data
data_cor <- left_join(data_cor, prot_long) %>% 
  filter(!is.na(abundance_median))


# protein level stats -----------------------------------------------------

## Functions used for ranked ANOVA
nested_stats_test <- function(data){
  lm((sign(t_half)*rank(abs(t_half))) ~ 1 + tissue,
     data = data)
}

# Counting the tissues for each peptide observed
tissue_count <- function(x){
  length(unique(x$tissue))
}


# Nesting the data and performing the ranked ANOVA
data_nest <- model_tidy %>% 
  mutate(tissue = factor(tissue, levels = c("Brain", "Liver"))) %>% 
  group_by(master_protein_accessions, master_protein_descriptions) %>% 
  nest() %>% 
  mutate(tissue_n = map_dbl(data, tissue_count)) %>% 
  filter(tissue_n == 2) %>% 
  mutate(stats_nest = map(data, nested_stats_test))

# Extracting the model coefficients
data_stats <- data_nest %>%
  select(-data) %>%
  mutate(stats_data = map(stats_nest, tidy)) %>%
  unnest(stats_data) %>%
  arrange(p.value) %>% 
  filter(term != "(Intercept)") %>%
  select(-stats_nest) %>% 
  group_by(master_protein_accessions) %>% 
  mutate(adj.p.value = p.adjust(p.value, method = "fdr")) %>% 
  ungroup()


# Plots -------------------------------------------------------------------

# Supplemental Figure 4A
ggplot(model_tidy) +
  geom_density(aes(x = estimate, color = tissue)) +
  scale_x_log10() +
  theme_bw(base_size = 14) +
  annotate("text", y = 0.6, x = 1e-3, label = paste0("Wilcox p.value = ",round(wilcox.test(estimate ~ tissue, data = model_tidy)$p.value, digits = 8)))


# Supplemental Figure 4B
# Volcano plot
ggplot(data_stats, aes(x = estimate, y = -log10(adj.p.value))) +
  geom_point() +
  geom_text_repel(aes(color = estimate > 0,label = ifelse(adj.p.value < 0.05, master_protein_descriptions, "")), force = 10, size = 3) +
  geom_hline(yintercept = -log10(0.05)) +
  theme_bw(base_size = 14) +
  labs(x = expression(Delta~half~life~(Liver~vs~Brain)))

# Supplemental Figure 4C
# Protein abundance vs protein mw
ggplot(data_cor %>% filter(p.value < 0.05,
                       term == "k"),
       aes(x = cut_number(log10(abundance_median), n = 12), y = log10(estimate))) +
  geom_boxplot(aes(color = tissue)) +
  scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Protein abundance Log10",
       y = expression(ksyn~Log10))


