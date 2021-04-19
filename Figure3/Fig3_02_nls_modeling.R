
# Document Setup ----------------------------------------------------------

library(tidyverse)
library(rio)
library(broom)
library(modelr)



# Data Import -------------------------------------------------------------


data_nls <- import("intermediate_data/02_out_Peptide_modeling_input_data.csv", setclass = "tibble")


# Formatting --------------------------------------------------------------


# Filtering data
data_nls <- data_nls %>% 
  filter(!is.na(fraction)) %>% 
  select(-heavy, -total, -light, 
         # -synthesized, 
         -has_fraction)

# Filtering peptides seen in all time points
temp <- data_nls %>% 
  filter(time > 0) %>%
  select(tissue, sequence, modifications, theo_mhplus_in_da, 
         master_protein_accessions, master_protein_descriptions) %>% 
  distinct() %>% 
  ungroup()

data_nls <- inner_join(data_nls, temp)
rm(temp);gc()



# Functions ---------------------------------------------------------------


# Function used to model turnover
nls_pulse_chase <- function(data, x = time, y = fraction){
  # Using a linear model to estimate the starting parameters
  start_val <- coef(lm(fraction ~ time, data = data))[[2]]
  tryCatch(
    nls(fraction ~ exp(-k*time)*(k*time),
        data = data,
        start = list(k = start_val),
        control = nls.control(maxiter = 100)),
    error = function(e) NA, 
    warning = function(w) NA
  )
  
}


# Function that creates a dataframe with model predictions
best_fit_line <- function(data){
  data.frame(time = seq(0, 6, 0.1)) %>% 
    add_predictions(data)
}


# Function to count the number of time points
time_count <- function(x){
  length(unique(x$time))
}


# Modeling ----------------------------------------------------------------


# Nesting the data
system.time(
  data_nest <- data_nls %>% 
    group_by(tissue, sequence, modifications, theo_mhplus_in_da,
             master_protein_accessions, master_protein_descriptions) %>% 
    nest() %>% 
    mutate(n_tp = map_dbl(data, time_count)) %>% 
    filter(n_tp == 5) %>% 
    mutate(n_runs = map_dbl(data, nrow))
)


# Performing the models
system.time(
  data_nest <- data_nest %>% 
    filter(n_runs >= 10) %>% 
    mutate(model_pc = map(data, nls_pulse_chase))
)


# Model summary statistics ------------------------------------------------


# Tidying the data
model_tidy <- data_nest %>% 
  filter(model_pc != "NA") %>%
  mutate(tidy_nls = map(model_pc, tidy)) %>%
  unnest(tidy_nls) %>%
  mutate(model = "pulse_chase") %>% 
  arrange(-n_runs, sequence)%>% 
  select(-data, -model_pc) %>% 
  filter(master_protein_accessions != "")


# Glancing ----------------------------------------------------------------


# Glancing the model
model_glance <- data_nest %>%
  filter(model_pc != "NA") %>%
  mutate(glance_nls = map(model_pc, glance)) %>%
  unnest(glance_nls) %>%
  mutate(model = "pulse_chase") %>% 
  select(-data, -model_pc) %>%
  arrange(-n_runs, sequence) %>% 
  filter(master_protein_accessions != "")



# Model Prediction - Best fit line ----------------------------------------


# Unnesting the model predictions
data_best_fit_line <- data_nest %>% 
  filter(model_pc != "NA") %>% 
  mutate(best_fit = map(model_pc, best_fit_line)) %>% 
  unnest(best_fit)%>% 
  select(-data, -model_pc) %>% 
  filter(master_protein_accessions != "")



data_best_fit_line <- data_best_fit_line %>% 
  unite(id, tissue, sequence, modifications, theo_mhplus_in_da, sep = "_", remove = FALSE)

ggplot(data_best_fit_line) +
  geom_line(aes(x = time, y = pred, group = id), alpha = 0.05) +
  facet_wrap(~tissue)


# Data Export -------------------------------------------------------------


export(model_glance, file = "intermediate_data/05_out_model_nls_peptide_fraction_broom_glance.csv")
export(model_tidy, file = "intermediate_data/05_out_model_nls_peptide_fraction_broom_tidy.csv")
export(data_best_fit_line, file = "intermediate_data/05_out_nls_model_BestFitLine_data.csv")

