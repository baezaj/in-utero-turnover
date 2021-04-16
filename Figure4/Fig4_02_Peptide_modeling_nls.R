
# libraries ---------------------------------------------------------------


library(tidyverse)
library(rio)
library(broom)
library(modelr)
library(ggthemes)


# data import -------------------------------------------------------------


data_lm <- import("intermediate_data/02_out_Liver_Peptide_modeling_input_data.csv")


# Formatting --------------------------------------------------------------

# calculating fraction_RIA
data_lm <- data_lm %>% 
  select(-heavy, -total, -light, -has_fraction)

# Filtering peptides seen in all time points
temp <- data_lm %>% 
  filter(time > 0) %>%
  select(sequence, modifications, theo_mhplus_in_da, 
         master_protein_accessions, master_protein_descriptions) %>% 
  distinct() %>% 
  ungroup()

data_lm <- inner_join(data_lm, temp)
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
  data_nest <- data_lm %>% 
    filter(dev_stage != "Batch") %>% 
    group_by(dev_stage, tissue, sequence, modifications, theo_mhplus_in_da,
             master_protein_accessions, master_protein_descriptions) %>% 
    nest() %>% 
    mutate(n_tp = map_dbl(data, time_count)) %>% 
    filter(n_tp == 4) %>% 
    mutate(n_runs = map_dbl(data, nrow))
)

ggplot(data_nest) +
  geom_bar(aes(x = n_runs))


# Performing the models
system.time(
  data_nest <- data_nest %>% 
    filter(n_runs >= 9) %>% 
    mutate(model_pc = map(data, nls_pulse_chase))
)


# Model summary statistics ------------------------------------------------


# Tidying the data
model_tidy <- data_nest %>% 
  filter(model_pc != "NA") %>% 
  mutate(tidy_nls = map(model_pc, tidy)) %>%
  unnest(tidy_nls) %>%
  select(-data, -model_pc)


# Glancing ----------------------------------------------------------------


# Glancing the model
model_glance <- data_nest %>%
  filter(model_pc != "NA") %>% 
  mutate(glance_nls = map(model_pc, glance)) %>%
  unnest(glance_nls) %>%
  select(-data, -model_pc)


# Augmenting --------------------------------------------------------------

# Glancing the model
model_augment <- data_nest %>%
  filter(model_pc != "NA") %>% 
  mutate(augment_nls = map(model_pc, augment)) %>%
  unnest(augment_nls) %>%
  select(-data, -model_pc)


# Model Prediction - Best fit line ----------------------------------------


# Unnesting the model predictions
data_best_fit_line <- data_nest %>% 
  filter(model_pc != "NA") %>% 
  mutate(best_fit = map(model_pc, best_fit_line)) %>% 
  unnest(best_fit)%>% 
  select(-data, -model_pc)


# Data Export -------------------------------------------------------------


export(model_glance, file = "intermediate_data/03_out_model_nls_peptide_fraction_broom_glance.csv")
export(model_tidy, file = "intermediate_data/03_out_model_nls_peptide_fraction_broom_tidy.csv")
export(data_best_fit_line, file = "intermediate_data/03_out_nls_model_BestFitLine_data.csv")


