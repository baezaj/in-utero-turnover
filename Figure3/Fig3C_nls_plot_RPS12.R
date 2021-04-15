
# setup -------------------------------------------------------------------


library(tidyverse)
library(rio)
library(broom)
library(modelr)


# data import -------------------------------------------------------------

data_lm <- import("intermediate_data/02_out_Peptide_modeling_input_data.csv")
model_tidy <- import("intermediate_data/05_out_model_nls_peptide_fraction_broom_tidy.csv")

# formatting --------------------------------------------------------------


# Selecting peptides for protein quant
data_filtered <- model_tidy %>% 
  filter(p.value < 0.05,
         n_runs >= 10) %>%
  select(tissue, sequence, modifications, theo_mhplus_in_da,
         master_protein_accessions, master_protein_descriptions,
         n_runs, n_tp, term, estimate, p.value) %>% 
  distinct()


data <- left_join(data_filtered, data_lm)


# peptide level data ------------------------------------------------------

# Brain
data_brain <- data %>% 
  filter(
    master_protein_accessions == "P63323",
    sequence == "ALIHDGLAR",
    tissue == "Brain"
  ) %>% 
  distinct()

# Liver
data_liver <- data %>% 
  filter(
    master_protein_accessions == "P63323",
    sequence == "ALIHDGLAR",
    tissue == "Liver"
  ) %>% 
  distinct()



# nls modeling ------------------------------------------------------------


# running model on brain data
mod_nls_brain <- nls(fraction ~ exp(-k*time)*(k*time),
               data = data_brain,
               start = list(
                 k = coef(lm(fraction ~ time, data = data_brain))[[2]]
               )
)

# running model on liver data
mod_nls_liver <- nls(fraction ~ exp(-k*time)*(k*time),
               data = data_liver,
               start = list(
                 k = coef(lm(fraction ~ time, data = data_liver))[[2]]
               )
)


# Combining data ----------------------------------------------------------

# nls model prediction
data_model_brain <- data.frame(time = seq(0, max(data_brain$time), 0.1)) %>%
  add_predictions(mod_nls_brain) %>% 
  mutate(tissue = "Brain")

# nls model predictions
data_model_liver <- data.frame(time = seq(0, max(data_liver$time), 0.1)) %>%
  add_predictions(mod_nls_liver) %>% 
  mutate(tissue = "Liver")

# combining data for plotting
data_model <- bind_rows(data_model_brain, data_model_liver)
data_tissue <- bind_rows(data_brain, data_liver)

# Plots -------------------------------------------------------------------

# Figure 3C
ggplot(data_model, aes(x = time)) +
  geom_line(aes(y = pred, color = tissue)) +
  geom_point(data = data_tissue, aes(y = fraction, color = tissue)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  theme_minimal(base_size = 14) +
  annotate("text", x = 0.5, y = 0.5, label = data_tissue$master_protein_accessions[1], size = 5) +
  annotate("text", x = 4, y = 0.1, label = paste0("Liver p = ", round(tidy(mod_nls_liver)$p.value, digits = 8))) +
  annotate("text", x = 4, y = 0.05, label = paste0("Brain p = ", round(tidy(mod_nls_brain)$p.value, digits = 6))) +
  labs(y = "fraction",
       x = "time (hours)")

