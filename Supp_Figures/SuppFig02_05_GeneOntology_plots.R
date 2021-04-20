

# Libraries ---------------------------------------------------------------


library(tidyverse)
library(rio)
library(janitor)


# data import -------------------------------------------------------------


data <- import("data/20200205_Volcano_plot_Gene_Ontology_GOrilla_results.csv")


# Formatting --------------------------------------------------------------


data <- data %>% 
  clean_names()


# Plots -------------------------------------------------------------------


# Heart
# 24 hr vs 48 hr
ggplot(data %>% filter(go == "BP",
                       time == 24,
                       tissue == "H",
                       fdr_q_value < 0.01)) +
  geom_point(aes(x = reorder(description, enrichment), y = enrichment, color = fdr_q_value), size = 5, shape = 17) +
  scale_color_distiller(palette = "YlOrRd") +
  theme_bw(base_size = 14) +
  expand_limits(y = c(1,7)) +
  labs(title = "Heart",
       subtitle = "24 vs 48 hours",
       x = NULL,
       color = "q-value") +
  coord_flip()


# Heart
# 8 hr vs 48 hr
ggplot(data %>% filter(go == "BP",
                       time == 8,
                       tissue == "H",
                       fdr_q_value < 0.01)) +
  geom_point(aes(x = reorder(description, enrichment), y = enrichment, color = fdr_q_value), size = 5, shape = 17) +
  scale_color_distiller(palette = "YlOrRd") +
  theme_bw(base_size = 14) +
  expand_limits(y = c(1,7)) +
  labs(title = "Heart",
       subtitle = "8 vs 48 hours",
       x = NULL,
       color = "q-value") +
  coord_flip()


# Liver
# 24 hr vs 48 hr
ggplot(data %>% filter(go == "BP",
                       time == 24,
                       tissue == "L",
                       fdr_q_value < 0.01)) +
  geom_point(aes(x = reorder(description, enrichment), y = enrichment, color = fdr_q_value), size = 5, shape = 17) +
  scale_color_distiller(palette = "YlOrRd") +
  theme_bw(base_size = 14) +
  expand_limits(y = c(1,7)) +
  labs(title = "Liver",
       subtitle = "24 vs 48 hours",
       x = NULL,
       color = "q-value") +
  coord_flip()


# Liver
# 8 hr vs 48 hr
ggplot(data %>% filter(go == "BP",
                       time == 8,
                       tissue == "L",
                       fdr_q_value < 0.01)) +
  geom_point(aes(x = reorder(description, enrichment), y = enrichment, color = fdr_q_value), size = 5, shape = 17) +
  scale_color_distiller(palette = "YlOrRd") +
  theme_bw(base_size = 14) +
  expand_limits(y = c(1,7)) +
  labs(title = "Liver",
       subtitle = "8 vs 48 hours",
       x = NULL,
       color = "q-value") +
  coord_flip()



