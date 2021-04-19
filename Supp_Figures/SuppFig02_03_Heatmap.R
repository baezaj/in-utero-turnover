
# Setup -------------------------------------------------------------------

library(tidyverse)
library(rio)
library(ggridges)
library(gplots)
library(RColorBrewer)


# Data Import -------------------------------------------------------------

data_tidy <- import("intermediate_data/SuppFig02_Peptide_abundance_tidy.csv", setclass = "tibble")
data_sig <- import("intermediate_data/SuppFig02_Anova_Liver_Heart_8hr_vs_24hr_48hr.csv", setclass = "tibble")

# Formatting --------------------------------------------------------------


data_sig <- data_sig %>% 
  filter(adj.p.value < 0.05) %>% 
  select(sequence, theo_m_hplus_in_da, master_protein_accessions, master_protein_descriptions) %>% 
  distinct()



# Hierarchical Clustering -------------------------------------------------


temp <- right_join(data_tidy, data_sig)


data_wide <- temp %>% 
  filter(isotope == "heavy",
         biorep < 33,
         !is.na(abundance_norm)) %>% 
  select(id, sequence, modifications, master_protein_accessions, master_protein_descriptions,
         abundance_norm) %>% 
  spread(id, abundance_norm)


data_hm <- data_wide %>%
  select(sequence, modifications, master_protein_accessions, contains("heart"), contains("liver")) %>% 
  unite(sequence_id, sequence, modifications, master_protein_accessions, sep = "_") %>% 
  column_to_rownames("sequence_id")

# Remove NAs
dat <- na.omit(data_hm)
data_hm <- na.omit(data_hm)

# Scale Rows
dat.n <- scale(t(dat))
dat.tn <- t(dat.n)

# Distance
# Calculate distance between experiments in rows
d1 <- dist(dat.n, method = "euclidean", diag = FALSE, upper = FALSE)
# Calculate distance between proteins in rows
d2 <- dist(dat.tn,method = "euclidean", diag = FALSE, upper = TRUE)

# hclust
# Clustering distance between experiments using Ward linkage
h1 <- hclust(d1, method = "ward.D2", members = NULL)
# Clustering distance between proteins using Ward linkage
h2 <- hclust(d2, method = "ward.D2", members = NULL)

# Plotting dendrogram
plot(h2, cex = 0.1, hang = -1)
abline(h = 30, col = "red")
clusters <- cutree(h2, 5)
data_hm$cluster <- clusters
# table(clusters)

# Row colorings
nofclust.height <-  length(unique(as.vector(clusters)))
selcol2 <- colorRampPalette(brewer.pal(5, "Dark2"))
clustcol.height <- selcol2(nofclust.height)

# Heatmap
heatmap.2(dat.tn,
          # Dendrogram
          Colv = as.dendrogram(h1),
          Rowv = as.dendrogram(h2),
          dendrogram = "both",
          
          # Data scaling
          scale = "none",
          
          # Level trace
          trace = "none",
          
          # Colors
          col = colorRampPalette(brewer.pal(10, "RdBu"))(50),
          
          
          # Row & Column labeling
          labRow = NA,
          cexCol = 0.9,
          margins = c(10,12),
          
          # Color key & density info
          density.info = "histogram",

          # Row colors
          RowSideColors=clustcol.height[clusters]
)

# Legend
legend("right",
       legend = unique(as.vector(clusters)),
       col = clustcol.height, 
       lty = 1,             
       lwd = 5,           
       cex = 0.8,
       title = "Cluster"
)



# Exporting data ----------------------------------------------------------

data_clust <- data_hm %>% 
  rownames_to_column("id") %>% 
  select(id, cluster, everything()) %>% 
  separate(id, c("sequence", "modifications", "master_protein_accessions"), sep = "_") %>% 
  arrange(cluster)

export(data_clust, file = "stats_data/03_Out_Heatmap_cluster_data.csv")
