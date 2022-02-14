library(tidyverse)
library(ggsci)
library(cowplot)

setwd("/zfs/analysis/paper_shiny/")

sn_vs_vta_results <- read_csv("input/markers/sn_vs_vta_da.csv")
sn_vs_vta_counts <- read_csv("input/markers/sn_vs_vta_da_counts.csv", 
                             col_names = c("sample_cell_id", "geneID", "SCT_count"))
sn_vs_vta_metadata <- read_csv("input/markers/sn_vs_vta_da_metadata.csv")
colnames(sn_vs_vta_metadata)[1] = "sample_cell_id"

sn_vs_vta_counts %>%
  filter(geneID == "Nrip3") %>%
  # filter(SCT_count > 0) %>%
  inner_join(sn_vs_vta_metadata) %>%
  mutate(SCT_count = as.double(SCT_count)) %>%
  ggplot(aes(x = region, 
             y = SCT_count, 
             fill = region)) +
  geom_violin() +
  geom_jitter(width = 0.1) +
  scale_fill_d3() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_cowplot() +
  theme(legend.position = "none")



sn_vs_vta_counts %>%
  filter(geneID == "Calb1") %>%
  inner_join(sn_vs_vta_metadata) %>%
  mutate(SCT_count = as.double(SCT_count)) %>%
  ggplot(aes(x = x, 
             y = y, 
             fill = SCT_count)) +
  geom_point(shape = 21, 
             size = 2) +
  facet_wrap(vars(mouse_id), 
             scales = "free") +
  scale_y_reverse() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_cowplot() +
  theme(axis.text.x = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(), 
        axis.title = element_blank()) +
  panel_border()

sn_vs_vta_results