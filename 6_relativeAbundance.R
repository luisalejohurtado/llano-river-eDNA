#Title: Family diversity and relative read abundance
#Author: Luis A. Hurtado
#Last updated: 19-May-2026

setwd("~/Desktop/alejo")
#rm(list=ls())

####Packages####

#install.packages("tidyverse")
library(tidyverse)
#install.packages("scatterpie")
library(scatterpie)
#install.packages("svglite")
library(svglite) 

####Load and process data####
df <- read.csv("taxaFinalwcounts_goodSampleNames_screenContam.csv", check.names = FALSE)
s4_raw <- read_excel("Supplementary_Tables.xlsx", sheet = "Table_S4", range = "Table_S4!A3:EL43")

site_order <- s4_raw %>%
  filter(!is.na(samp_name), samp_name != "") %>%
  mutate(CleanSite = str_extract(samp_name, "^[^.]+")) %>%
  distinct(CleanSite) %>%
  pull(CleanSite)

df_long <- df %>%
  pivot_longer(
    cols = matches(paste0("^(", paste(site_order, collapse="|"), ")\\.")),
    names_to = "Sample",
    values_to = "Count"
  ) %>%
  mutate(Site = str_extract(Sample, "^[^.]+"),
         Site = factor(Site, levels = site_order)) %>% 
  filter(!is.na(Site), Count > 0)

##Retain only families that represent at least 5% of Viridiplantae reads at ANY site##
vip_families <- df_long %>%
  group_by(Site, family) %>%
  summarise(SumCount = sum(Count), .groups = 'drop_last') %>%
  mutate(Pct = SumCount / sum(SumCount)) %>%
  filter(Pct >= 0.05) %>% 
  ungroup() %>%
  distinct(family) %>% 
  pull(family) %>% sort()

##Set custom color palette##
all_legend_items <- c(vip_families, "Other")
user_colors <- c("#4E79A7", "#A0CBE8", "#F28E2B", "#FFBE7D", "#59A14F", 
                 "#8CD17D", "#B6992D", "#F1CE63", "#499894", "#86BCB6", 
                 "#E15759", "#FF9D9A", "#B07AA1", "#D4A6C8", "#9D7660", 
                 "#D7B5A6", "#D37295", "white")
master_colors <- setNames(colorRampPalette(user_colors)(length(all_legend_items)), all_legend_items)
master_colors["Other"] <- "#D3D3D3"

####Specify data by plot row####
data_r1 <- df_long %>% filter(phylum == "Streptophyta") %>% mutate(RowLabel = "Streptophyta")
data_r2 <- df_long %>% filter(phylum == "Chlorophyta") %>% mutate(RowLabel = "Chlorophyta")
data_r3 <- df_long %>% mutate(RowLabel = "Viridiplantae")

unified_df <- bind_rows(data_r1, data_r2, data_r3) %>%
  mutate(RowLabel = factor(RowLabel, levels = c("Streptophyta", "Chlorophyta", "Viridiplantae")))

div_data <- unified_df %>%
  group_by(RowLabel, Site) %>%
  summarise(DivCount = n_distinct(family), .groups = 'drop')

max_div <- max(div_data$DivCount, na.rm = TRUE)

pie_data <- unified_df %>%
  group_by(RowLabel, Site, family) %>%
  summarise(SumCount = sum(Count), .groups = 'drop_last') %>%
  mutate(FinalTaxon = ifelse(family %in% vip_families, family, "Other")) %>%
  group_by(RowLabel, Site, FinalTaxon) %>%
  summarise(SumCount = sum(SumCount), .groups = 'drop') %>%
  group_by(RowLabel, Site) %>%
  mutate(Pct = SumCount / sum(SumCount)) %>%
  pivot_wider(names_from = FinalTaxon, values_from = Pct, values_fill = 0)

####Proportional scaling & dynamic positioning####
final_pie_data <- left_join(pie_data, div_data, by = c("RowLabel", "Site")) %>%
  mutate(X = as.numeric(Site), Y = 1,
         Radius = (sqrt(DivCount) / sqrt(max_div)) * 0.45,
         # Pushing counts a bit further down (0.2 units) to avoid overlap
         CountY = 1 - Radius - 0.20) 

site_label_df <- data.frame(
  X = 1:length(site_order), 
  Y = 1.75, # Fixed height for site names
  label = site_order,
  RowLabel = factor("Streptophyta", levels = levels(unified_df$RowLabel))
)

####Build and save final plot####
p <- ggplot() +
  # COLOR = NA removes the white outlines
  geom_scatterpie(aes(x=X, y=Y, r=Radius), data=final_pie_data, 
                  cols=intersect(all_legend_items, names(final_pie_data)), 
                  color=NA) +
  geom_text(data=final_pie_data, aes(x=X, y=CountY, label=DivCount), size=5) + 
  geom_text(data=site_label_df, aes(x=X, y=Y, label=label), size=5, fontface="bold") +
  scale_fill_manual(values = master_colors, name = NULL, breaks = all_legend_items) +
  facet_grid(RowLabel ~ ., switch = "y") +
  # Expanded ylim to c(0.2, 1.8) to give the labels room
  coord_fixed(xlim = c(0.5, 10.5), ylim = c(0.2, 1.85), clip = "off") + 
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "lines"),
    panel.spacing = unit(1.5, "lines"), 
    strip.text.y.left = element_text(angle = 0, size = 16, face = "bold", 
                                     hjust = 0.5, margin = margin(r = 40)),
    plot.margin = margin(20, 20, 100, 180) 
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# View and Save
print(p)

# Save the plot as an SVG
#ggsave(
  #filename = "taxon_distances_barplot_final_colored.svg", 
  #plot = distance_plot_final, 
  #width = 20, 
  #height = 12, 
  #device = svglite::svglite # Explicitly use the svglite graphics device
#)




