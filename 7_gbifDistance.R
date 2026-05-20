#Title: Distance between taxa at each eDNA site and the nearest corresponding occurrence record
#Author: Luis A. Hurtado
#Last updated: 19-May-2026

setwd("~/Desktop/alejo")
#rm(list=ls())

####Packages####
#install.packages("tidyverse")
library(tidyverse)
#install.packages("sf")
library(sf)
#install.packages("data.table")
library(data.table)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("dplyr")
library(dplyr)
#install.packages("scales")
library(scales)

####Load eDNA detection data and metadata####
taxa_counts <- read.csv("taxaFinalwcounts_goodSampleNames_screenContam.csv", check.names=F)

sample_metadata <- read_excel("Supplementary_Tables.xlsx", sheet = "Table_S4", range = "Table_S4!A3:EL43") %>% 
  mutate(
    # If longitude is positive (> 0), we know lat and lon were accidentally swapped.
    # This automatically corrects the 1ST260, 2ND260, and R307 samples!
    corrected_lon = ifelse(decimalLongitude > 0, decimalLatitude, decimalLongitude),
    corrected_lat = ifelse(decimalLongitude > 0, decimalLongitude, decimalLatitude)
  ) %>%
  mutate(
    decimalLongitude = corrected_lon,
    decimalLatitude = corrected_lat
  ) %>%
  select(-corrected_lon, -corrected_lat)

####Load GBIF data and convert to simple features object####
gbif_cols_to_keep <- c("species", "genus", "family", "decimalLatitude", "decimalLongitude")
print("Loading GBIF data...")
gbif_data <- fread("gbif_data_13May2026.csv", select = gbif_cols_to_keep)

gbif_clean <- gbif_data %>%
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
  filter(
    decimalLongitude >= -107.0 & decimalLongitude <= -93.0,
    decimalLatitude >= 25.0 & decimalLatitude <= 37.0
  )

# Free up memory
rm(gbif_data)
gc()

# Convert GBIF to spatial object
gbif_sf <- st_as_sf(gbif_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

####Prepare eDNA detection data and convert to simple features object####
sample_cols <- sample_metadata$samp_name

edna_long <- taxa_counts %>%
  mutate(lowest_taxon = coalesce(scientificName, genus, family)) %>%
  filter(!is.na(lowest_taxon)) %>%
  select(lowest_taxon, all_of(intersect(names(taxa_counts), sample_cols))) %>%
  pivot_longer(
    cols = -lowest_taxon, 
    names_to = "samp_name", 
    values_to = "read_count"
  ) %>%
  filter(read_count > 0)

edna_sites <- edna_long %>%
  left_join(sample_metadata %>% select(samp_name, decimalLongitude, decimalLatitude), by = "samp_name") %>%
  distinct(lowest_taxon, samp_name, decimalLongitude, decimalLatitude) %>%
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
  filter(
    decimalLongitude >= -107.0 & decimalLongitude <= -93.0,
    decimalLatitude >= 25.0 & decimalLatitude <= 37.0
  )

edna_sf <- st_as_sf(edna_sites, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

####Calculate shortest distance between eDNA detection location and corresponding GBIF record####
unique_taxa <- unique(edna_sites$lowest_taxon)
results <- list()

for (taxon in unique_taxa) {
  
  taxon_edna <- edna_sf %>% filter(lowest_taxon == taxon)
  taxon_gbif <- gbif_sf %>% 
    filter(species == taxon | genus == taxon | family == taxon)
  
  if (nrow(taxon_gbif) > 0) {
    dist_matrix <- st_distance(taxon_edna, taxon_gbif)
    absolute_min_dist <- min(apply(dist_matrix, 1, min))
    
    results[[taxon]] <- data.frame(
      lowest_taxon = taxon,
      min_distance_meters = as.numeric(absolute_min_dist)
    )
  } else {
    results[[taxon]] <- data.frame(
      lowest_taxon = taxon,
      min_distance_meters = NA
    )
  }
}

####Combine and save results
final_distances <- bind_rows(results) %>%
  mutate(min_distance_km = min_distance_meters / 1000)

#write_csv(final_distances, "taxon_gbif_min_distances_lowest_taxa.csv")


####Create and save plot####

# 1. Turn off scientific notation
options(scipen = 999)

# 2. Get the taxonomic hierarchy, adding the clade column (Streptophyta and Chlorophyta), and create plotmath expressions
taxon_hierarchy <- taxa_counts %>%
  mutate(lowest_taxon = coalesce(scientificName, genus, family)) %>%
  # --- NEW: Added 'phylum' to the select statement ---
  select(lowest_taxon, scientificName, genus, family, specificEpithet, phylum) %>%
  distinct(lowest_taxon, .keep_all = TRUE) %>%
  mutate(
    is_species = !is.na(specificEpithet) | grepl(" ", scientificName),
    is_family = is.na(scientificName) & is.na(genus) & !is.na(family),
    
    parse_string = case_when(
      is_species ~ paste0("italic('", lowest_taxon, "')"),                 
      is_family ~ paste0("plain('", lowest_taxon, "') ~ plain('sp.')"),                     
      TRUE ~ paste0("italic('", lowest_taxon, "') ~ plain('sp.')")                     
    )
  )

# 3. Join hierarchy to distances and filter >= 5 km
base_plot_data <- final_distances %>%
  mutate(min_distance_km = as.numeric(min_distance_km)) %>%
  filter(!is.na(min_distance_km) & min_distance_km >= 5) %>%
  left_join(taxon_hierarchy, by = "lowest_taxon")

# 4. Find all the families that have subordinate genera/species in THIS plot
families_to_drop <- base_plot_data %>%
  filter(!is_family) %>% 
  pull(family) %>%       
  unique()               

# 5. Filter out redundant families, sort, and explicitly factor the Y-axis
plot_data <- base_plot_data %>%
  filter(!(is_family & lowest_taxon %in% families_to_drop)) %>%
  arrange(min_distance_km) %>%
  mutate(lowest_taxon = factor(lowest_taxon, levels = unique(lowest_taxon)))

# 6. Extract and parse the expressions for the axis
label_map <- plot_data %>% distinct(lowest_taxon, parse_string)

my_parsed_labels <- parse(text = as.character(label_map$parse_string))
names(my_parsed_labels) <- as.character(label_map$lowest_taxon)

# 7. Create the deeply customized plot
# --- NEW: Added fill = phylum to the aesthetic mapping ---
distance_plot_final <- ggplot(plot_data, aes(y = lowest_taxon, x = min_distance_km, fill = phylum)) +
  
  # Removed the manual fill="gray" from here so it listens to our new 'phylum' mapping
  geom_col(color = "black") +
  
  # --- NEW: Define custom colors for the two phyla ---
  scale_fill_manual(
    name = "", # Sets the legend title
    values = c(
      "Chlorophyta" = "#1E8E99", 
      "Streptophyta" = "#99F9FF"
    )
  ) +
  
  scale_x_continuous(
    labels = function(x) format(x, scientific = FALSE),
    breaks = scales::pretty_breaks(n = 12), 
    expand = expansion(mult = c(0, 0.05))   
  ) +
  
  scale_y_discrete(labels = my_parsed_labels) +
  
  labs(
    title = NULL,
    y = NULL, 
    x = "Distance (km)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 22, color = "black", margin = margin(r = 5)),
    
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    panel.grid.major.x = element_line(color = "darkgray", linewidth = 0.5),
    panel.grid.minor.x = element_line(color = "lightgray", linewidth = 0.25),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    
    axis.text.x = element_text(size = 20, color = "black"),
    axis.title.x = element_text(size = 32, face = "bold", margin = margin(t = 15)),
    
    # --- NEW: Format the legend so it is large and readable ---
    legend.position = "right",
    #legend.title = element_text(size = 22, face = "bold", color = "black"),
    legend.text = element_text(size = 22, color = "black")
  )

# Display the plot
print(distance_plot_final)

# Save the plot
#ggsave("taxon_distances_barplot_final_colored_blues.png", plot = distance_plot_final, width = 10, height = 12, dpi = 300)


# Save the plot as an SVG
#ggsave(
 # filename = "taxon_distances_barplot_final_colored.svg", 
  #plot = distance_plot_final, 
  #width = 20, 
  #height = 12, 
 # device = svglite::svglite # Explicitly use the svglite graphics device
#)



