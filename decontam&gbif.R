setwd("C:/Users/luisa/Desktop/bryos/dm")

#rm(list=ls())

#Title: Contamination Screening
#Author: Luis A. Hurtado
#Last updated: 29-August-2026

####Packages####

#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("biomformat", force = T) #update all worked
library(BiocManager)
#install.packages('pak')
library(pak)
#install.packages("rbiom")
library(rbiom)
#install.packages("dplyr")
library(dplyr)
#install.packages("readxl")
library(readxl)
#install.packages("tidyverse")
library(tidyverse)

####Join QIIME ASV table and taxaRaw list from MEGAN (in FAIR metadata format)####

#Convert QIIME .biom to dataframe
biom <- read_biom("C:/Users/luisa/Desktop/bryos/dm/feature-table.biom")
matrix <- biom$counts %>% as.matrix
counts <- as.data.frame(matrix)
counts_w_seqid <- tibble::rownames_to_column(counts, var = "seq_id")

#Join new QIIME dataframe with taxaRaw dataframe, which already excludes any non-Viridiplantae ASVs
taxaRaw <- readxl::read_excel("C:/Users/luisa/Desktop/bryos/dm/Supplementary_Tables.xlsx",
                              sheet = "Table_S1",
                              range = "Table_S1!A3:U1840")
taxaRawcounts <- full_join(taxaRaw, counts_w_seqid, by = "seq_id")
taxaRawcounts_noNA <- taxaRawcounts %>% drop_na(dna_sequence)

####Filter contaminants####

df <- taxaRawcounts_noNA

#Identify site negative control columns (ending in "0")
control_cols <- grep("*0$", colnames(df), value = TRUE)

#Identify sample columns (ending in "1", "2", or "3")
samp_cols <- grep("*[1-3]$", colnames(df), value=T)

#Calculate the maximum read count per taxon across all negative controls
max_ctrl <- apply(df[, control_cols, drop = FALSE], 1, max)

#Zero out sample counts if any negative control count exceeds the sample count
for (samp in samp_cols) {
  rows_to_zero <- which(max_ctrl > df[[samp]])
  df[rows_to_zero, samp] <- 0
}

#Zero out control columns entirely for downstream analysis
df[, control_cols] <- 0

#Remove ASVs with 0 cumulative reads across all sites post-filtering
df_noZero <- df[rowSums(df[, samp_cols]) > 0, ]

#Remove ASVs that did not BLAST to any kingdom
df_noZero_noKingdom <- df_noZero[!is.na(df_noZero$kingdom), ]

#Apply 96% identity and 90% query cover cutoffs
df_noZero_noKingdom_query90 <- df_noZero_noKingdom[df_noZero_noKingdom$percent_query_cover >= 90, ]
df_noZero_noKingdom_query90_id96 <- df_noZero_noKingdom_query90[df_noZero_noKingdom_query90$percent_match >= 96, ]

write.csv(df_noZero_noKingdom_query90_id96, file = "df_noZero_noKingdom_query90_id96.csv")

####Keep family level assignments or lower####

df_noZero_noKingdom_query90_id96_famOnly <- df_noZero_noKingdom_query90_id96[!is.na(df_noZero_noKingdom_query90_id96$family), ]

write.csv(df_noZero_noKingdom_query90_id96_famOnly, file = "df_noZero_noKingdom_query90_id96_famOnly.csv")

#Count total number of families
n_distinct(df_noZero_noKingdom_query90_id96_famOnly$family)

#Title: Distance between taxa at each eDNA site and the nearest corresponding occurrence record
#Author: Luis A. Hurtado
#Last updated: 29-August-2026

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
#install.packages("data.table")
library(data.table)

####Load eDNA detection data and metadata####
taxa_counts <- df_noZero_noKingdom_query90_id96_famOnly

sample_metadata <- read_excel("Supplementary_Tables.xlsx", sheet = "Table_S3", range = "Table_S3!A3:EL43") %>% 
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
#sample_cols <- sample_metadata$samp_name

edna_long <- taxa_counts %>%
  mutate(lowest_taxon = coalesce(scientificName, genus, family)) %>%
  filter(!is.na(lowest_taxon)) %>%
  select(lowest_taxon, all_of(intersect(names(taxa_counts), samp_cols))) %>%
  pivot_longer(
    cols = -lowest_taxon, 
    names_to = "samp_name", 
    values_to = "read_count"
  ) %>%
  filter(read_count > 0)

sample_metadata <- sample_metadata %>%
  mutate(samp_name = str_remove_all(samp_name, "\\."))
sample_metadata <- sample_metadata[!endsWith(sample_metadata$samp_name, "0"), ]
sample_metadata <- sample_metadata[sample_metadata$samp_name !="2ND260C2", ]

edna_sites <- edna_long %>%
  full_join(sample_metadata %>% select(samp_name, decimalLongitude, decimalLatitude), by = "samp_name") %>% 
  distinct(lowest_taxon, samp_name, decimalLongitude, decimalLatitude, .keep_all = TRUE)

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

write_csv(final_distances, "taxon_gbif_min_distances_lowest_taxa.csv")
