# Title: R Analyses for Hurtado et al. 2026 'Waterborne eDNA metabarcoding
# of the upper Llano River watershed (Texas, USA) detects regionally unreported
# green plant taxa'

setwd("C:/Users/luisa/Desktop/bryos/dm")
rm(list=ls())

#1. Initial Data Processing and Contamination Screening ####

#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
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

#Join QIIME ASV table and taxaRaw list from MEGAN (in FAIR metadata format)
#Convert QIIME .biom to dataframe
biom <- read_biom("C:/Users/luisa/Desktop/bryos/dm/feature-table.biom")
matrix <- biom$counts %>% as.matrix
counts <- as.data.frame(matrix)
counts_w_seqid <- tibble::rownames_to_column(counts, var = "seq_id")
#Join new QIIME dataframe with taxaRaw dataframe (originally from MEGAN)
taxaRaw <- readxl::read_excel("C:/Users/luisa/Desktop/bryos/dm/Supplementary_Tables.xlsx",
                              sheet = "Table_S1",
                              range = "Table_S1!A3:U1840")
taxaRawcounts <- full_join(taxaRaw, counts_w_seqid, by = "seq_id")
taxaRawcounts_noNA <- taxaRawcounts %>% drop_na(dna_sequence)

#Filter contaminants
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
#Zero out negative control columns entirely for downstream analysis
df[, control_cols] <- 0
#Remove ASVs with 0 cumulative reads across all sites post-filtering
df_noZero <- df[rowSums(df[, samp_cols]) > 0, ]
#Remove ASVs that did not BLAST to any kingdom
df_noZero_noKingdom <- df_noZero[!is.na(df_noZero$kingdom), ]
#Apply 96% identity and 90% query cover cutoffs
df_noZero_noKingdom_query90 <- df_noZero_noKingdom[df_noZero_noKingdom$percent_query_cover >= 90, ]
df_noZero_noKingdom_query90_id96 <- df_noZero_noKingdom_query90[df_noZero_noKingdom_query90$percent_match >= 96, ]
write.csv(df_noZero_noKingdom_query90_id96, file = "df_noZero_noKingdom_query90_id96.csv") #this is the final ASV data table

#Count the number of unique families
df_noZero_noKingdom_query90_id96_famOnly <- df_noZero_noKingdom_query90_id96[!is.na(df_noZero_noKingdom_query90_id96$family), ]
write.csv(df_noZero_noKingdom_query90_id96_famOnly, file = "df_noZero_noKingdom_query90_id96_famOnly.csv") #this is the final family data table
n_distinct(df_noZero_noKingdom_query90_id96_famOnly$family) #should be 45

#2. Calculation of Shortest Distance Between eDNA Detection and GBIF Occurrence ####

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

#Load eDNA data
taxa_counts <- df_noZero_noKingdom_query90_id96_famOnly #from previous section
sample_metadata <- read_excel("Supplementary_Tables.xlsx", sheet = "Table_S3", range = "Table_S3!A3:EL43") %>% 
  mutate(
    #If longitude is positive (> 0), we know lat and lon were accidentally swapped.
    #This automatically corrects the 1ST260, 2ND260, and R307 samples!
    #i.e. these are just safeguards against erroneous coordinate inputs that have since been fixed
    corrected_lon = ifelse(decimalLongitude > 0, decimalLatitude, decimalLongitude),
    corrected_lat = ifelse(decimalLongitude > 0, decimalLongitude, decimalLatitude)
  ) %>%
  mutate(
    decimalLongitude = corrected_lon,
    decimalLatitude = corrected_lat
  ) %>%
  select(-corrected_lon, -corrected_lat)

#Load GBIF data and convert to simple features object#
gbif_cols_to_keep <- c("species", "genus", "family", "decimalLatitude", "decimalLongitude")
print("Loading GBIF data...")
gbif_data <- fread("gbif_data_13May2026.csv", select = gbif_cols_to_keep)
#Clean it
gbif_clean <- gbif_data %>%
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
  filter(
    decimalLongitude >= -107.0 & decimalLongitude <= -93.0,
    decimalLatitude >= 25.0 & decimalLatitude <= 37.0
  )
#Free up memory
rm(gbif_data)
gc()
#Convert GBIF data to simple features (spatial) object
gbif_sf <- st_as_sf(gbif_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

#Prepare eDNA detection data and convert it to a simple features object
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
#Metadata cleaning and prepping
sample_metadata <- sample_metadata %>%
  mutate(samp_name = str_remove_all(samp_name, "\\."))
sample_metadata <- sample_metadata[!endsWith(sample_metadata$samp_name, "0"), ]
sample_metadata <- sample_metadata[sample_metadata$samp_name !="2ND260C2", ]
#Join eDNA data with relevant metadata
edna_sites <- edna_long %>%
  full_join(sample_metadata %>% select(samp_name, decimalLongitude, decimalLatitude), by = "samp_name") %>% 
  distinct(lowest_taxon, samp_name, decimalLongitude, decimalLatitude, .keep_all = TRUE)
#Convert eDNA data to simple features (spatial) object
edna_sf <- st_as_sf(edna_sites, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

#Calculate shortest distance between eDNA detection location and corresponding GBIF record
#Set up for loop
unique_taxa <- unique(edna_sites$lowest_taxon)
results <- list()
#Run for loop
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

#Combine and save results
final_distances <- bind_rows(results) %>%
  mutate(min_distance_km = min_distance_meters / 1000)

write_csv(final_distances, "taxon_gbif_min_distances_lowest_taxa.csv")


#3. Accumulation Curves ####

#install.packages("readxl")
library(readxl)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("vegan")
library(vegan)
#install.packages("svglite")
library(svglite)

#Load cleaned ASV read counts and total sample read counts
taxa <- df_noZero_noKingdom_query90_id96
#Load Table S5 and clean sample names BEFORE creating map
table_s5 <- readxl::read_excel("C:/Users/luisa/Desktop/bryos/dm/Supplementary_Tables.xlsx",
                               sheet = "Table_S5",
                               range = "Table_S5!A3:T43") %>%
  mutate(output_read_count = as.numeric(output_read_count))
#Clean dirty names (shouldn't be an issue anymore but keeping just in case)
table_s5$samp_name <- gsub("\\.", "", table_s5$samp_name)
#Extract total reads per sample as lookup map
total_reads_map <- setNames(table_s5$output_read_count, table_s5$samp_name)
#Identify overlapping sample columns
sample_cols <- intersect(colnames(taxa), table_s5$samp_name)

#Prepare accumulation curves
ordered_sites_from_s5 <- unique(gsub("[0-9]+$", "", table_s5$samp_name))
site_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
                  "#E6AB02", "#A65628", "#F781BF", "#999999", "#00CED1")
color_mapping <- setNames(site_palette[seq_along(ordered_sites_from_s5)], ordered_sites_from_s5)
#Define function for pre-calculating accumulation curves
prep_accumulation_data <- function(taxa_data, target_phyla, rank_col) {
  
  matrix_data <- taxa_data %>% filter(phylum %in% target_phyla)
  
  if (rank_col != "seq_id") {
    matrix_data <- matrix_data %>% filter(!is.na(.data[[rank_col]]) & .data[[rank_col]] != "")
  }
  
  matrix_data <- matrix_data %>%
    group_by(.data[[rank_col]]) %>%
    summarise(across(all_of(sample_cols), sum, na.rm = TRUE)) %>%
    column_to_rownames(rank_col) %>% 
    t() %>% as.data.frame()
  
  matrix_data <- matrix_data[rowSums(matrix_data) > 0, , drop = FALSE]
  
  sample_names_in_matrix <- rownames(matrix_data)
  target_reads <- rowSums(matrix_data)
  total_reads <- total_reads_map[sample_names_in_matrix]
  scale_factors <- total_reads / target_reads
  
  #Dummy device for rarecurve calculation
  tmp_img <- tempfile(fileext = ".png")
  png(tmp_img)
  rc_data <- rarecurve(matrix_data, step = 500)
  dev.off()
  if (file.exists(tmp_img)) unlink(tmp_img) 
  
  max_x <- 0
  max_y <- 0
  for(i in seq_along(rc_data)) {
    sf <- scale_factors[i]
    if (is.na(sf) || is.infinite(sf)) sf <- 1
    
    scaled_x <- attr(rc_data[[i]], "Subsample") * sf
    attr(rc_data[[i]], "Subsample_Scaled") <- scaled_x
    
    if(max(scaled_x, na.rm = TRUE) > max_x) max_x <- max(scaled_x, na.rm = TRUE)
    if(max(rc_data[[i]], na.rm = TRUE) > max_y) max_y <- max(rc_data[[i]], na.rm = TRUE)
  }
  
  site_names <- gsub("[0-9]+$", "", sample_names_in_matrix)
  plot_colors <- color_mapping[site_names]
  
  return(list(
    rc_data = rc_data,
    max_x = max_x,
    max_y = max_y,
    plot_colors = plot_colors
  ))
}
#Run accumulation curve calculation
cat("Calculating rarecurve data (this may take a few seconds)...\n")
data_a <- prep_accumulation_data(taxa, c("Streptophyta", "Chlorophyta"), "family")
data_b <- prep_accumulation_data(taxa, c("Streptophyta", "Chlorophyta"), "seq_id")

#Prepare plot
draw_panel <- function(prep_obj, plot_title, y_label) {
  
  plot(1, type = "n", 
       xlim = c(0, prep_obj$max_x), ylim = c(0, prep_obj$max_y),
       xlab = "", ylab = y_label,
       main = "", 
       cex.axis = 1, cex.lab = 1.2)
  
  usr <- par("usr") #Get the plot boundaries (x1, x2, y1, y2)
  x_range <- usr[2] - usr[1]
  y_range <- usr[4] - usr[3]
  
  #Top-Right position calculation
  #usr[2] = right edge, usr[4] = top edge
  text(x = usr[2] - (x_range * 0.02), 
       y = usr[4] - (y_range * 0.03), 
       labels = plot_title, 
       adj = c(1, 1), #'1' right-aligns horizontally, '1' top-aligns vertically
       font = 2, 
       cex = 1.8)
  
  for(i in seq_along(prep_obj$rc_data)) {
    lines(x = attr(prep_obj$rc_data[[i]], "Subsample_Scaled"), 
          y = prep_obj$rc_data[[i]], 
          col = prep_obj$plot_colors[i], lwd = 2)
  }
}

#Generate and export final figure
library(svglite)
#Reset graphics device
graphics.off()
#Open SVG device
svglite("Accumulation_Curves_Final.svg", width = 11, height = 6.5)
#Define outer margins:
#Left space  = oma[2] (1.0) + mar[2] (4.5) = 5.5 lines total
#Right space = mar[4] (1.5) + oma[4] (4.0) = 5.5 lines total
#This guarantees the figure midpoint is exactly x = 0.5
par(mfrow = c(1, 2), oma = c(7.2, 1.0, 0.5, 4.0))
#Set up panel A ("ASVs")
par(mar = c(2.2, 4.5, 0.5, 1.5))
draw_panel(data_b, "a", "ASVs")
#Set up panel B ("Families"): full 4.5 left margin restored
par(mar = c(2.2, 4.5, 0.5, 1.5))
draw_panel(data_a, "b", "Families")
#Ensure full-figure overlays coordinate space (0.0 to 1.0)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
#Set up x-axis Title: Bold, centered at x = 0.5, y = 0.105
text(x = 0.5, y = 0.105, labels = "Total Reads", font = 1, cex = 1.3)
#Set up legend (vector/SVG optimized)
legend("bottom", 
       legend = names(color_mapping), 
       title = expression(bold("Site")), 
       col = color_mapping, 
       lty = 1, 
       lwd = 2.0, 
       bty = "n", 
       ncol = 5,       
       cex = 0.75,          
       y.intersp = 0.85,    
       x.intersp = 0.9,     
       seg.len = 0.8,       
       text.width = max(strwidth(names(color_mapping), cex = 0.75)) * 1.15,
       inset = c(0, 0.005))
# 9. Write and close SVG file
dev.off()

#4. Pi Charts ####

#Open SVG device
svg("taxonomic_rank_pie_charts.svg", width = 10, height = 5)

#Load data
df <- read.csv("df_noZero_noKingdom_query90_id96.csv")

#Prepare plot
rank_colors <- c(
  "domain"  = "#CFE4A1",
  "kingdom" = "#F69779",
  "phylum"  = "#96BBE5",
  "class"   = "#D89EC7",
  "order"   = "#FFE3A4",
  "family"  = "#9BD2B3",
  "genus"   = "#FCC99C",
  "species" = "#9CA9D5"
)
#Set 1x2 panel layout
par(mfrow = c(1, 2), mar = c(0.5, 0.5, 2.5, 0.5), oma = c(1, 1, 1, 1))

#Plot Chlorophyta
chloro_tbl <- table(df$taxonRank[df$phylum == "Chlorophyta"])
chloro_pct <- round(chloro_tbl / sum(chloro_tbl) * 100, 1)
chloro_labels <- paste0(names(chloro_tbl), "\n", chloro_pct, "%")
pie(
  chloro_tbl, 
  labels = chloro_labels, 
  col = rank_colors[tolower(names(chloro_tbl))],
  radius = 0.85,
  border = "white",
  cex = 0.9,
  main = "a",
  cex.main = 1.8    # Increases title size (default is 1.2)
)
#Plot Streptophyta
strep_tbl <- table(df$taxonRank[df$phylum == "Streptophyta"])
strep_pct <- round(strep_tbl / sum(strep_tbl) * 100, 1)
strep_labels <- paste0(names(strep_tbl), "\n", strep_pct, "%")
pie(
  strep_tbl, 
  labels = strep_labels, 
  col = rank_colors[tolower(names(strep_tbl))],
  radius = 0.85,
  border = "white",
  cex = 0.9,
  main = "b",
  cex.main = 1.8    # Increases title size (default is 1.2)
)

#Close and save the SVG file
dev.off()

#5. GLMMs ####

#install.packages("tidyverse")
library(tidyverse)
#install.packages("glmmTMB")
library(glmmTMB)
#install.packages("readxl")
library(readxl)
#install.packages("performance")
library(performance)

#Load data
table_s3 <- readxl::read_excel("Supplementary_Tables.xlsx",
                               sheet = "Table_S3",
                               range = "Table_S3!A3:EL43")
table_s3$samp_name <- gsub("\\.", "", table_s3$samp_name)
table_s4 <- readxl::read_excel("Supplementary_Tables.xlsx",
                               sheet = "Table_S4",)
table_s4$samp_name <- gsub("\\.", "", table_s4$samp_name)
table_s5 <- readxl::read_excel("Supplementary_Tables.xlsx",
                               sheet = "Table_S5",
                               range = "Table_S5!A3:T43") %>% 
  mutate(output_read_count = as.numeric(output_read_count))
table_s5$samp_name <- gsub("\\.", "", table_s5$samp_name)
asv_table <- read.csv("df_noZero_noKingdom_query90_id96.csv", check.names = F)
asv_table <- asv_table[,-1] #remove first column to avoid downstream issues

#Add site code column
table_s3 <- table_s3 %>%
  mutate(Site = str_sub(samp_name, 1, -2))

#Calculate number of families per sample
asv_long <- asv_table %>%
  pivot_longer(
    cols = -c(seq_id:identificationRemarks), 
    names_to = "samp_name", 
    values_to = "read_count"
  ) %>%
  filter(!grepl("\\.0$", samp_name)) #Remove negative controls
n_families_df <- asv_long %>%
  filter(read_count > 0) %>%
  filter(!is.na(family) & family != "") %>%
  group_by(samp_name) %>%
  summarize(n_families = n_distinct(family), .groups = 'drop')

#Standardize data for merging
meta_s3 <- table_s3 %>% select(samp_name, samp_size, Site) %>%
  filter(!grepl("*0$", samp_name)) #Remove negative controls
meta_s3 <- meta_s3 %>% filter(samp_name !='2ND260C.2') #Remove sample that was not sequenced
meta_s4 <- table_s4 %>% select(samp_name, pcr_dna_vol, pcr_cycles)%>%
  filter(!grepl("*0$", samp_name)) #Remove negative controls
meta_s4 <- meta_s4 %>% filter(samp_name !='2ND260C.2') #Remove sample that was not sequenced
meta_s5 <- table_s5 %>% select(samp_name, output_read_count)%>%
  filter(!grepl("*0$", samp_name)) #Remove negative controls
meta_s5 <- drop_na(meta_s5) #Remove sample that was not sequenced

#Merge data
final_fam_df <- meta_s3 %>%
  left_join(meta_s4, by = "samp_name") %>%
  left_join(meta_s5, by = "samp_name") %>%
  left_join(n_families_df, by = "samp_name") %>%
  mutate(n_families = replace_na(n_families, 0))
final_fam_df <- drop_na(final_fam_df)

#Transform and scale GLMM predictors
final_fam_df <- final_fam_df %>%
  mutate(
    log_read_count = log(output_read_count),
    s_log_reads    = scale(log_read_count)[,1],
    s_samp_size    = scale(samp_size)[,1],
    s_pcr_vol      = scale(pcr_dna_vol)[,1],
  )

#Run all-samples Poisson GLMM
model_fam_poisson <- glmmTMB(
  n_families ~ s_samp_size + s_log_reads + s_pcr_vol + (1|Site),
  data = final_fam_df,
  family = poisson())
#Summarize and diagnose all-samples Poisson GLMM
summary(model_fam_poisson)
diagnose(model_fam_poisson)
VarCorr(model_fam_poisson)
drop1(model_fam_poisson, test = "Chisq")
check_overdispersion(model_fam_poisson)
exp(cbind(
  IRR = fixef(model_fam_poisson)$cond,
  confint(model_fam_poisson, parm = "beta_")
))

#Set up and run Poisson GLMM excluding the four samples that underwent only one round of PCR
final_fam_df_2rd <- final_fam_df %>% filter(pcr_cycles >= 35)
model_fam_poisson_2rd <- glmmTMB(
  n_families ~ s_samp_size + s_log_reads + s_pcr_vol + (1|Site),
  data = final_fam_df_2rd,
  family = poisson())
#Summarize and diagnose two-rounds-only (2rd) Poisson GLMM
summary(model_fam_poisson_2rd)
diagnose(model_fam_poisson_2rd)
VarCorr(model_fam_poisson_2rd)
drop1(model_fam_poisson_2rd, test = "Chisq")
check_overdispersion(model_fam_poisson_2rd)
exp(cbind(
  IRR = fixef(model_fam_poisson_2rd)$cond,
  confint(model_fam_poisson_2rd, parm = "beta_")
))

#Set up ASV Poisson GLMM
n_asvs_df <- asv_long %>%
  filter(read_count > 0) %>%
  group_by(samp_name) %>%
  summarize(n_asvs = n_distinct(seq_id), .groups = 'drop')
#Merge data
final_asv_df <- meta_s3 %>%
  left_join(meta_s4, by = "samp_name") %>%
  left_join(meta_s5, by = "samp_name") %>%
  left_join(n_asvs_df, by = "samp_name") %>%
  mutate(n_asvs = replace_na(n_asvs, 0))
final_asv_df <- drop_na(final_asv_df)
#Transform and scale GLMM predictors
final_asv_df <- final_asv_df %>%
  mutate(
    log_read_count = log(output_read_count),
    s_log_reads    = scale(log_read_count)[,1],
    s_samp_size    = scale(samp_size)[,1],
    s_pcr_vol      = scale(pcr_dna_vol)[,1],
  )
#Run all-samples Poisson GLMM
model_asv_poisson <- glmmTMB(
  n_asvs ~ s_samp_size + s_log_reads + s_pcr_vol + (1|Site),
  data = final_asv_df,
  family = poisson())
#Summarize and diagnose all-samples Poisson GLMM
summary(model_asv_poisson)
diagnose(model_asv_poisson)
VarCorr(model_asv_poisson)
drop1(model_asv_poisson, test = "Chisq")
check_overdispersion(model_fam_poisson)
exp(cbind(
  IRR = fixef(model_asv_poisson)$cond,
  confint(model_asv_poisson, parm = "beta_")
))
#Set up and run Poisson GLMM excluding the four samples that underwent only one round of PCR
final_asv_df_2rd <- final_asv_df %>% filter(pcr_cycles >= 35)
model_asv_poisson_2rd <- glmmTMB(
  n_asvs ~ s_samp_size + s_log_reads + s_pcr_vol + (1|Site),
  data = final_asv_df_2rd,
  family = poisson())
#Summarize and diagnose two-rounds-only (2rd) Poisson GLMM
summary(model_asv_poisson_2rd)
diagnose(model_asv_poisson_2rd)
VarCorr(model_asv_poisson_2rd)
drop1(model_asv_poisson_2rd, test = "Chisq")
check_overdispersion(model_asv_poisson_2rd)

#Overdispersion is present so switch to nbinom2 instead of Poisson
model_asv_nb2_2rd <- glmmTMB(
  n_asvs ~ s_samp_size + s_log_reads + s_pcr_vol + (1 | Site),
  data = final_asv_df_2rd,
  family = nbinom2()
)
summary(model_asv_poisson_2rd)
diagnose(model_asv_poisson_2rd)
VarCorr(model_asv_poisson_2rd)
drop1(model_asv_nb2_2rd, test = "Chisq")
check_overdispersion(model_asv_nb2_2rd)
#Calculate incidence rate ratios (IRRs)
exp(cbind(
  IRR = fixef(model_asv_nb2_2rd)$cond,
  confint(model_asv_nb2_2rd, parm = "beta_")
))

