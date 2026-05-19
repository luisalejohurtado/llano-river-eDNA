#Title: Accumulation curves
#Author: Luis A. Hurtado
#Last updated: 19-May-2026

setwd("~/Desktop/alejo")
#rm(list=ls())

####Packages####

# Install and load readxl
#install.packages("readxl")
library(readxl)
#install.packages("tidyverse")
library(tidyverse)
#install.packages("vegan")
library(vegan)
#install.packages("svglite")
library(svglite)

####Load cleaned ASV read counts and total sample read counts####
taxa <- read.csv("taxaFinalwcounts_goodSampleNames_screenContam.csv", check.names = FALSE)

# Force output_read_count to numeric during the load to prevent errors later
table_s3 <- readxl::read_excel("~/Desktop/alejo/Supplementary_Tables.xlsx",
                               sheet = "Table_S3",
                               range = "Table_S3!A3:T43") %>%
  mutate(output_read_count = as.numeric(output_read_count))

# Extract total reads per sample from Table S3 as a lookup dictionary
total_reads_map <- setNames(table_s3$output_read_count, table_s3$samp_name)

# Identify overlapping sample columns securely
sample_cols <- intersect(colnames(taxa), table_s3$samp_name)

####Specify site colors as ordered in Table S3####
ordered_sites_from_s3 <- unique(gsub("\\.[0-9]+$", "", table_s3$samp_name))
unique_sites <- ordered_sites_from_s3[ordered_sites_from_s3 %in% unique(gsub("\\.[0-9]+$", "", sample_cols))]

site_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
                  "#E6AB02", "#A65628", "#F781BF", "#999999", "#00CED1")

color_mapping <- setNames(site_palette[1:length(unique_sites)], unique_sites)

####Pre-calculate accumulation curves####
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
  
  # Dummy device for rarecurve calculation
  tmp_img <- tempfile(fileext = ".png")
  png(tmp_img)
  rc_data <- rarecurve(matrix_data, step = 500)
  dev.off()
  unlink(tmp_img) 
  
  max_x <- 0
  max_y <- 0
  for(i in seq_along(rc_data)) {
    scaled_x <- attr(rc_data[[i]], "Subsample") * scale_factors[i]
    attr(rc_data[[i]], "Subsample_Scaled") <- scaled_x
    
    if(max(scaled_x, na.rm = TRUE) > max_x) max_x <- max(scaled_x, na.rm = TRUE)
    if(max(rc_data[[i]], na.rm = TRUE) > max_y) max_y <- max(rc_data[[i]], na.rm = TRUE)
  }
  
  site_names <- gsub("\\.[0-9]+$", "", sample_names_in_matrix)
  plot_colors <- color_mapping[site_names]
  
  return(list(
    rc_data = rc_data,
    max_x = max_x,
    max_y = max_y,
    plot_colors = plot_colors
  ))
}

cat("Calculating rarecurve data (this may take a few seconds)...\n")
data_a <- prep_accumulation_data(taxa, "Streptophyta", "family")
data_b <- prep_accumulation_data(taxa, "Chlorophyta", "genus")
data_c <- prep_accumulation_data(taxa, c("Streptophyta", "Chlorophyta"), "family")
data_d <- prep_accumulation_data(taxa, c("Streptophyta", "Chlorophyta"), "seq_id")

####Draw panel function####
draw_panel <- function(prep_obj, plot_title, y_label) {
  
  plot(1, type = "n", 
       xlim = c(0, prep_obj$max_x), ylim = c(0, prep_obj$max_y),
       xlab = "", ylab = y_label,
       main = "", 
       cex.axis = 1, cex.lab = 1.2)
  
  usr <- par("usr") # Get the plot boundaries
  x_range <- usr[2] - usr[1]
  y_range <- usr[4] - usr[3]
  
  # Calculate Bottom-Right position
  # usr[2] is the right edge, usr[3] is the bottom edge
  # We subtract 1.5% from the right and add 3% from the bottom for padding
  text(x = usr[2] - (x_range * 0.015), 
       y = usr[3] + (y_range * 0.03), 
       labels = plot_title, 
       adj = c(1, 0), # '1' right-aligns text, '0' bottom-aligns it
       font = 2, 
       cex = 1.8)
  
  for(i in seq_along(prep_obj$rc_data)) {
    lines(x = attr(prep_obj$rc_data[[i]], "Subsample_Scaled"), 
          y = prep_obj$rc_data[[i]], 
          col = prep_obj$plot_colors[i], lwd = 2)
  }
}

####Generate and export final figure####
#svglite("~/Desktop/Accumulation_Curves_Final_bottomRight.svg", width = 12, height = 10)

# Keep the 8.5 bottom margin for the legend/title spacing
par(mfrow = c(2, 2), oma = c(8.5, 1, 1, 1), mar = c(3, 4.5, 1.5, 1.5))

draw_panel(data_a, "a", "Families")
draw_panel(data_b, "b", "Genera")
draw_panel(data_c, "c", "Families")
draw_panel(data_d, "d", "ASVs")

# X-axis title (large, plain, close to plots)
mtext("Total Reads", side = 1, outer = TRUE, line = 0.5, cex = 1.6, font = 1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

# Legend at the very bottom
legend("bottom", 
       legend = names(color_mapping), 
       col = color_mapping, 
       lty = 1, lwd = 3, 
       bty = "n", 
       ncol = 5,       
       cex = 1.2,      
       title = "Site",
       inset = c(0, 0)) 

dev.off()



