#Title: Beta diversity vs distance
#Author: Luis A. Hurtado
#Last updated: 19-May-2026

setwd("~/Desktop/alejo")
#rm(list=ls())

####Packages####

#install.packages("tidyverse")
library(tidyverse)
#install.packages("vegan")
library(vegan)
#install.packages("geosphere")
library(geosphere)
#install.packages("patchwork")
library(patchwork)
#install.packages("readxl")
library(readxl)

####Load and process data####
s4_data <- read_excel("Supplementary_Tables.xlsx", sheet = "Table_S4", range = "Table_S4!A3:EL43")
taxa_data <- read.csv("taxaFinalwcounts_goodSampleNames_screenContam.csv", check.names = F)

ordered_sites <- unique(sub("\\..*", "", s4_data$samp_name))
site_info <- s4_data %>%
  filter(samp_category != "negative control") %>%
  mutate(site = sub("\\..*", "", samp_name),
         lat = ifelse(decimalLatitude < 0, decimalLongitude, decimalLatitude),
         lon = ifelse(decimalLatitude < 0, decimalLatitude, decimalLongitude)) %>%
  group_by(site) %>%
  summarize(lat = mean(lat, na.rm = TRUE), lon = mean(lon, na.rm = TRUE), .groups = "drop") %>%
  mutate(River = factor(case_when(site %in% ordered_sites[1:6] ~ "South Llano", 
                                  TRUE ~ "North Llano"), levels = c("South Llano", "North Llano")))

####Mantel tests####
get_dist_and_stats <- function(df, group_col, filter_phylum = NULL) {
  df_filt <- df
  if(!is.null(filter_phylum)) df_filt <- df %>% filter(phylum %in% filter_phylum)
  obs_wide <- df_filt %>%
    pivot_longer(cols = matches("\\."), names_to = "samp_name", values_to = "count") %>%
    mutate(site = sub("\\..*", "", samp_name)) %>%
    filter(!grepl("\\.0$", samp_name)) %>%
    group_by(site, !!sym(group_col)) %>%
    summarize(total_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = !!sym(group_col), values_from = total_count, values_fill = 0)
  
  run_mantel <- function(sites, metric) {
    sub_sites <- site_info %>% filter(site %in% sites)
    sub_bio <- obs_wide %>% filter(site %in% sites) %>% arrange(match(site, sub_sites$site)) %>% column_to_rownames("site")
    if(nrow(sub_bio) < 4) return(1.0) 
    d_geo <- as.dist(distm(sub_sites[,c("lon","lat")], fun=distHaversine)/1000)
    d_bio <- vegdist(sub_bio, method = ifelse(metric == "Bray_Curtis", "bray", "jaccard"), binary = (metric == "Jaccard"))
    return(mantel(d_geo, d_bio, permutations = 999)$signif)
  }
  
  plot_df <- map_dfr(unique(site_info$River), function(riv) {
    riv_sites <- site_info$site[site_info$River == riv]
    pairs <- t(combn(riv_sites, 2)); map_dfr(1:nrow(pairs), function(i) {
      s1 <- pairs[i,1]; s2 <- pairs[i,2]
      dist_km <- distHaversine(site_info[site_info$site==s1, c("lon","lat")], site_info[site_info$site==s2, c("lon","lat")]) / 1000
      bio_pair <- obs_wide %>% filter(site %in% c(s1, s2)) %>% column_to_rownames("site")
      data.frame(River = riv, Distance_km = dist_km, Bray_Curtis = as.numeric(vegdist(bio_pair, "bray")), Jaccard = as.numeric(vegdist(bio_pair, "jaccard", binary = TRUE)))
    })
  })
  list(data = plot_df, stats = list(
    Bray_Curtis_Global = run_mantel(site_info$site, "Bray_Curtis"),
    Bray_Curtis_South = run_mantel(site_info$site[site_info$River=="South Llano"], "Bray_Curtis"),
    Bray_Curtis_North = run_mantel(site_info$site[site_info$River=="North Llano"], "Bray_Curtis"),
    Jaccard_Global = run_mantel(site_info$site, "Jaccard"),
    Jaccard_South = run_mantel(site_info$site[site_info$River=="South Llano"], "Jaccard"),
    Jaccard_North = run_mantel(site_info$site[site_info$River=="North Llano"], "Jaccard")
  ))
}

data_master <- list(
  "All Families" = get_dist_and_stats(taxa_data, "family", c("Streptophyta", "Chlorophyta")),
  "Streptophyta Families" = get_dist_and_stats(taxa_data, "family", "Streptophyta"),
  "Chlorophyta Genera" = get_dist_and_stats(taxa_data, "genus", "Chlorophyta"),
  "All ASVs" = get_dist_and_stats(taxa_data, "seq_id", c("Streptophyta", "Chlorophyta"))
)

####Final plotting function####
create_panel <- function(data_idx, metric, letter, show_y_labels = FALSE, show_x = FALSE, header_text = NULL) {
  obj <- data_master[[data_idx]]
  df <- obj$data %>% rename(Value = !!sym(metric))
  p_vals <- obj$stats
  
  global_sig <- isTRUE(p_vals[[paste0(metric, "_Global")]] < 0.05)
  south_sig  <- isTRUE(p_vals[[paste0(metric, "_South")]] < 0.05)
  north_sig  <- isTRUE(p_vals[[paste0(metric, "_North")]] < 0.05)
  
  ggplot(df, aes(x = Distance_km, y = Value, fill = River)) +
    theme_bw() + 
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_line(color = "gray92"),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank()
    ) +
    # Manual Axes drawn BEFORE points
    geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
    # Regression lines
    {if(global_sig) geom_smooth(aes(group = 1), method = "lm", color = "#014D4E", linetype = "dashed", se = FALSE, linewidth = 0.8)} +
    {if(south_sig) geom_smooth(data = filter(df, River == "South Llano"), aes(color = River), method = "lm", se = FALSE, linewidth = 1)} +
    {if(north_sig) geom_smooth(data = filter(df, River == "North Llano"), aes(color = River), method = "lm", se = FALSE, linewidth = 1)} +
    # Points drawn LAST
    geom_point(size = 4, shape = 21, color = "black", stroke = 0.8) +
    annotate("text", x = Inf, y = -Inf, label = letter, hjust = 1.4, vjust = -0.5, size = 8, fontface = "bold") +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    scale_x_continuous(expand = c(0,0)) +
    # UPDATED: Custom labeling function to remove the 0.00 label
    scale_y_continuous(
      breaks = seq(0, 1, 0.25), 
      labels = function(x) ifelse(x == 0, "", sprintf("%.2f", x)),
      expand = c(0,0)
    ) +
    scale_fill_manual(values = c("South Llano" = "lightgray", "North Llano" = "gray20")) +
    scale_color_manual(values = c("South Llano" = "gray40", "North Llano" = "gray10")) +
    labs(subtitle = header_text, x = NULL, y = NULL) +
    theme(
      plot.subtitle = element_text(size = 14, face = "italic", hjust = 0.5, margin = margin(b = 10)),
      axis.text.y = element_text(size = 12, color = if(show_y_labels) "black" else "transparent"),
      axis.ticks.y = element_line(color = if(show_y_labels) "black" else "transparent"),
      axis.text.x = if(show_x) element_text(size = 12) else element_blank(),
      axis.ticks.x = if(show_x) element_line() else element_blank(),
      legend.position = "none", 
      plot.margin = margin(5, 5, 5, 5)
    )
}

####Final stitching####
y_bray <- ggplot() + annotate("text", x = 0, y = 0, label = "Bray-Curtis", angle = 90, size = 7, fontface = "bold") + theme_void()
y_jacc <- ggplot() + annotate("text", x = 0, y = 0, label = "Jaccard", angle = 90, size = 7, fontface = "bold") + theme_void()
x_lab  <- ggplot() + annotate("text", x = 0, y = 0, label = "Geographic Distance (km)", size = 7, fontface = "bold") + theme_void()

leg_plot <- ggplot(site_info, aes(x = lat, y = lon, fill = River)) + 
  geom_point(shape = 21, size = 4) + 
  scale_fill_manual(values = c("South Llano" = "lightgray", "North Llano" = "gray20")) +
  theme_minimal() + theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 14))
shared_legend <- cowplot::get_legend(leg_plot)

p_a <- create_panel(1, "Bray_Curtis", "a", T, F, "All Families")
p_b <- create_panel(2, "Bray_Curtis", "b", F, F, "Streptophyta Families")
p_c <- create_panel(3, "Bray_Curtis", "c", F, F, "Chlorophyta Genera")
p_d <- create_panel(4, "Bray_Curtis", "d", F, F, "All ASVs")
p_e <- create_panel(1, "Jaccard", "e", T, T)
p_f <- create_panel(2, "Jaccard", "f", F, T)
p_g <- create_panel(3, "Jaccard", "g", F, T)
p_h <- create_panel(4, "Jaccard", "h", F, T)

####Build and save final plot####

final_figure <- wrap_plots(
  K = y_bray, A = (p_a|p_b|p_c|p_d),
  L = y_jacc, B = (p_e|p_f|p_g|p_h),
  C = x_lab,
  D = wrap_elements(shared_legend),
  design = "
    KAAAA
    LBBBB
    #CCCC
    #DDDD
  "
) + plot_layout(widths = c(0.04, 1), heights = c(1, 1, 0.08, 0.08))


print(final_figure) #borderline significance (p = 0.047) for autocorrelation of Jaccard Global ASVs so trendline may or may not be present

#ggsave("llano_diversity_no_zero_19May2026.svg", final_figure, width = 16, height = 10, device = "svg")



