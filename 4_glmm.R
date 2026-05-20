#Title: Generalized Linear Model for predicting observed diversity
#Author: Luis A. Hurtado
#Last updated: 19-May-2026

setwd("~/Desktop/alejo")
#rm(list=ls())

#Packages####
#install.packages("tidyverse")
library(tidyverse)
#install.packages("glmmTMB")
library(glmmTMB)
#install.packages("readxl")
library(readxl)

####Load metadata and ASV counts####
# Make sure your working directory is set to where these files are located
table_s3 <- readxl::read_excel("~/Desktop/alejo/Supplementary_Tables.xlsx",
                               sheet = "Table_S3",
                               range = "Table_S3!A3:T43") %>%
  mutate(output_read_count = as.numeric(output_read_count))
table_s4 <- readxl::read_excel("~/Desktop/alejo/Supplementary_Tables.xlsx",
                               sheet = "Table_S4",
                               range = "Table_S4!A3:EL43")
table_s5 <- readxl::read_excel("~/Desktop/alejo/Supplementary_Tables.xlsx",
                               sheet = "Table_S5")
asv_table <- read.csv("taxaFinalwcounts_goodSampleNames_screenContam.csv", check.names = F)

####Family GLM####

##Derive variables##

# Extract Site name (everything before the ".")
table_s4 <- table_s4 %>%
  mutate(Site = sub("\\..*", "", samp_name))

#Calculate Number of Unique Families per Sample
asv_long <- asv_table %>%
  pivot_longer(
    cols = -c(seq_id:identificationRemarks), 
    names_to = "samp_name", 
    values_to = "read_count"
  ) %>%
  filter(!grepl("\\.0$", samp_name)) # Remove negative controls

n_families_df <- asv_long %>%
  filter(read_count > 0) %>%
  filter(!is.na(family) & family != "") %>%
  group_by(samp_name) %>%
  summarize(n_families = n_distinct(family), .groups = 'drop')

#Calculate Proportion of Streptophyta vs. Chlorophyta reads
prop_strep_df <- asv_long %>%
  filter(phylum %in% c("Streptophyta", "Chlorophyta")) %>%
  group_by(samp_name, phylum) %>%
  summarize(total_reads = sum(read_count), .groups = 'drop') %>%
  pivot_wider(names_from = phylum, values_from = total_reads, values_fill = 0) %>%
  mutate(prop_strepto = Streptophyta / (Streptophyta + Chlorophyta)) %>%
  select(samp_name, prop_strepto)

##Merge data for GLM##
meta_s3 <- table_s3 %>% select(samp_name, output_read_count)
meta_s4 <- table_s4 %>% select(samp_name, samp_size, Site)
meta_s5 <- table_s5 %>% select(samp_name, pcr_dna_vol)

final_fam_df <- meta_s4 %>%
  filter(!grepl("\\.0$", samp_name)) %>% 
  left_join(meta_s3, by = "samp_name") %>%
  left_join(meta_s5, by = "samp_name") %>%
  left_join(n_families_df, by = "samp_name") %>%
  left_join(prop_strep_df, by = "samp_name") %>%
  mutate(n_families = replace_na(n_families, 0))

final_fam_df <- drop_na(final_fam_df)

##Transform and scale predictors##
final_fam_df <- final_fam_df %>%
  mutate(
    log_read_count = log(output_read_count),
    s_log_reads    = scale(log_read_count)[,1],
    s_samp_size    = scale(samp_size)[,1],
    s_pcr_vol      = scale(pcr_dna_vol)[,1],
    s_prop_strep   = scale(prop_strepto)[,1]
  )

##Run negative binomial GLM##

model_fam_nbinom2 <- glmmTMB(
  n_families ~ s_samp_size + s_log_reads + s_pcr_vol + s_prop_strep + (1|Site),
  data = final_fam_df,
  family = nbinom2())

# Check the summary
summary(model_fam_nbinom2)

# (Optional) Run diagnose again just to be absolutely sure everything is clean!
diagnose(model_fam_nbinom2) # Diagnosis suggests use of Poisson model due to large values of nbinom2 dispersion

##Run poisson GLM##

model_fam_poisson <- glmmTMB(
  n_families ~ s_samp_size + s_log_reads + s_pcr_vol + s_prop_strep + (1|Site),
  data = final_fam_df,
  family = poisson())

summary(model_fam_poisson)

diagnose(model_fam_poisson)

#Likelihood ratio test
drop1(model_fam_poisson, test = "Chisq")

####ASV GLM####

##Derive variables##

# Calculate Number of Unique ASVs per Sample
n_asvs_df <- asv_long %>%
  filter(read_count > 0) %>%
  group_by(samp_name) %>%
  # We count the unique sequence IDs instead of families
  summarize(n_asvs = n_distinct(seq_id), .groups = 'drop')

##Merge and scale data##
meta_s3 <- table_s3 %>% select(samp_name, output_read_count)
meta_s4 <- table_s4 %>% select(samp_name, samp_size, Site)
meta_s5 <- table_s5 %>% select(samp_name, pcr_dna_vol)

final_asv_df <- meta_s4 %>%
  filter(!grepl("\\.0$", samp_name)) %>% 
  left_join(meta_s3, by = "samp_name") %>%
  left_join(meta_s5, by = "samp_name") %>%
  left_join(n_asvs_df, by = "samp_name") %>%
  left_join(prop_strep_df, by = "samp_name") %>%
  mutate(n_asvs = replace_na(n_asvs, 0))

final_asv_df <- drop_na(final_asv_df)

final_asv_df <- final_asv_df %>%
  mutate(
    log_read_count = log(output_read_count),
    s_log_reads    = scale(log_read_count)[,1],
    s_samp_size    = scale(samp_size)[,1],
    s_pcr_vol      = scale(pcr_dna_vol)[,1],
    s_prop_strep   = scale(prop_strepto)[,1]
  )

##Run negative binomial GLM##
model_asv_nbinom2 <- glmmTMB(
  n_asvs ~ s_samp_size + s_log_reads + s_pcr_vol + s_prop_strep + (1|Site),
  data = final_asv_df,
  family = nbinom2()
)

summary(model_asv_nbinom2)

diagnose(model_asv_nbinom2)

# Check for significance using robust LRTs
drop1(model_asv_nbinom2, test = "Chisq")
