#Title: Adding read counts to taxaFinal (Supplementary Table S2) and changing site names for downstrean analyses
#Author: Luis A. Hurtado
#Last updated: 19-May-2026

setwd("~/Desktop/alejo")
#rm(list=ls())

####Packages####


##Finicky stuff##

#Install and load BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomformat", force = T) #update all worked
library(BiocManager)

#pak and h5lite
#install.packages('pak')
library(pak)
#pak::pkg_install('h5lite')
library(h5lite)

##Not usually finicky##

#install.packages("rbiom")
library(rbiom)
#install.packages("dplyr")
library(dplyr)
#install.packages("readxl")
library(readxl)
#install.packages("tidyverse")
library(tidyverse)

####Load and reformat QIIME ASV read count data####

#Load asv-count-per-sample data biom file
biom <- read_biom("~/Desktop/alejo/qiime/exported-feature-table/feature-table.biom")

#Extract matrix
matrix <- biom$counts %>% as.matrix

#Transform to data frame
counts <- as.data.frame(matrix)

#Transform row names into a column called seq_id
counts_w_seqid <- tibble::rownames_to_column(counts, var = "seq_id")

####Load clean MEGAN taxa and use QIIME data from above to add read counts####

#Load curated taxaFinal table
taxaFinal <- readxl::read_excel("~/Desktop/alejo/Supplementary_Tables.xlsx",
                                sheet = "Table_S2",
                                range = "Table_S2!A3:U1068")

#Join taxaFinal and counts_w_seqid
taxaFinalwcounts <- full_join(taxaFinal, counts_w_seqid, by = "seq_id")

#Remove all taxa with NA in the dna_sequence column
taxaFinalwcounts_noNA <- taxaFinalwcounts %>% drop_na(dna_sequence)

#Write csv
#write.csv(taxaFinalwcounts_noNA, "taxaFinalwcounts.csv")

####Replace sample names so that they are cleaner and match metadata files - modified from Gemini output####
# 1. Load the data
taxa_data <- read.csv("taxaFinalwcounts.csv", check.names = FALSE)
table_s3 <- readxl::read_excel("~/Desktop/alejo/Supplementary_Tables.xlsx",
                               sheet = "Table_S3",
                               range = "Table_S3!A3:T43")

# 2. Extract the target sample names from Table S3
# Make sure to remove any leading/trailing whitespace just in case
target_names <- trimws(table_s3$samp_name)

# 3. Create a matching key (the target names without the dot)
# This mimics the format currently in taxaFinalwcounts (e.g. "TT.0" becomes "TT0")
key_names <- gsub("\\.", "", target_names)

# 4. Get the current column names of your taxa dataset
current_cols <- colnames(taxa_data)

# 5. Loop through and replace the matching column names
for (i in seq_along(target_names)) {
  # Find where the column name matches the key (no-dot version)
  match_index <- which(current_cols == key_names[i])
  
  # If a match is found, replace it with the target name (with dot)
  if (length(match_index) > 0) {
    current_cols[match_index] <- target_names[i]
  }
}

# 6. Assign the corrected column names back to the dataframe
colnames(taxa_data) <- current_cols

# 7. (Optional) Save the updated dataframe to a new CSV file
#write.csv(taxa_data, "taxaFinalwcounts_goodSampleNames.csv", row.names = FALSE)






