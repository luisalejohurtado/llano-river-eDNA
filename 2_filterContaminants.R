#Title: Screening and filtering contaminants
#Author: Luis A. Hurtado
#Last updated: 19-May-2026

setwd("~/Desktop/alejo")
#rm(list=ls())

#Load taxaFinalwcounts_goodSampleNames.csv (from first script)
all <- read.csv("taxaFinalwcounts_goodSampleNames.csv", row.names = 1)

#Use for-loop to set read count=0 for samples that have fewer reads than their corresponding negative control (for each taxon) - modified from Gemini output
# 1. Load the data 
# (row.names = 1 ensures we don't get that extra column of numbers)
df <- read.csv("taxaFinalwcounts_goodSampleNames.csv", row.names = 1, check.names = FALSE)

# 2. Automatically find all site negative control columns (ending in ".0")
control_cols <- grep("\\.0$", colnames(df), value = TRUE)

# 3. Extract the base site names (e.g., "TT.0" becomes "TT")
sites <- gsub("\\.0$", "", control_cols)

# 4. Loop through each site and apply the pairwise correction
for (site in sites) {
  
  # Identify the exact control column name for this site
  ctrl_col <- paste0(site, ".0")
  
  # Identify the corresponding sample columns (e.g., ending in .1, .2, .3)
  samp_cols <- grep(paste0("^", site, "\\.[1-9]$"), colnames(df), value = TRUE)
  
  # Loop through each individual sample replicate
  for (samp in samp_cols) {
    
    # Identify rows where the control count is strictly greater than THIS sample's count
    rows_to_zero <- which(df[[ctrl_col]] > df[[samp]])
    
    # Zero out the counts for those specific cells in the sample column
    df[rows_to_zero, samp] <- 0
  }
  
  # 5. Zero out the control column entirely
  # This is standard practice so that contamination counts don't artificially inflate your total abundances in downstream analyses.
  df[, ctrl_col] <- 0
}

# 6. Save the filtered dataset to a new CSV file
#write.csv(df, "taxaFinalwcounts_goodSampleNames_screenContam.csv", row.names = FALSE)

