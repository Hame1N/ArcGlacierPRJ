# Load necessary packages
library(tidyverse)
library(readr)

# Read MAGs TPM data
mag_data <- read.delim("MAGs-TPM.txt", check.names = FALSE, row.names = 1)

# Read environmental data
env_data <- read.delim("env.txt", check.names = FALSE)

# Transpose MAGs data so that sample names are row names and MAGs are column names
mag_t <- as.data.frame(t(mag_data))

# Add sample name column
mag_t$SampleID <- rownames(mag_t)

# Merge environmental data
colnames(env_data)[1] <- "SampleID"  # 将第一列列名改为SampleID

#Merge data
combined_data <- left_join(mag_t, env_data, by = "SampleID")

# Extract MAG column names (excluding SampleID and environmental variable columns)
mag_cols <- setdiff(colnames(combined_data), 
                     c("SampleID", "MAT", "Glacier"))

# Remove rows that may contain missing values
combined_data_complete <- combined_data[complete.cases(combined_data[, c("MAT", "Glacier")]), ]

# Initialize results data frame
results <- data.frame(
  MAG = mag_cols,
  rho_MAT = NA,
  p_MAT = NA,
  rho_Glacier = NA,
  p_Glacier = NA,
  Category = NA
)

# Calculate Spearman correlation between each MAG and MAT/Glacier
for (i in 1:length(mag_cols)) {
  mag_name <- mag_cols[i]
  
  # Check if the column is all zeros or NA
  if (sum(combined_data_complete[[mag_name]] > 0, na.rm = TRUE) < 3) {
    # Skip if more than half of the samples are zeros
    next
  }
  
  # Correlation with MAT
  test_mat <- tryCatch(
    cor.test(combined_data_complete[[mag_name]], 
             combined_data_complete$MAT, 
             method = "spearman", 
             use = "complete.obs"),
    error = function(e) return(list(estimate = NA, p.value = NA))
  )
  
  # Correlation with Glacier
  test_glacier <- tryCatch(
    cor.test(combined_data_complete[[mag_name]], 
             combined_data_complete$Glacier, 
             method = "spearman", 
             use = "complete.obs"),
    error = function(e) return(list(estimate = NA, p.value = NA))
  )
  
  # Store results
  results$rho_MAT[i] <- test_mat$estimate
  results$p_MAT[i] <- test_mat$p.value
  results$rho_Glacier[i] <- test_glacier$estimate
  results$p_Glacier[i] <- test_glacier$p.value
}

# Remove MAGs for which correlation could not be calculated
results <- results[complete.cases(results[, c("rho_MAT", "rho_Glacier", "p_MAT", "p_Glacier")]), ]

# Classify based on significance
results$Category <- ifelse(
  results$p_MAT < 0.05 & results$rho_MAT > 0 & 
    results$p_Glacier < 0.05 & results$rho_Glacier > 0,
  "Winner",
  ifelse(
    results$p_MAT < 0.05 & results$rho_MAT < 0 & 
      results$p_Glacier < 0.05 & results$rho_Glacier < 0,
    "Loser",
    ifelse(
      results$p_MAT >= 0.05 & results$p_Glacier >= 0.05,
      "Stabilizer",
      "Unclassified"
    )
  )
)

# View classification results
table(results$Category)

# Filter MAGs for each category
winners <- results %>% filter(Category == "Winner")
losers <- results %>% filter(Category == "Loser")
stabilizers <- results %>% filter(Category == "Stabilizer")

# Save results 
write.csv(results, "Ecological response groups.csv", row.names = FALSE)



