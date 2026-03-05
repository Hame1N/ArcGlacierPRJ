# Load necessary packages
library(tidyverse)

# 1. Read data
env_data <- read.delim("env.txt", sep = "\t", header = TRUE, check.names = FALSE)
mags_data <- read.delim("MAGs-TPM.txt", sep = "\t", header = TRUE, check.names = FALSE)

# 2. Correctly handle column names of MAGs data
# The first column is MAG ID, other columns are sample names
mags_names <- names(mags_data)

# Extract MAG ID column name
mag_id_col <- mags_names[1]

# Extract sample column names (excluding the first column)
sample_names <- mags_names[-1]

# 3. Convert MAGs data to presence/absence matrix
mags_matrix <- as.matrix(mags_data[, -1])  # Exclude the first column (MAG ID)
rownames(mags_matrix) <- mags_data[[1]]    # Set row names to MAG ID

# Convert to presence/absence matrix (assuming TPM > 0 indicates presence)
mags_presence_matrix <- ifelse(mags_matrix > 0, 1, 0)

# Convert to data frame
mags_presence <- as.data.frame(mags_presence_matrix)

# 4. Filter MAGs: remove MAGs that appear in only one sample
# Calculate occurrence count for each MAG
occurrence_counts <- rowSums(mags_presence)

print(table(occurrence_counts))

# Filter MAGs with occurrence count greater than 1
mags_presence_filtered <- mags_presence[occurrence_counts > 1, ]

# 5. Standardize environmental data (0-1 normalization)
env_vars <- c("Evaporation", "MAP", "NSR", "MAT", 'Glacier','Temperature','Conductivity','pH')  env_standardized <- env_data

# Perform 0-1 normalization for each environmental variable
for(var in env_vars) {
  min_val <- min(env_data[[var]], na.rm = TRUE)
  max_val <- max(env_data[[var]], na.rm = TRUE)
  # Calculate standardized values
  env_standardized[[paste0(var, "_std")]] <- (env_data[[var]] - min_val) / (max_val - min_val)
}

# 6. Create a lookup table for environmental data (using SiteID as key)
env_lookup <- env_standardized %>%
  select(SiteID, ends_with("_std")) %>%
  column_to_rownames("SiteID")


# 7. Calculate niche breadth separately for each environmental parameter
calculate_niche_breadth_by_env <- function(mag_name, mag_presence_row, env_lookup_df, env_vars) {
  present_samples <- names(mag_presence_row)[mag_presence_row > 0]
  if(length(present_samples) < 2) {
    result <- list(MAG_ID = mag_name, Samples = length(present_samples))
    for(var in env_vars) {
      result[[paste0(var, "_breadth")]] <- NA
    }
    result[["Overall_breadth"]] <- NA
    return(result)
  }
  
  # Check which samples are present in the environmental data
  valid_samples <- present_samples[present_samples %in% rownames(env_lookup_df)]
  
  if(length(valid_samples) < 2) {
    result <- list(MAG_ID = mag_name, Samples = length(valid_samples))
    for(var in env_vars) {
      result[[paste0(var, "_breadth")]] <- NA
    }
    result[["Overall_breadth"]] <- NA
    return(result)
  }
  
  # Get environmental variable values for samples where this MAG is present
  env_values <- env_lookup_df[valid_samples, ]
  
  # Initialize result list
  result <- list(MAG_ID = mag_name, Samples = length(valid_samples))
  breadth_values <- numeric(length(env_vars))
  names(breadth_values) <- env_vars
  
  # Calculate niche breadth (range of standardized values) separately for each environmental variable
  for(i in seq_along(env_vars)) {
    var <- env_vars[i]
    std_var_name <- paste0(var, "_std")
    
    if(std_var_name %in% colnames(env_values)) {
      values <- env_values[[std_var_name]]
      breadth_values[i] <- max(values) - min(values)
      result[[paste0(var, "_breadth")]] <- breadth_values[i]
    } else {
      breadth_values[i] <- NA
      result[[paste0(var, "_breadth")]] <- NA
    }
  }
  
  # Calculate average niche breadth (mean of ranges across all environmental variables)
  overall_breadth <- mean(breadth_values, na.rm = TRUE)
  result[["Overall_breadth"]] <- overall_breadth
  
  # Add range values for each environmental variable
  result[["Breadth_values"]] <- breadth_values
  
  return(result)
}

# 8. Calculate niche breadth for each environmental parameter for each filtered MAG
niche_results <- list()

for(i in 1:nrow(mags_presence_filtered)) {
  mag_name <- rownames(mags_presence_filtered)[i]
  mag_row <- mags_presence_filtered[i, ]
  
  result <- calculate_niche_breadth_by_env(mag_name, mag_row, env_lookup, env_vars)
  
  niche_results[[i]] <- result
  
  # Display progress
  if(i %% 50 == 0) {
    cat(sprintf("Processed %d/%d 个MAGs\n", i, nrow(mags_presence_filtered)))
  }
}

# 9. Extract niche breadth for each environmental variable separately
#  Create results data frame
niche_breadth_df <- data.frame(
  MAG_ID = sapply(niche_results, function(x) x$MAG_ID),
  Samples = sapply(niche_results, function(x) x$Samples),
  stringsAsFactors = FALSE
)

# Add a column for each environmental variable
for(var in env_vars) {
  col_name <- paste0(var, "_breadth")
  niche_breadth_df[[col_name]] <- sapply(niche_results, function(x) x[[col_name]])
}

# Add overall niche breadth
niche_breadth_df$Overall_breadth <- sapply(niche_results, function(x) x$Overall_breadth)

# Remove NA values
niche_breadth_df_complete <- niche_breadth_df %>%
  filter(complete.cases(across(all_of(paste0(env_vars, "_breadth")))))


# 10. Add occurrence information for each MAG
occurrence_info <- data.frame(
  MAG_ID = rownames(mags_presence),
  Occurrence = rowSums(mags_presence),
  stringsAsFactors = FALSE
)

# Merge results
niche_breadth_final <- niche_breadth_df_complete %>%
  left_join(occurrence_info, by = "MAG_ID") %>%
  arrange(desc(Overall_breadth))

# 16. Save results 
write.csv(niche_breadth_final, "MAGs_niche_breadth_by_environment.csv", row.names = FALSE)








