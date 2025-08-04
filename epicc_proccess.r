```{r}
# Extract bins data

# List all files ending with "_bins"
file_list <- list.files(pattern = "_bins")

# Read the first file to get bin coordinates
first_df <- read.table(file_list[1], header = TRUE)
combined_df <- first_df[, 1:4]  # Keep feature, chromosome, start, end

# Loop through all files and bind their last column
for (file in file_list) {
  df <- read.table(file, header = TRUE)
  
  # Get the name of the last column (original column name)
  last_col_name <- colnames(df)[ncol(df)]
  
  # Extract the last column and preserve its name
  combined_df[[last_col_name]] <- df[[ncol(df)]]
}

# Reorder columns: chromosome, start, end, feature, then value columns
reordered_cols <- c("chromosome", "start", "end", "feature",
                    setdiff(colnames(combined_df), c("chromosome", "start", "end", "feature")))
combined_df <- combined_df[, reordered_cols]

# View result
head(combined_df)

# Optionally, save it
write.table(combined_df, "colorectal_bins.txt", sep = "\t", row.names = FALSE, quote = FALSE)
```

```{r}
# Extract group, samplename, purity info

name_list <- list.files(pattern = "_bins")
purity_list <- list.files(pattern = '_metrics')
purity_list <- purity_list[!grepl("_multiregion_metrics", purity_list)]

result_df <- data.frame()

for (i in seq_along(name_list)) {
  df_name <- read.table(name_list[i], header = TRUE)
  df_purity <- read.table(purity_list[i], header = TRUE)
  
  result_df <- rbind(result_df,data.frame(
    # Get the name of the last column (original column name)
    samplenames = colnames(df_name)[ncol(df_name)],
    
    # Get the group name
    patient_id =  df_purity$Patient,
    
    # Get the purity
    purities = df_purity$Purity,
    
    # Get the ploidy
    ploidies = df_purity$psit
  ))
  
write.csv(result_df, "EPICC_qdna_ref.csv", row.names = FALSE)
}
```

